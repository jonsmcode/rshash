#pragma once

#include <algorithm>
#include <bit>
#include <cstdint>
#include <span>
#include <vector>
#include <gtl/phmap.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "kxsort.h"
#include "compact_vector.hpp"


// (key, value) pair of a structure built by sorting instead of hashing every pair
struct [[gnu::packed]] KeyValue {
    uint64_t key;
    uint32_t value;
};

struct RadixTraitsKeyValue {
    static const int nBytes = 12;
    int kth_byte(const KeyValue &x, int k) {
        if (k >= 4) return x.key >> ((k - 4) * 8) & 0xFF;
        return x.value >> (k * 8) & 0xFF;
    }
    bool compare(const KeyValue &x, const KeyValue &y) {
        return x.key < y.key || (x.key == y.key && x.value < y.value);
    }
};

// sorts the pairs by key and value and removes duplicates (overlapping super-k-mers can add a pair twice)
inline void sort_unique(std::vector<KeyValue> &pairs) {
    kx::radix_sort(pairs.begin(), pairs.end(), RadixTraitsKeyValue());
    pairs.erase(std::unique(pairs.begin(), pairs.end(),
        [](const KeyValue &x, const KeyValue &y) { return x.key == y.key && x.value == y.value; }), pairs.end());
}


// the values (text positions) are stored in a compact vector with the bit width of the largest value; compiled with
// FLATMAP_STD_VECTOR, they are stored as 32-bit integers in a std::vector (for comparison); compiled with
// FLATMAP_ENCODED, the compact vector stores the positional information as the u2o table of piscem-style rshash:
// (sequence id << (ref_len_bits + 1)) | (position in sequence << 1) | forward, see from_sorted_encoded
struct FlatMap {
    using Key = uint64_t;
    using Value = uint32_t;
#ifdef FLATMAP_STD_VECTOR
    using Values = std::vector<Value>;
    using Span = std::span<const Value>;
#else
    using Values = bits::compact_vector;

    // the values of one key: values[offset, offset + size), each read with random access() (one unaligned load,
    // no loop-carried state), which is faster than the sequential enumerator of the compact vector
    struct Span {
        struct iterator {
            bits::compact_vector const* values_;
            uint64_t i_;

            uint64_t operator*() const { return values_->access(i_); }
            iterator& operator++() { ++i_; return *this; }
            bool operator==(iterator const& other) const { return i_ == other.i_; }
            bool operator!=(iterator const& other) const { return i_ != other.i_; }
        };

        bits::compact_vector const* values_ = nullptr;
        uint64_t start_ = 0;
        size_t len_ = 0;

        bool empty() const { return len_ == 0; }
        size_t size() const { return len_; }
        iterator begin() const { return {values_, start_}; }
        iterator end() const { return {values_, start_ + len_}; }
    };
#endif

    Span find(Key key) const {
        auto it = index.find(key);
        if (it == index.end()) {
#ifdef FLATMAP_STD_VECTOR
            return {};
#else
            return {};
#endif
        }

        const auto& e = it->second;
#ifdef FLATMAP_STD_VECTOR
        return {values.data() + e.offset, e.size};
#else
        return {&values, e.offset, e.size};
#endif
    }

    bool contains(Key key) const {
        return index.contains(key);
    }

    Span operator[](Key key) const {
        return find(key);
    }

    size_t size() const {
        return index.size();
    }

    // from pairs sorted by sort_unique: every key is hashed once, the values are copied in order
#ifdef FLATMAP_ENCODED
    // number of bits of the position in a sequence (the u2o table's ref_len_bits), needed to decode an entry
    uint64_t ref_len_bits = 0;

    // from pairs sorted by sort_unique whose keys are shifted left by key_shift bits (the shifted out bits are
    // used by encode, e.g. the orientation); encode(pair) returns the entry of width bits
    template <typename Encode>
    static FlatMap from_sorted_encoded(const std::vector<KeyValue> &pairs, const uint64_t width, const int key_shift,
        const uint64_t ref_len_bits, Encode &&encode) {
        FlatMap result;
        result.ref_len_bits = ref_len_bits;
        size_t no_keys = 0;
        for (size_t i = 0; i < pairs.size(); i++)
            no_keys += (i == 0 || (pairs[i].key >> key_shift) != (pairs[i-1].key >> key_shift));

        result.index.reserve(no_keys);
        bits::compact_vector::builder builder(pairs.size(), width);
        for (size_t i = 0; i < pairs.size();) {
            const Key key = pairs[i].key >> key_shift;
            const uint32_t offset = static_cast<uint32_t>(i);
            for (; i < pairs.size() && (pairs[i].key >> key_shift) == key; i++)
                builder.push_back(encode(pairs[i]));
            result.index.emplace(key, Entry{.offset = offset, .size = static_cast<uint32_t>(i - offset)});
        }
        builder.build(result.values);

        return result;
    }
#endif

    static FlatMap from_sorted(const std::vector<KeyValue> &pairs) {
        FlatMap result;
        size_t no_keys = 0;
        for (size_t i = 0; i < pairs.size(); i++)
            no_keys += (i == 0 || pairs[i].key != pairs[i-1].key);

        result.index.reserve(no_keys);
#ifdef FLATMAP_STD_VECTOR
        result.values.reserve(pairs.size());
#else
        Value max_value = 0;
        for (const KeyValue &pair : pairs)
            max_value = std::max<Value>(max_value, pair.value);
        bits::compact_vector::builder builder(pairs.size(), std::max<uint64_t>(1, std::bit_width(max_value)));
#endif
        for (size_t i = 0; i < pairs.size();) {
            const Key key = pairs[i].key;
            const uint32_t offset = static_cast<uint32_t>(i);
            for (; i < pairs.size() && pairs[i].key == key; i++) {
#ifdef FLATMAP_STD_VECTOR
                result.values.push_back(pairs[i].value);
#else
                builder.push_back(pairs[i].value);
#endif
            }
            result.index.emplace(key, Entry{.offset = offset, .size = static_cast<uint32_t>(i - offset)});
        }
#ifndef FLATMAP_STD_VECTOR
        builder.build(result.values);
#endif

        return result;
    }

    size_t memory_bits() const {
        size_t total = 0;

        total += sizeof(index);
        total += index.capacity() * sizeof(decltype(index)::value_type);

#ifdef FLATMAP_STD_VECTOR
        total += sizeof(values);
        total += values.capacity() * sizeof(Value);
#else
        total += values.bytes();
#endif

        return total * 8;
    }

private:
    struct Entry {
        uint32_t offset;
        uint32_t size;
    };

    gtl::flat_hash_map<Key, Entry> index;
    Values values;

    template <class Archive>
    friend void save(Archive&, const FlatMap&);

    template <class Archive>
    friend void load(Archive&, FlatMap&);
};

