#pragma once

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>
#include <gtl/phmap.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "kxsort.h"


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


struct FlatMap {
    using Key = uint64_t;
    using Value = uint32_t;
    using Span = std::span<const Value>;

    Span find(Key key) const {
        auto it = index.find(key);
        if (it == index.end())
            return {};

        const auto& e = it->second;
        return {values.data() + e.offset, e.size};
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
    static FlatMap from_sorted(const std::vector<KeyValue> &pairs) {
        FlatMap result;
        size_t no_keys = 0;
        for (size_t i = 0; i < pairs.size(); i++)
            no_keys += (i == 0 || pairs[i].key != pairs[i-1].key);

        result.index.reserve(no_keys);
        result.values.reserve(pairs.size());
        for (size_t i = 0; i < pairs.size();) {
            const Key key = pairs[i].key;
            const uint32_t offset = static_cast<uint32_t>(i);
            for (; i < pairs.size() && pairs[i].key == key; i++)
                result.values.push_back(pairs[i].value);
            result.index.emplace(key, Entry{.offset = offset, .size = static_cast<uint32_t>(i - offset)});
        }

        return result;
    }

    size_t memory_bits() const {
        size_t total = 0;

        total += sizeof(index);
        total += index.capacity() * sizeof(decltype(index)::value_type);

        total += sizeof(values);
        total += values.capacity() * sizeof(Value);

        return total * 8;
    }

private:
    struct Entry {
        uint32_t offset;
        uint32_t size;
    };

    gtl::flat_hash_map<Key, Entry> index;
    std::vector<Value> values;

    template <class Archive>
    friend void save(Archive&, const FlatMap&);

    template <class Archive>
    friend void load(Archive&, FlatMap&);
};

