#pragma once

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>
#include <gtl/phmap.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>


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

    friend struct FlatMapBuilder;

    template <class Archive>
    friend void save(Archive&, const FlatMap&);

    template <class Archive>
    friend void load(Archive&, FlatMap&);
};


struct FlatMapBuilder {
    using Key = FlatMap::Key;
    using Value = FlatMap::Value;
    using Span = FlatMap::Span;

    void reserve(size_t keys) {
        map.reserve(keys);
    }

    void add(Key key, Value value) {
        map[key].push_back(value);
    }

    void add(Key key, Span values) {
        auto& v = map[key];
        v.insert(v.end(), values.begin(), values.end());
    }

    template<typename Container>
    void add(Key key, const Container& values) {
        add(key, Span(values.data(), values.size()));
    }

    FlatMap build() && {
        FlatMap result;

        // Compute total storage needed; overlapping (extended) super-k-mers can add the same value twice.
        size_t total_values = 0;
        for (auto& [_, v] : map) {
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            total_values += v.size();
        }

        result.index.reserve(map.size());
        result.values.reserve(total_values);

        for (auto& [key, v] : map) {
            FlatMap::Entry e{
                .offset = static_cast<uint32_t>(result.values.size()),
                .size   = static_cast<uint32_t>(v.size())
            };

            result.values.insert(
                result.values.end(),
                std::make_move_iterator(v.begin()),
                std::make_move_iterator(v.end()));

            result.index.emplace(key, e);
        }

        return result;
    }

private:
    gtl::flat_hash_map<Key, std::vector<Value>> map;
};
