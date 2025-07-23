// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>


static uint64_t MurmurHash2_64(void const* key, size_t len, uint64_t seed) {
    const uint64_t m = 0xc6a4a7935bd1e995ULL;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

#if defined(__arm) || defined(__arm__)
    const size_t ksize = sizeof(uint64_t);
    const unsigned char* data = (const unsigned char*)key;
    const unsigned char* end = data + (std::size_t)(len / 8) * ksize;
#else
    const uint64_t* data = (const uint64_t*)key;
    const uint64_t* end = data + (len / 8);
#endif

    while (data != end) {
#if defined(__arm) || defined(__arm__)
        uint64_t k;
        memcpy(&k, data, ksize);
        data += ksize;
#else
        uint64_t k = *data++;
#endif

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char* data2 = (const unsigned char*)data;

    switch (len & 7) {
        // fall through
        case 7:
            h ^= uint64_t(data2[6]) << 48;
        // fall through
        case 6:
            h ^= uint64_t(data2[5]) << 40;
        // fall through
        case 5:
            h ^= uint64_t(data2[4]) << 32;
        // fall through
        case 4:
            h ^= uint64_t(data2[3]) << 24;
        // fall through
        case 3:
            h ^= uint64_t(data2[2]) << 16;
        // fall through
        case 2:
            h ^= uint64_t(data2[1]) << 8;
        // fall through
        case 1:
            h ^= uint64_t(data2[0]);
            h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

struct murmurhash2_64 {
    // generic range of bytes
    // static inline uint64_t hash(byte_range range, uint64_t seed) {
    //     return MurmurHash2_64(range.begin, range.end - range.begin, seed);
    // }

    // specialization for std::string
    static inline uint64_t hash(std::string const& val, uint64_t seed) {
        return MurmurHash2_64(val.data(), val.size(), seed);
    }

    // specialization for uint64_t
    static inline uint64_t hash(uint64_t val, uint64_t seed) {
        return MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
    }
};


namespace srindex
{

struct two_minimisers_and_window_hash_parameters
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed1{};
    uint64_t seed2{};
};

struct two_minimisers_and_window_hash_result
{
    uint64_t minimiser1_value;
    uint64_t minimiser2_value;
    uint64_t window_value;
    uint64_t window_value_rev;
};

}

namespace srindex::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class two_minimisers_and_window_hash : public std::ranges::view_interface<two_minimisers_and_window_hash<range_t>>
{
private:
    range_t range{};
    two_minimisers_and_window_hash_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    two_minimisers_and_window_hash()
        requires std::default_initializable<range_t>
    = default;
    two_minimisers_and_window_hash(two_minimisers_and_window_hash const & rhs) = default;
    two_minimisers_and_window_hash(two_minimisers_and_window_hash && rhs) = default;
    two_minimisers_and_window_hash & operator=(two_minimisers_and_window_hash const & rhs) = default;
    two_minimisers_and_window_hash & operator=(two_minimisers_and_window_hash && rhs) = default;
    ~two_minimisers_and_window_hash() = default;

    explicit two_minimisers_and_window_hash(range_t range, two_minimisers_and_window_hash_parameters params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class two_minimisers_and_window_hash<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = two_minimisers_and_window_hash_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t window_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t seed1{};
    uint64_t seed2{};
    uint8_t minimisers_in_window{};
    uint64_t minimiser_size{};
    uint64_t window_size{};

    size_t minimiser1_position{};
    size_t minimiser2_position{};

    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};

    std::deque<uint64_t> kmer_values1_in_window{};
    std::deque<uint64_t> kmer_values2_in_window{};

    static inline constexpr uint64_t compute_mask(uint64_t const size)
    {
        assert(size > 0u);
        assert(size <= 64u);

        if(size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        window_mask{it.window_mask},
        kmer_value{it.kmer_value},
        kmer_value_rev{it.kmer_value_rev},
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current},
        kmer_values1_in_window{it.kmer_values1_in_window},
        kmer_values2_in_window{it.kmer_values2_in_window},
        minimiser1_position{it.minimiser1_position},
        minimiser2_position{it.minimiser2_position}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   two_minimisers_and_window_hash_parameters const & params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * params.minimiser_size)},
        window_mask{compute_mask(2u * params.window_size)},
        range_size{range_size}
    {
        if (range_size < params.window_size)
            range_position = range_size;
        else
            init(params);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position == lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        while (!next_minimiser_is_new())
        {}
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        while (!next_minimiser_is_new())
        {}
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return current;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
    {
        uint64_t const new_rank = seqan3::to_rank(*range_it);

        current.window_value = ((current.window_value << 2) | new_rank) & window_mask;
        // kmer_value = current.window_value & kmer_mask;
        kmer_value = ((kmer_value << 2) | new_rank) & kmer_mask;
        current.window_value_rev = (current.window_value_rev >> 2) | ((new_rank^0b11) << 2*(window_size-1));
        // kmer_value_rev = current.window_value_rev >> 2*(window_size - minimiser_size);
        kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes) {
            kmer_values1_in_window.pop_front();
            kmer_values2_in_window.pop_front();
        }

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash2 = murmurhash2_64::hash(canonical_kmer, seed2) & kmer_mask;
        
        kmer_values1_in_window.push_back(kmerhash1);
        kmer_values2_in_window.push_back(kmerhash2);
    }

    void find_minimiser1_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values1_in_window, std::less_equal<uint64_t>{});
        current.minimiser1_value = *minimiser_it;
        minimiser1_position = std::distance(std::begin(kmer_values1_in_window), minimiser_it);
    }

    void find_minimiser2_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values2_in_window, std::less_equal<uint64_t>{});
        current.minimiser2_value = *minimiser_it;
        minimiser2_position = std::distance(std::begin(kmer_values2_in_window), minimiser_it);
    }

    void init(two_minimisers_and_window_hash_parameters const & params)
    {
        minimiser_size = params.minimiser_size;
        window_size = params.window_size;
        minimisers_in_window = window_size - minimiser_size;
        seed1 = params.seed1 & kmer_mask;
        seed2 = params.seed2;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        current.window_value <<= 2;
        current.window_value |= new_rank;
        current.window_value_rev >>= 2;
        current.window_value_rev |= ((new_rank^0b11) << 2*(window_size-1));
        kmer_value = ((kmer_value << 2) | new_rank) & kmer_mask;
        kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        for (size_t i = 1u; i < params.minimiser_size; ++i) {
            ++range_position;
            ++range_it;
            new_rank = seqan3::to_rank(*range_it);
            current.window_value <<= 2;
            current.window_value |= new_rank;
            current.window_value_rev >>= 2;
            current.window_value_rev |= ((new_rank^0b11) << 2*(window_size-1));
            kmer_value = ((kmer_value << 2) | new_rank) & kmer_mask;
            kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        }
        // kmer_value = current.window_value & kmer_mask;
        // kmer_value_rev = current.window_value_rev >> 2*(window_size - minimiser_size);

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash2 = murmurhash2_64::hash(canonical_kmer, seed2) & kmer_mask;

        kmer_values1_in_window.push_back(kmerhash1);
        kmer_values2_in_window.push_back(kmerhash2);

        for (size_t i = params.minimiser_size; i < params.window_size; ++i)
            next_window<pop_first::no>();

        find_minimiser1_in_window();
        find_minimiser2_in_window();
    }

    bool next_minimiser_is_new()
    {
        // assert(minimiser1_position != minimiser2_position);

        // If we reached the end of the range, we are done.
        if (range_position + 1 == range_size)
            return ++range_position; // Return true, but also increment range_position

        next_window<pop_first::yes>();

        bool new_min1 = false;
        bool new_min2 = false;

        // The minimiser left the window.
        if (minimiser1_position == 0) {
            find_minimiser1_in_window();
            new_min1 = true;
            // return true;
        }

        if (minimiser2_position == 0) {
            find_minimiser2_in_window();
            new_min2 = true;
            // return true;
        }

        if (uint64_t new_kmer_value = kmer_values1_in_window.back(); new_kmer_value < current.minimiser1_value) {
            current.minimiser1_value = new_kmer_value;
            minimiser1_position = minimisers_in_window;
            new_min1 = true;
            // return true;
        }

        if (uint64_t new_kmer_value = kmer_values2_in_window.back(); new_kmer_value < current.minimiser2_value) {
            current.minimiser2_value = new_kmer_value;
            minimiser2_position = minimisers_in_window;
            new_min2 = true;
            // return true;
        }

        if(!new_min1)
            --minimiser1_position;
        if(!new_min2)
            --minimiser2_position;

        return true;
    }
};


template <std::ranges::viewable_range rng_t>
two_minimisers_and_window_hash(rng_t &&, two_minimisers_and_window_hash_parameters &&)
    -> two_minimisers_and_window_hash<std::views::all_t<rng_t>>;

struct two_minimisers_and_window_hash_fn
{
    constexpr auto operator()(two_minimisers_and_window_hash_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, two_minimisers_and_window_hash_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.minimiser_size == 0u)
            throw std::invalid_argument{"minimiser_size must be > 0."};
        if (params.minimiser_size > 32u)
            throw std::invalid_argument{"minimiser_size must be <= 32."};
        if (params.window_size == 0u)
            throw std::invalid_argument{"window_size must be > 0."};
        if (params.window_size > 32u)
            throw std::invalid_argument{"window_size must be <= 32."};
        if (params.window_size < params.minimiser_size)
            throw std::invalid_argument{"window_size must be >= minimiser_size."};

        return two_minimisers_and_window_hash{range, std::move(params)};
    }
};

}

namespace srindex::views
{

inline constexpr auto two_minimisers_and_window_hash = srindex::detail::two_minimisers_and_window_hash_fn{};

}


namespace srindex
{

struct xor_minimiser_and_positions_parameters
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed{};
};

struct xor_minimiser_and_positions_result
{
    uint64_t minimiser_value;
    size_t range_position;
    size_t occurrences;
};

}

namespace srindex::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class xor_minimiser_and_positions : public std::ranges::view_interface<xor_minimiser_and_positions<range_t>>
{
private:
    range_t range{};
    xor_minimiser_and_positions_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    xor_minimiser_and_positions()
        requires std::default_initializable<range_t>
    = default;
    xor_minimiser_and_positions(xor_minimiser_and_positions const & rhs) = default;
    xor_minimiser_and_positions(xor_minimiser_and_positions && rhs) = default;
    xor_minimiser_and_positions & operator=(xor_minimiser_and_positions const & rhs) = default;
    xor_minimiser_and_positions & operator=(xor_minimiser_and_positions && rhs) = default;
    ~xor_minimiser_and_positions() = default;

    explicit xor_minimiser_and_positions(range_t range, xor_minimiser_and_positions_parameters params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class xor_minimiser_and_positions<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = xor_minimiser_and_positions_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};
    uint64_t seed{0x8F'3F'73'B5'CF'1C'9A'DE};

    uint64_t minimiser_size{};

    size_t range_size{};
    size_t range_position{};

    xor_minimiser_and_positions_parameters params{};
    value_type current{}; // range_position -> position in the window
    value_type cached{};  // range_position -> position in the range

    std::deque<uint64_t> kmer_values_in_window{};

    static inline constexpr uint64_t compute_mask(uint64_t const size) {
        assert(size > 0u);
        assert(size <= 64u);

        if (size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        kmer_value{it.kmer_value},
        range_size{it.range_size},
        range_position{it.range_position},
        params{it.params},
        current{it.current},
        cached{it.cached},
        kmer_values_in_window{it.kmer_values_in_window}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   xor_minimiser_and_positions_parameters params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * params.minimiser_size)},
        minimiser_size{params.minimiser_size},
        range_size{range_size},
        params{std::move(params)}
    {
        if (range_size < params.window_size)
            range_position = range_size + 1u;
        else
            init();
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position > lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        while (!next_minimiser_is_new())
        {}
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        while (!next_minimiser_is_new())
        {}
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return cached;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
    {
        uint64_t const base = seqan3::to_rank(*range_it);
        kmer_value = ((kmer_value << 2) | base) & kmer_mask;
        kmer_value_rev = (kmer_value_rev >> 2) | ((base^0b11) << 2*(minimiser_size-1));
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes)
            kmer_values_in_window.pop_front();

        const uint64_t kmerhash = std::min<uint64_t>(kmer_value, kmer_value_rev) ^ seed;
        
        kmer_values_in_window.push_back(kmerhash);
    }

    void find_minimiser_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values_in_window, std::less_equal<uint64_t>{});
        current.minimiser_value = *minimiser_it;
        current.range_position = std::distance(std::begin(kmer_values_in_window), minimiser_it);
    }

    void init()
    {
        seed = params.seed & kmer_mask;
        minimiser_size = params.minimiser_size;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        kmer_value = ((kmer_value << 2) | new_rank);
        kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        for (size_t i = 1u; i < minimiser_size; ++i) {
            ++range_position;
            ++range_it;

            new_rank = seqan3::to_rank(*range_it);
            kmer_value = ((kmer_value << 2) | new_rank);
            kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        }

        const uint64_t kmerhash = std::min<uint64_t>(kmer_value, kmer_value_rev) ^ seed;
        
        kmer_values_in_window.push_back(kmerhash);

        for (size_t i = minimiser_size; i < params.window_size; ++i)
            next_window<pop_first::no>();

        find_minimiser_in_window();

        while (!next_minimiser_is_new())
        {}
    }

    void update_cache()
    {
        cached.minimiser_value = current.minimiser_value;
        cached.range_position = range_position - params.window_size - current.occurrences;
        cached.occurrences = current.occurrences + 1;
    }

    bool next_minimiser_is_new()
    {
        // If we reached the end of the range, we are done.
        if (range_position + 1 >= range_size)
        {
            ++range_position;
            update_cache();
            return true;
        }

        next_window<pop_first::yes>();

        // The minimiser left the window.
        if (current.range_position == 0)
        {
            update_cache();
            find_minimiser_in_window();
            bool const same_value = current.minimiser_value == cached.minimiser_value;
            current.occurrences *= same_value;
            current.occurrences += same_value;
            return !same_value;
        }

        // Update minimiser if the new kmer value is smaller than the current minimiser.
        if (uint64_t new_kmer_value = kmer_values_in_window.back(); new_kmer_value < current.minimiser_value)
        {
            update_cache();
            current.minimiser_value = new_kmer_value;
            current.range_position = kmer_values_in_window.size() - 1u;
            current.occurrences = 0;
            return true;
        }

        --current.range_position;
        ++current.occurrences;
        return false;
    }
};

template <std::ranges::viewable_range rng_t>
xor_minimiser_and_positions(rng_t &&, xor_minimiser_and_positions_parameters &&)
    -> xor_minimiser_and_positions<std::views::all_t<rng_t>>;

struct xor_minimiser_and_positions_fn
{
    constexpr auto operator()(xor_minimiser_and_positions_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, xor_minimiser_and_positions_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.minimiser_size == 0u)
            throw std::invalid_argument{"minimiser_size must be > 0."};
        if (params.minimiser_size > 32u)
            throw std::invalid_argument{"minimiser_size must be <= 32."};
        if (params.window_size < params.minimiser_size)
            throw std::invalid_argument{"window_size must be >= minimiser_size."};

        return xor_minimiser_and_positions{std::forward<range_t>(range), std::move(params)};
    }
};

}

namespace srindex::views
{

inline constexpr auto xor_minimiser_and_positions = srindex::detail::xor_minimiser_and_positions_fn{};

}





namespace srindex
{

struct two_minimisers_and_occurence_hash_parameters
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed1{};
    uint64_t seed2{};
};

struct two_minimisers_and_occurence_hash_result
{
    uint64_t new_minimiser1_value;
    uint64_t new_minimiser2_value;
    uint64_t minimiser1_value;
    uint64_t minimiser2_value;
    size_t occurrences2;
    size_t new_occurrences2;
};

}

namespace srindex::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class two_minimisers_and_occurence_hash : public std::ranges::view_interface<two_minimisers_and_occurence_hash<range_t>>
{
private:
    range_t range{};
    two_minimisers_and_occurence_hash_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    two_minimisers_and_occurence_hash()
        requires std::default_initializable<range_t>
    = default;
    two_minimisers_and_occurence_hash(two_minimisers_and_occurence_hash const & rhs) = default;
    two_minimisers_and_occurence_hash(two_minimisers_and_occurence_hash && rhs) = default;
    two_minimisers_and_occurence_hash & operator=(two_minimisers_and_occurence_hash const & rhs) = default;
    two_minimisers_and_occurence_hash & operator=(two_minimisers_and_occurence_hash && rhs) = default;
    ~two_minimisers_and_occurence_hash() = default;

    explicit two_minimisers_and_occurence_hash(range_t range, two_minimisers_and_occurence_hash_parameters params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class two_minimisers_and_occurence_hash<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = two_minimisers_and_occurence_hash_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t window_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t seed1{};
    uint64_t seed2{};
    uint8_t minimisers_in_window{};
    uint64_t minimiser_size{};
    uint64_t window_size{};

    size_t minimiser1_position{};
    size_t minimiser2_position{};

    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};

    std::deque<uint64_t> kmer_values1_in_window{};
    std::deque<uint64_t> kmer_values2_in_window{};

    static inline constexpr uint64_t compute_mask(uint64_t const size)
    {
        assert(size > 0u);
        assert(size <= 64u);

        if(size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        window_mask{it.window_mask},
        kmer_value{it.kmer_value},
        kmer_value_rev{it.kmer_value_rev},
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current},
        kmer_values1_in_window{it.kmer_values1_in_window},
        kmer_values2_in_window{it.kmer_values2_in_window},
        minimiser1_position{it.minimiser1_position},
        minimiser2_position{it.minimiser2_position}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   two_minimisers_and_occurence_hash_parameters const & params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * params.minimiser_size)},
        window_mask{compute_mask(2u * params.window_size)},
        range_size{range_size}
    {
        if (range_size < params.window_size)
            range_position = range_size;
        else
            init(params);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position == lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        while (!next_minimiser_is_new())
        {}
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        while (!next_minimiser_is_new())
        {}
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return current;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
    {
        uint64_t const base = seqan3::to_rank(*range_it);
        kmer_value = ((kmer_value << 2) | base) & kmer_mask;
        kmer_value_rev = (kmer_value_rev >> 2) | ((base^0b11) << 2*(minimiser_size-1));
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes) {
            kmer_values1_in_window.pop_front();
            kmer_values2_in_window.pop_front();
        }

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash2 = murmurhash2_64::hash(canonical_kmer, seed2) & kmer_mask;
        
        kmer_values1_in_window.push_back(kmerhash1);
        kmer_values2_in_window.push_back(kmerhash2);
    }

    void find_minimiser1_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values1_in_window, std::less_equal<uint64_t>{});
        current.new_minimiser1_value = *minimiser_it;
        minimiser1_position = std::distance(std::begin(kmer_values1_in_window), minimiser_it);
    }

    void find_minimiser2_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values2_in_window, std::less_equal<uint64_t>{});
        current.new_minimiser2_value = *minimiser_it;
        minimiser2_position = std::distance(std::begin(kmer_values2_in_window), minimiser_it);
    }

    void init(two_minimisers_and_occurence_hash_parameters const & params)
    {
        minimiser_size = params.minimiser_size;
        window_size = params.window_size;
        minimisers_in_window = window_size - minimiser_size;
        seed1 = params.seed1 & kmer_mask;
        seed2 = params.seed2;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        kmer_value = ((kmer_value << 2) | new_rank);
        kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        for (size_t i = 1u; i < minimiser_size; ++i) {
            ++range_position;
            ++range_it;

            new_rank = seqan3::to_rank(*range_it);
            kmer_value = ((kmer_value << 2) | new_rank);
            kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        }

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash2 = murmurhash2_64::hash(canonical_kmer, seed2) & kmer_mask;

        kmer_values1_in_window.push_back(kmerhash1);
        kmer_values2_in_window.push_back(kmerhash2);

        for (size_t i = params.minimiser_size; i < params.window_size; ++i)
            next_window<pop_first::no>();

        find_minimiser1_in_window();
        find_minimiser2_in_window();

        while (!next_minimiser_is_new()) {}
    }
    

    bool next_minimiser_is_new()
    {
        if(range_position + 1 == range_size) {
            current.minimiser1_value = current.new_minimiser1_value;
            current.minimiser2_value = current.new_minimiser2_value;
            current.occurrences2 = current.new_occurrences2 + 1;
            return ++range_position;
        }

        next_window<pop_first::yes>();

        bool new_min1 = false;
        bool new_min2 = false;

        if(minimiser1_position == 0) {
            current.minimiser1_value = current.new_minimiser1_value;
            find_minimiser1_in_window();
            new_min1 = current.minimiser1_value != current.new_minimiser1_value;
        }

        if(minimiser2_position == 0) {
            current.minimiser2_value = current.new_minimiser2_value;
            current.occurrences2 = current.new_occurrences2 + 1;
            find_minimiser2_in_window();
            bool const same_value = current.minimiser2_value == current.new_minimiser2_value;
            current.new_occurrences2 *= same_value;
            current.new_occurrences2 += same_value;
            new_min2 = !same_value;
        }

        if(uint64_t new_kmer_value = kmer_values1_in_window.back(); new_kmer_value < current.minimiser1_value) {
            current.minimiser1_value = current.new_minimiser1_value;
            current.new_minimiser1_value = new_kmer_value;
            minimiser1_position = minimisers_in_window;
            new_min1 = true;
        }

        if(uint64_t new_kmer_value = kmer_values2_in_window.back(); new_kmer_value < current.minimiser2_value) {
            current.minimiser2_value = current.new_minimiser2_value;
            current.new_minimiser2_value = new_kmer_value;
            current.occurrences2 = current.new_occurrences2 + 1;
            current.new_occurrences2 = 0;
            minimiser2_position = minimisers_in_window;
            new_min2 = true;
        }

        if(!new_min1)
            --minimiser1_position;
        if(!new_min2) {
            ++current.new_occurrences2;
            --minimiser2_position;
        }

        return new_min1 || new_min2;
    }

};


template <std::ranges::viewable_range rng_t>
two_minimisers_and_occurence_hash(rng_t &&, two_minimisers_and_occurence_hash_parameters &&)
    -> two_minimisers_and_occurence_hash<std::views::all_t<rng_t>>;

struct two_minimisers_and_occurence_hash_fn
{
    constexpr auto operator()(two_minimisers_and_occurence_hash_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, two_minimisers_and_occurence_hash_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.minimiser_size == 0u)
            throw std::invalid_argument{"minimiser_size must be > 0."};
        if (params.minimiser_size > 32u)
            throw std::invalid_argument{"minimiser_size must be <= 32."};
        if (params.window_size == 0u)
            throw std::invalid_argument{"window_size must be > 0."};
        if (params.window_size > 32u)
            throw std::invalid_argument{"window_size must be <= 32."};
        if (params.window_size < params.minimiser_size)
            throw std::invalid_argument{"window_size must be >= minimiser_size."};

        return two_minimisers_and_occurence_hash{range, std::move(params)};
    }
};

}

namespace srindex::views
{

inline constexpr auto two_minimisers_and_occurence_hash = srindex::detail::two_minimisers_and_occurence_hash_fn{};

}




namespace srindex
{

struct two_minimisers_hash_parameters
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed1{};
    uint64_t seed2{};
};

struct two_minimisers_hash_result
{
    uint64_t minimiser1_value;
    uint64_t minimiser2_value;
};

}

namespace srindex::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class two_minimisers_hash : public std::ranges::view_interface<two_minimisers_hash<range_t>>
{
private:
    range_t range{};
    two_minimisers_hash_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    two_minimisers_hash()
        requires std::default_initializable<range_t>
    = default;
    two_minimisers_hash(two_minimisers_hash const & rhs) = default;
    two_minimisers_hash(two_minimisers_hash && rhs) = default;
    two_minimisers_hash & operator=(two_minimisers_hash const & rhs) = default;
    two_minimisers_hash & operator=(two_minimisers_hash && rhs) = default;
    ~two_minimisers_hash() = default;

    explicit two_minimisers_hash(range_t range, two_minimisers_hash_parameters params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class two_minimisers_hash<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = two_minimisers_hash_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t window_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t seed1{};
    uint64_t seed2{};
    uint8_t minimisers_in_window{};
    uint64_t minimiser_size{};
    uint64_t window_size{};

    size_t minimiser1_position{};
    size_t minimiser2_position{};

    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};

    std::deque<uint64_t> kmer_values1_in_window{};
    std::deque<uint64_t> kmer_values2_in_window{};

    static inline constexpr uint64_t compute_mask(uint64_t const size)
    {
        assert(size > 0u);
        assert(size <= 64u);

        if(size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        window_mask{it.window_mask},
        kmer_value{it.kmer_value},
        kmer_value_rev{it.kmer_value_rev},
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current},
        kmer_values1_in_window{it.kmer_values1_in_window},
        kmer_values2_in_window{it.kmer_values2_in_window},
        minimiser1_position{it.minimiser1_position},
        minimiser2_position{it.minimiser2_position}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   two_minimisers_hash_parameters const & params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * params.minimiser_size)},
        window_mask{compute_mask(2u * params.window_size)},
        range_size{range_size}
    {
        if (range_size < params.window_size)
            range_position = range_size;
        else
            init(params);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position == lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        while (!next_minimiser_is_new())
        {}
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        while (!next_minimiser_is_new())
        {}
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return current;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
    {
        uint64_t const base = seqan3::to_rank(*range_it);
        kmer_value = ((kmer_value << 2) | base) & kmer_mask;
        kmer_value_rev = (kmer_value_rev >> 2) | ((base^0b11) << 2*(minimiser_size-1));
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes) {
            kmer_values1_in_window.pop_front();
            kmer_values2_in_window.pop_front();
        }

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash2 = murmurhash2_64::hash(canonical_kmer, seed2) & kmer_mask;
        
        kmer_values1_in_window.push_back(kmerhash1);
        kmer_values2_in_window.push_back(kmerhash2);
    }

    void find_minimiser1_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values1_in_window, std::less_equal<uint64_t>{});
        current.minimiser1_value = *minimiser_it;
        minimiser1_position = std::distance(std::begin(kmer_values1_in_window), minimiser_it);
    }

    void find_minimiser2_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values2_in_window, std::less_equal<uint64_t>{});
        current.minimiser2_value = *minimiser_it;
        minimiser2_position = std::distance(std::begin(kmer_values2_in_window), minimiser_it);
    }

    void init(two_minimisers_hash_parameters const & params)
    {
        minimiser_size = params.minimiser_size;
        window_size = params.window_size;
        minimisers_in_window = window_size - minimiser_size;
        seed1 = params.seed1 & kmer_mask;
        seed2 = params.seed2;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        kmer_value = ((kmer_value << 2) | new_rank);
        kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        for (size_t i = 1u; i < minimiser_size; ++i) {
            ++range_position;
            ++range_it;

            new_rank = seqan3::to_rank(*range_it);
            kmer_value = ((kmer_value << 2) | new_rank);
            kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        }

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash2 = murmurhash2_64::hash(canonical_kmer, seed2) & kmer_mask;

        kmer_values1_in_window.push_back(kmerhash1);
        kmer_values2_in_window.push_back(kmerhash2);

        for (size_t i = params.minimiser_size; i < params.window_size; ++i)
            next_window<pop_first::no>();

        find_minimiser1_in_window();
        find_minimiser2_in_window();
    }
    

    bool next_minimiser_is_new()
    {
        if (range_position + 1 == range_size)
            return ++range_position;

        next_window<pop_first::yes>();

        bool new_min1 = false;
        bool new_min2 = false;

        if(minimiser1_position == 0) {
            find_minimiser1_in_window();
            new_min1 = true;
        }

        if(minimiser2_position == 0) {
            find_minimiser2_in_window();
            new_min2 = true;
        }

        if(uint64_t new_kmer_value = kmer_values1_in_window.back(); new_kmer_value < current.minimiser1_value) {
            current.minimiser1_value = new_kmer_value;
            minimiser1_position = minimisers_in_window;
            new_min1 = true;
        }

        if(uint64_t new_kmer_value = kmer_values2_in_window.back(); new_kmer_value < current.minimiser2_value) {
            current.minimiser2_value = new_kmer_value;
            minimiser2_position = minimisers_in_window;
            new_min2 = true;
        }

        if(!new_min1)
            --minimiser1_position;
        if(!new_min2) {
            --minimiser2_position;
        }

        return new_min1 || new_min2;
    }

};


template <std::ranges::viewable_range rng_t>
two_minimisers_hash(rng_t &&, two_minimisers_hash_parameters &&)
    -> two_minimisers_hash<std::views::all_t<rng_t>>;

struct two_minimisers_hash_fn
{
    constexpr auto operator()(two_minimisers_hash_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, two_minimisers_hash_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.minimiser_size == 0u)
            throw std::invalid_argument{"minimiser_size must be > 0."};
        if (params.minimiser_size > 32u)
            throw std::invalid_argument{"minimiser_size must be <= 32."};
        if (params.window_size == 0u)
            throw std::invalid_argument{"window_size must be > 0."};
        if (params.window_size > 32u)
            throw std::invalid_argument{"window_size must be <= 32."};
        if (params.window_size < params.minimiser_size)
            throw std::invalid_argument{"window_size must be >= minimiser_size."};

        return two_minimisers_hash{range, std::move(params)};
    }
};

}

namespace srindex::views
{

inline constexpr auto two_minimisers_hash = srindex::detail::two_minimisers_hash_fn{};

}
