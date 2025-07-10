// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/debug_stream.hpp>


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

namespace srindex::util {

static inline uint64_t crc(uint64_t x, uint64_t k) {
    // assert(k <= 32);
    uint64_t c = ~x;

    /* swap byte order */
    uint64_t res = __builtin_bswap64(c);

    /* Swap nuc order in bytes */
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;              // ...0000.1111.0000.1111
    const uint64_t c2 = 0x3333333333333333;              // ...0011.0011.0011.0011
    res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4);  // swap 2-nuc order in bytes
    res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2);  // swap nuc order in 2-nuc

    /* Realign to the right */
    res >>= 64 - 2 * k;

    return res;
}

}


namespace srindex
{

struct minimiser_and_window_hash_parameters
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed{};
};

struct minimiser_and_window_hash_result
{
    uint64_t minimiser_value;
    uint64_t window_value;
    uint64_t window_value_rev;
};

}

namespace srindex::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class minimiser_and_window_hash : public std::ranges::view_interface<minimiser_and_window_hash<range_t>>
{
private:
    range_t range{};
    minimiser_and_window_hash_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    minimiser_and_window_hash()
        requires std::default_initializable<range_t>
    = default;
    minimiser_and_window_hash(minimiser_and_window_hash const & rhs) = default;
    minimiser_and_window_hash(minimiser_and_window_hash && rhs) = default;
    minimiser_and_window_hash & operator=(minimiser_and_window_hash const & rhs) = default;
    minimiser_and_window_hash & operator=(minimiser_and_window_hash && rhs) = default;
    ~minimiser_and_window_hash() = default;

    explicit minimiser_and_window_hash(range_t range, minimiser_and_window_hash_parameters params) :
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
class minimiser_and_window_hash<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = minimiser_and_window_hash_result;
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
    // uint64_t seed{0x8F'3F'73'B5'CF'1C'9A'DE};
    uint64_t seed{1};
    uint8_t minimisers_in_window{};
    uint64_t minimiser_size{};
    uint64_t window_size{};

    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};
    size_t minimiser_position{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};

    std::deque<uint64_t> kmer_values_in_window{};

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
        minimiser_position{it.minimiser_position},
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current},
        kmer_values_in_window{it.kmer_values_in_window}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   minimiser_and_window_hash_parameters const & params) :
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
        kmer_value = current.window_value & kmer_mask;
        // kmer_value = ((kmer_value << 2) | new_rank) & kmer_mask;
        current.window_value_rev = (current.window_value_rev >> 2) | ((new_rank^0b11) << 2*(window_size-1));
        kmer_value_rev = current.window_value_rev >> 2*(window_size - minimiser_size);
        // kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes)
            kmer_values_in_window.pop_front();

        const uint64_t kmerhash = murmurhash2_64::hash(std::min<uint64_t>(kmer_value, kmer_value_rev), seed) & kmer_mask;
        kmer_values_in_window.push_back(kmerhash);
    }

    void find_minimiser_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values_in_window, std::less_equal<uint64_t>{});
        current.minimiser_value = *minimiser_it;
        minimiser_position = std::distance(std::begin(kmer_values_in_window), minimiser_it);
    }

    void init(minimiser_and_window_hash_parameters const & params)
    {
        // seed = params.seed & kmer_mask;
        seed = params.seed;
        minimiser_size = params.minimiser_size;
        window_size = params.window_size;
        minimisers_in_window = window_size - minimiser_size;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        current.window_value <<= 2;
        current.window_value |= new_rank;
        current.window_value_rev >>= 2;
        current.window_value_rev |= ((new_rank^0b11) << 2*(window_size-1));
        // kmer_value = ((kmer_value << 2) | new_rank) & kmer_mask;
        // kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        for (size_t i = 1u; i < params.minimiser_size; ++i) {
            ++range_position;
            ++range_it;
            new_rank = seqan3::to_rank(*range_it);
            current.window_value <<= 2;
            current.window_value |= new_rank;
            current.window_value_rev >>= 2;
            current.window_value_rev |= ((new_rank^0b11) << 2*(window_size-1));
            // kmer_value = ((kmer_value << 2) | new_rank) & kmer_mask;
            // kmer_value_rev = (kmer_value_rev >> 2) | ((new_rank^0b11) << 2*(minimiser_size-1));
        }
        kmer_value = current.window_value & kmer_mask;
        kmer_value_rev = current.window_value_rev >> 2*(window_size - minimiser_size);

        const uint64_t kmerhash = murmurhash2_64::hash(std::min<uint64_t>(kmer_value, kmer_value_rev), seed) & kmer_mask;
        kmer_values_in_window.push_back(kmerhash);

        for (size_t i = minimiser_size; i < window_size; ++i)
            next_window<pop_first::no>();

        find_minimiser_in_window();
    }

    bool next_minimiser_is_new()
    {
        // If we reached the end of the range, we are done.
        if (range_position + 1 == range_size)
            return ++range_position;

        next_window<pop_first::yes>();

        // The minimiser left the window.
        if (minimiser_position == 0)
        {
            find_minimiser_in_window();
            return true;
        }

        // Update minimiser if the new kmer value is smaller than the current minimiser.
        if (uint64_t new_kmer_value = kmer_values_in_window.back(); new_kmer_value < current.minimiser_value)
        {
            current.minimiser_value = new_kmer_value;
            minimiser_position = minimisers_in_window;
            return true;
        }

        --minimiser_position;
        return true; // workaround to retrieve all kmers always return true
    }
};


template <std::ranges::viewable_range rng_t>
minimiser_and_window_hash(rng_t &&, minimiser_and_window_hash_parameters &&)
    -> minimiser_and_window_hash<std::views::all_t<rng_t>>;

struct minimiser_and_window_hash_fn
{
    constexpr auto operator()(minimiser_and_window_hash_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, minimiser_and_window_hash_parameters params) const
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

        return minimiser_and_window_hash{range, std::move(params)};
    }
};

}

namespace srindex::views
{

inline constexpr auto minimiser_and_window_hash = srindex::detail::minimiser_and_window_hash_fn{};

}



namespace srindex
{

struct minimiser_hash_and_positions_parameters
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed{};
};

struct minimiser_hash_and_positions_result
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
class minimiser_hash_and_positions : public std::ranges::view_interface<minimiser_hash_and_positions<range_t>>
{
private:
    range_t range{};
    minimiser_hash_and_positions_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    minimiser_hash_and_positions()
        requires std::default_initializable<range_t>
    = default;
    minimiser_hash_and_positions(minimiser_hash_and_positions const & rhs) = default;
    minimiser_hash_and_positions(minimiser_hash_and_positions && rhs) = default;
    minimiser_hash_and_positions & operator=(minimiser_hash_and_positions const & rhs) = default;
    minimiser_hash_and_positions & operator=(minimiser_hash_and_positions && rhs) = default;
    ~minimiser_hash_and_positions() = default;

    explicit minimiser_hash_and_positions(range_t range, minimiser_hash_and_positions_parameters params) :
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
class minimiser_hash_and_positions<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = minimiser_hash_and_positions_result;
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
    // uint64_t seed{0x8F'3F'73'B5'CF'1C'9A'DE};
    uint64_t seed{1};

    uint64_t minimiser_size{};

    size_t range_size{};
    size_t range_position{};

    minimiser_hash_and_positions_parameters params{};
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
                   minimiser_hash_and_positions_parameters params) :
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

        const uint64_t kmerhash = murmurhash2_64::hash(std::min<uint64_t>(kmer_value, kmer_value_rev), seed) & kmer_mask;
        
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
        // seed = params.seed & kmer_mask;
        seed = params.seed;
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

        const uint64_t kmerhash = murmurhash2_64::hash(std::min<uint64_t>(kmer_value, kmer_value_rev), seed) & kmer_mask;
        
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
minimiser_hash_and_positions(rng_t &&, minimiser_hash_and_positions_parameters &&)
    -> minimiser_hash_and_positions<std::views::all_t<rng_t>>;

struct minimiser_hash_and_positions_fn
{
    constexpr auto operator()(minimiser_hash_and_positions_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, minimiser_hash_and_positions_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.minimiser_size == 0u)
            throw std::invalid_argument{"minimiser_size must be > 0."};
        if (params.minimiser_size > 32u)
            throw std::invalid_argument{"minimiser_size must be <= 32."};
        if (params.window_size < params.minimiser_size)
            throw std::invalid_argument{"window_size must be >= minimiser_size."};

        return minimiser_hash_and_positions{std::forward<range_t>(range), std::move(params)};
    }
};

}

namespace srindex::views
{

inline constexpr auto minimiser_hash_and_positions = srindex::detail::minimiser_hash_and_positions_fn{};

}


