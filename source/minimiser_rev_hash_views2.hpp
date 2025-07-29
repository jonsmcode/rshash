// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>


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
    // uint64_t minimiser_mask{std::numeric_limits<uint64_t>::max()};
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
        // seed1{params.seed1},
        // seed2{params.seed2},
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

        if constexpr (pop == pop_first::yes) {
            kmer_values1_in_window.pop_front();
            kmer_values2_in_window.pop_front();
        }

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        // const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash1 = murmurhash2_64::hash(canonical_kmer, seed1) & kmer_mask;
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
        // seed1 = params.seed1 & kmer_mask;
        seed1 = params.seed1;
        // seed2 = params.seed2 & kmer_mask;
        seed2 = params.seed2;

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

        const uint64_t canonical_kmer = std::min<uint64_t>(kmer_value, kmer_value_rev);
        // const uint64_t kmerhash1 = canonical_kmer ^ seed1;
        const uint64_t kmerhash1 = murmurhash2_64::hash(canonical_kmer, seed1) & kmer_mask;
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




