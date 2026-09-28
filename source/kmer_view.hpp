#pragma once
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>


namespace rshash
{

struct kmerview_parameters
{
    size_t window_size{};
};

struct kmerview_result
{
    uint64_t value{};
    uint64_t value_hi{};
    uint64_t value_rev{};
    uint64_t value_rev_hi{};
};

}

namespace rshash::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class kmerview : public std::ranges::view_interface<kmerview<range_t>>
{
private:
    range_t range{};
    kmerview_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    kmerview()
        requires std::default_initializable<range_t>
    = default;
    kmerview(kmerview const & rhs) = default;
    kmerview(kmerview && rhs) = default;
    kmerview & operator=(kmerview const & rhs) = default;
    kmerview & operator=(kmerview && rhs) = default;
    ~kmerview() = default;

    explicit kmerview(range_t range, kmerview_parameters params) :
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
class kmerview<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = kmerview_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t window_size{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};

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
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   kmerview_parameters const & params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * params.window_size)},
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
        ++range_position;
        ++range_it;
        rolling_hash();
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        ++range_position;
        ++range_it;
        rolling_hash();
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
        current.value = (current.value >> 2) | (new_rank << 2*(window_size-1));
        current.value_rev = ((current.value_rev << 2) | (new_rank^0b11)) & kmer_mask;
    }

    void init(kmerview_parameters const & params)
    {
        window_size = params.window_size;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        current.value = new_rank << 2*(window_size-1);
        current.value_rev = new_rank^0b11;
        for (size_t i = 1u; i < params.window_size; ++i) {
            ++range_position;
            ++range_it;
            rolling_hash();
        }
        
    }

};


template <std::ranges::viewable_range rng_t>
kmerview(rng_t &&, kmerview_parameters &&)
    -> kmerview<std::views::all_t<rng_t>>;

struct kmerview_fn
{
    constexpr auto operator()(kmerview_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, kmerview_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.window_size == 0u)
            throw std::invalid_argument{"window_size must be > 0."};
        if (params.window_size > 32u)
            throw std::invalid_argument{"window_size must be <= 32."};

        return kmerview{range, std::move(params)};
    }
};

}

namespace rshash::views
{

inline constexpr auto kmerview = rshash::detail::kmerview_fn{};

}



namespace rshash::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class longkmerview : public std::ranges::view_interface<longkmerview<range_t>>
{
private:
    range_t range{};
    kmerview_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    longkmerview()
        requires std::default_initializable<range_t>
    = default;
    longkmerview(longkmerview const & rhs) = default;
    longkmerview(longkmerview && rhs) = default;
    longkmerview & operator=(longkmerview const & rhs) = default;
    longkmerview & operator=(longkmerview && rhs) = default;
    ~longkmerview() = default;

    explicit longkmerview(range_t range, kmerview_parameters params) :
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
class longkmerview<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = kmerview_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{};
    uint64_t high_mask{};
    uint64_t high_shift{};
    uint64_t window_size{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};

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
        window_size{it.window_size},
        high_mask{it.high_mask},
        high_shift{it.high_shift},
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   kmerview_parameters const & params) :
        range_it{std::move(range_iterator)},
        high_mask{compute_mask(2 * params.window_size - 64)},
        high_shift{2 * (params.window_size - 33)},
        range_size{range_size},
        window_size{params.window_size}
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
        ++range_position;
        ++range_it;
        rolling_hash();
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        ++range_position;
        ++range_it;
        rolling_hash();
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

    void rolling_hash() noexcept {
        uint64_t const new_rank = seqan3::to_rank(*range_it);

        current.value = (current.value >> 2) | (current.value_hi << 62);
        current.value_hi = (current.value_hi >> 2) | (new_rank << high_shift);

        current.value_rev_hi = ((current.value_rev_hi << 2) | (current.value_rev >> 62)) & high_mask;
        current.value_rev =  (current.value_rev << 2) | (new_rank ^ 0b11);
    }

    void init(kmerview_parameters const & params)
    {
        window_size = params.window_size;

        uint64_t new_rank = seqan3::to_rank(*range_it);

        current.value = 0;
        current.value_hi = new_rank << high_shift;
        current.value_rev = new_rank ^ 0b11;
        current.value_rev_hi = 0;

        for (size_t i = 1u; i < params.window_size; ++i) {
            ++range_position;
            ++range_it;
            rolling_hash();
        }
    }

};


template <std::ranges::viewable_range rng_t>
longkmerview(rng_t &&, kmerview_parameters &&) -> longkmerview<std::views::all_t<rng_t>>;

struct longkmerview_fn
{
    constexpr auto operator()(kmerview_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, kmerview_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.window_size <= 32u)
            throw std::invalid_argument{"window_size must be > 32."};
        if (params.window_size >= 64u)
            throw std::invalid_argument{"window_size must be <= 64."};

        return longkmerview{range, std::move(params)};
    }
};

}

namespace rshash::views
{

inline constexpr auto longkmerview = rshash::detail::longkmerview_fn{};

}


namespace rshash::detail
{

// ============================================================================
// Combined view
// ============================================================================

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t>
          && std::ranges::sized_range<range_t>
class kmer_view :
    public std::ranges::view_interface<kmer_view<range_t>>
{
private:

    using short_view_t = kmerview<range_t>;
    using long_view_t  = longkmerview<range_t>;

    using short_iterator_t =
        std::ranges::iterator_t<short_view_t>;

    using long_iterator_t =
        std::ranges::iterator_t<long_view_t>;

    using iterator_t =
        std::variant<
            short_iterator_t,
            long_iterator_t
        >;

    std::variant<
        short_view_t,
        long_view_t
    > view;


    class basic_iterator
    {
    public:

        iterator_t iterator;


        // IMPORTANT:
        // The underlying kmerview/longkmerview iterators return
        // kmerview_result by value, so we must also return by value.
        kmerview_result operator*() const
        {
            return std::visit(
                [](auto const & it) -> kmerview_result
                {
                    return *it;
                },
                iterator
            );
        }


        basic_iterator & operator++()
        {
            std::visit(
                [](auto & it)
                {
                    ++it;
                },
                iterator
            );

            return *this;
        }


        void operator++(int)
        {
            ++(*this);
        }


        friend bool operator==(
            basic_iterator const & it,
            std::default_sentinel_t sentinel)
        {
            return std::visit(
                [sentinel](auto const & i)
                {
                    return i == sentinel;
                },
                it.iterator
            );
        }


        friend bool operator!=(
            basic_iterator const & it,
            std::default_sentinel_t sentinel)
        {
            return !(it == sentinel);
        }


        friend bool operator==(
            std::default_sentinel_t sentinel,
            basic_iterator const & it)
        {
            return it == sentinel;
        }


        friend bool operator!=(
            std::default_sentinel_t sentinel,
            basic_iterator const & it)
        {
            return !(it == sentinel);
        }
    };


    static auto make_view(
        range_t range,
        kmerview_parameters params)
        -> std::variant<
            short_view_t,
            long_view_t>
    {
        if (params.window_size <= 32u)
        {
            return std::variant<
                short_view_t,
                long_view_t>{
                    std::in_place_type<short_view_t>,
                    std::move(range),
                    params
                };
        }

        return std::variant<
            short_view_t,
            long_view_t>{
                std::in_place_type<long_view_t>,
                std::move(range),
                params
            };
    }


public:

    kmer_view(
        range_t range,
        kmerview_parameters params) :
        view{make_view(std::move(range), params)}
    {}


    basic_iterator begin()
    {
        return std::visit(
            [](auto & v) -> basic_iterator
            {
                return basic_iterator{v.begin()};
            },
            view
        );
    }


    std::default_sentinel_t end() const noexcept
    {
        return {};
    }
};


// ============================================================================
// Combined view deduction guide
// ============================================================================

template <std::ranges::viewable_range rng_t>
kmer_view(rng_t &&, kmerview_parameters &&)
    -> kmer_view<std::views::all_t<rng_t>>;


// ============================================================================
// Combined view adaptor
// ============================================================================

struct kmer_view_fn
{
    constexpr auto operator()(kmerview_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{
            *this,
            std::move(params)
        };
    }


    template <std::ranges::range range_t>
    constexpr auto operator()(
        range_t && range,
        kmerview_parameters params
    ) const
    {
        static_assert(
            std::same_as<
                std::ranges::range_value_t<range_t>,
                seqan3::dna4
            >,
            "Only dna4 supported."
        );

        static_assert(
            std::ranges::sized_range<range_t>,
            "Input range must be a std::ranges::sized_range."
        );


        if (params.window_size == 0u)
        {
            throw std::invalid_argument{
                "window_size must be > 0."
            };
        }


        if (params.window_size >= 64u)
        {
            throw std::invalid_argument{
                "window_size must be < 64."
            };
        }


        return kmer_view{
            std::views::all(std::forward<range_t>(range)),
            std::move(params)
        };
    }
};

} // namespace rshash::detail


// ============================================================================
// Public adaptor
// ============================================================================

namespace rshash::views
{

inline constexpr auto kmer_view =
    rshash::detail::kmer_view_fn{};

} // namespace rshash::views