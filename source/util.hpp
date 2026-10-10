struct mask128_t {
    uint64_t lo;
    uint64_t hi;
};


static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if (size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << size) - 1u;
}

static inline constexpr uint64_t compute_mask(unsigned x, unsigned z)
{
    assert(x <= z);
    assert(z < 64);

    return compute_mask(z - x + 1) << x;
}

static inline constexpr mask128_t compute_mask128(unsigned n)
{
    if (n == 0)
        return {0, 0};
    if (n < 64)
        return {(1ULL << n) - 1, 0};
    if (n == 64)
        return {UINT64_MAX, 0};
    if (n < 128)
        return {UINT64_MAX, (1ULL << (n - 64)) - 1};
    return {UINT64_MAX, UINT64_MAX};
}

static inline constexpr uint64_t rc64(uint64_t x) {
    x = __builtin_bswap64(~x);

    constexpr uint64_t m4 = 0x0f0f0f0f0f0f0f0f;
    constexpr uint64_t m2 = 0x3333333333333333;

    x = ((x & m4) << 4) | ((x & ~m4) >> 4);
    x = ((x & m2) << 2) | ((x & ~m2) >> 2);

    return x;
}

static inline constexpr uint64_t crc(uint64_t x, uint64_t k) {
    return rc64(x) >> (64 - 2 * k);
}

static inline constexpr mask128_t crc128(mask128_t x, uint64_t k) {
    const uint64_t n = k - 32;

    const uint64_t rhi = rc64(x.lo);
    const uint64_t rlo = crc(x.hi, n);

    return {
        .lo = (rhi << (2 * n)) | rlo,
        .hi = rhi >> (64 - 2 * n)
    };
}
