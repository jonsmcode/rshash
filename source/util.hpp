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

static inline constexpr uint64_t crc(uint64_t x, uint64_t k) {
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
