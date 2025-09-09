

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


std::vector<uint64_t> text_kmers(RSIndexComp &index) {
    constexpr uint64_t n = 1000000;

    std::uniform_int_distribution<uint32_t> distr;
    std::mt19937 m_rand(1);
    std::vector<std::uint64_t> kmers;
    kmers.reserve(n);

    const uint64_t no_unitigs = index.number_unitigs();
    const uint64_t k = index.getk();
    for (uint64_t i = 0; i < n; ++i) {
        const uint64_t unitig_id = distr(m_rand) % no_unitigs;
        const uint64_t offset = distr(m_rand) % index.unitig_size(unitig_id);
        // const uint64_t offset = 0;
        // std::cout << "endpoint: " << endpoint << " next endpoint: " << next_endpoint << " unitig_size: " << next_endpoint-endpoint-index.getk()+1 << " offset: " << offset << '\n';
        // std::cout << "unitig_id: " << unitig_id << " unitig_size: " << index.unitig_size(unitig_id) << " offset: " << offset << '\n';
        const uint64_t kmer = index.access(unitig_id, offset);

        /* transform 50% of the kmers into their reverse complements */
        if ((i & 1) == 0)
            kmers.push_back(crc(kmer, k));
        else
            kmers.push_back(kmer);
    }

    return kmers;
}


static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if(size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << (size)) - 1u;
}


std::vector<uint64_t> rand_kmers(const uint64_t k) {
    constexpr uint64_t n = 1000000;
    const uint64_t mask = compute_mask(2u * k);

    std::uniform_int_distribution<uint64_t> distr;
    std::mt19937_64 m_rand(1);
    std::vector<uint64_t> kmers;
    kmers.reserve(n);

    for (uint64_t i = 0; i < n; ++i) {
        const uint64_t kmer = distr(m_rand) & mask;
        kmers.push_back(kmer);
    }

    return kmers;
}
