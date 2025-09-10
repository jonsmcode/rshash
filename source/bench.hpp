

static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if(size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << (size)) - 1u;
}


std::vector<uint64_t> rand_kmers(const uint64_t n, const uint64_t k) {
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
