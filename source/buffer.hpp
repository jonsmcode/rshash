

inline void UnitigsDictionaryHash2::fill_buffer(std::vector<uint64_t> &buffer, const uint64_t mask, const size_t p, const size_t q)
{
    // Precompute min(next_endpoint, o+k+k) for the range [p, q)
    std::vector<size_t> endpoints(q - p);
    size_t total_elements = 0;
    for (size_t i = 0; i < q - p; ++i) {
        size_t o = offsets1[p + i];
        size_t next_endpoint = endpoints_select(endpoints_rank(o + 1) + 1);
        size_t e = o + span + k;
        if (e > next_endpoint)
            e = next_endpoint;
        endpoints[i] = e;
        total_elements += e - o;
    }

    buffer.resize(total_elements);
    size_t buf_idx = 0;

    for (size_t i = 0; i < q - p; ++i) {
        uint64_t hash = 0;
        size_t o = offsets1[p + i];
        size_t e = endpoints[i];
        for (uint64_t j = o; j < o + k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
        }
        buffer[buf_idx++] = hash;
        for (size_t j = o + k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
            buffer[buf_idx++] = hash;
        }
    }
}
