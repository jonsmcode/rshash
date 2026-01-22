#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <gtl/phmap.hpp>

#include "compact_vector.hpp"
#include "EliasFano.hpp"

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;


// #pragma once


static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if (size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << size) - 1u;
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

inline std::vector<uint64_t> pack_dna4_to_uint64(const std::vector<std::vector<seqan3::dna4>> & input)
{
    auto ranks = input | std::views::join | seqan3::views::to_rank;

    std::vector<uint64_t> packed;
    uint64_t word = 0;
    size_t shift = 0;

    for (uint8_t r : ranks)
    {
        word |= uint64_t(r) << shift; // pack 2 bits per base
        shift += 2;

        if (shift == 64)
        {
            packed.push_back(word);
            word = 0;
            shift = 0;
        }
    }

    if (shift != 0) // last partial word
        packed.push_back(word);

    packed.push_back(0);

    return packed;
}


// const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed1 = 1;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed3 = 0xE5'9A'38'5F'03'76'C9'F6;


struct SkmerInfo {
    uint64_t position;
    uint64_t unitig_begin;
    uint64_t unitig_end;
};


class RSHash1
{
private:
    uint64_t k, m1, m_thres1;
    uint64_t span;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r1;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    pthash::compact_vector offsets1;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    inline bool check(uint64_t*, std::array<uint64_t, 2>*, const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, const uint64_t);
    inline bool check_minimiser_pos(uint64_t *, const SkmerInfo&, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, size_t &, size_t &, size_t &);
    inline void refill_buffer(uint64_t *, uint64_t*, SkmerInfo*, size_t, size_t, const uint64_t, const uint64_t);
    inline bool lookup_buffer(uint64_t *, SkmerInfo *, const size_t, const uint64_t, const uint64_t, size_t &, const size_t, const size_t, bool &, size_t &, size_t &);
    inline bool extend_in_text(size_t&, size_t, size_t, bool, const uint64_t, const uint64_t, const uint64_t);
    const inline uint64_t get_word64(uint64_t pos);
    const inline uint64_t get_base(uint64_t pos);


public:
    RSHash1();
    RSHash1(uint8_t const k, uint8_t const m1, uint8_t const m_thres1);
    uint8_t getk() { return k; }
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    std::vector<uint64_t> rand_text_kmers(const uint64_t);
    uint64_t access(const uint64_t, const size_t);
    uint64_t lookup(const std::vector<uint64_t>&, bool verbose);
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSHash2
{
private:
    uint64_t k, m1, m_thres1, m2, m_thres2;
    uint64_t span1, span2;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r1, r2;
    bit_vector s1, s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select, s2_select;
    pthash::compact_vector offsets1, offsets2;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    template<int level>
    inline bool check(uint64_t*, std::array<uint64_t, 2>*, const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, const uint64_t);
    template<int level>
    inline void refill_buffer(uint64_t *, SkmerInfo *, size_t, size_t, const uint64_t, const uint64_t);
    inline bool check_minimiser_pos(uint64_t *, const SkmerInfo&, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, size_t &, size_t &, size_t &);
    template<int level>
    inline bool lookup_buffer(uint64_t *, SkmerInfo *, const size_t, const uint64_t,  const uint64_t, size_t &, const size_t, const size_t, bool &, size_t &, size_t &);
    inline bool extend_in_text(size_t&, size_t, size_t, bool, const uint64_t, const uint64_t);
    const inline uint64_t get_word64(uint64_t pos);
    const inline uint64_t get_base(uint64_t pos);


public:
    RSHash2();
    RSHash2(uint8_t const k, uint8_t const m1, uint8_t const m_thres1,
                 uint8_t const m2, uint8_t const m_thres2);
    uint8_t getk() { return k; }
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    std::vector<uint64_t> rand_text_kmers(const uint64_t);
    uint64_t access(const uint64_t, const size_t);
    uint64_t lookup(const std::vector<uint64_t>&, bool verbose);
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSHash3
{
private:
    uint64_t k, m1, m_thres1, m2, m_thres2, m3, m_thres3;
    uint64_t span1, span2, span3;
    bit_vector r1;
    rank_support_v<1> r1_rank;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r2, r3;
    bit_vector s1, s2, s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select, s2_select, s3_select;
    pthash::compact_vector offsets1, offsets2, offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    template<int level>
    inline bool check(uint64_t*, std::array<uint64_t, 2>*, const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, const uint64_t);
    template<int level>
    inline void refill_buffer(uint64_t *, SkmerInfo *, size_t, size_t, const uint64_t, const uint64_t);
    inline bool check_minimiser_pos(uint64_t *, const SkmerInfo&, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, size_t &, size_t &, size_t &);
    template<int level>
    inline bool lookup_buffer(uint64_t *, SkmerInfo *, const size_t, const uint64_t,  const uint64_t, size_t &, const size_t, const size_t, bool &, size_t &, size_t &);
    inline bool extend_in_text(size_t&, size_t, size_t, bool, const uint64_t, const uint64_t);
    const inline uint64_t get_word64(uint64_t pos);
    const inline uint64_t get_base(uint64_t pos);


public:
    RSHash3();
    RSHash3(uint8_t const k, uint8_t const m1, uint8_t const m_thres1,
                 uint8_t const m2, uint8_t const m_thres2, uint8_t const m3, uint8_t const m_thres3);
    uint8_t getk() { return k; }
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    std::vector<uint64_t> rand_text_kmers(const uint64_t);
    uint64_t access(const uint64_t, const size_t);
    uint64_t lookup(const std::vector<uint64_t>&, bool verbose);
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSHash3C
{
private:
    uint64_t k, m1, m_thres1, m2, m_thres2, m3, m_thres3;
    uint64_t span1, span2, span3;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r1, r2, r3;
    bit_vector s1, s2, s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select, s2_select, s3_select;
    pthash::compact_vector offsets1, offsets2, offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    template<int level>
    inline bool check(uint64_t*, std::array<uint64_t, 2>*, const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, const uint64_t);
    template<int level>
    inline void refill_buffer(uint64_t *, SkmerInfo *, size_t, size_t, const uint64_t, const uint64_t);
    inline bool check_minimiser_pos(uint64_t *, const SkmerInfo&, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, size_t &, size_t &, size_t &);
    template<int level>
    inline bool lookup_buffer(uint64_t *, SkmerInfo *, const size_t, const uint64_t,  const uint64_t, size_t &, const size_t, const size_t, bool &, size_t &, size_t &);
    inline bool extend_in_text(size_t&, size_t, size_t, bool, const uint64_t, const uint64_t);
    const inline uint64_t get_word64(uint64_t pos);
    const inline uint64_t get_base(uint64_t pos);


public:
    RSHash3C();
    RSHash3C(uint8_t const k, uint8_t const m1, uint8_t const m_thres1,
                 uint8_t const m2, uint8_t const m_thres2, uint8_t const m3, uint8_t const m_thres3);
    uint8_t getk() { return k; }
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    std::vector<uint64_t> rand_text_kmers(const uint64_t);
    uint64_t access(const uint64_t, const size_t);
    uint64_t lookup(const std::vector<uint64_t>&, bool verbose);
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};