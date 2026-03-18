#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

#include <sux/bits/SimpleSelect.hpp>
#include <gtl/phmap.hpp>

#include "compact_vector.hpp"
#include "EliasFano.hpp"
#include "minimiser_views.hpp"

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;


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

inline std::vector<uint64_t> pack_dna4_to_uint64(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    auto ranks = input | std::views::join | seqan3::views::to_rank;

    std::vector<uint64_t> packed;
    packed.push_back(UINT64_MAX); // padding for left extension

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

    packed.push_back(UINT64_MAX); // padding for right extension

    return packed;
}


const uint64_t seed1 = 1;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed3 = 0xE5'9A'38'5F'03'76'C9'F6;


class RSHash1
{
private:
    uint64_t k, m1, m_thres1;
    uint64_t span;
    uint64_t kmermask, mmermask;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r1;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bits::compact_vector offsets1;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    uint64_t no_text_kmers;
    mixer_64 m_hasher;
    size_t mark_sequences(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &);
    uint64_t get_unfrequent_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &, std::vector<uint64_t> &, std::vector<uint8_t> &);
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> get_frequent_skmers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &);
    void mark_minimizer_occurences(const size_t, const std::vector<uint8_t> &);
    void fill_minimizer_offsets(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &, std::vector<uint8_t> &, const size_t, const size_t, const size_t);
    inline uint64_t find_minimiser(const uint64_t, const uint64_t, size_t &, size_t &);
    inline void update_minimiser(const uint64_t, const uint64_t, uint64_t&, size_t &, size_t &);
    inline bool check(const uint64_t, const uint64_t, uint64_t*, const size_t, const size_t, const size_t, const size_t);
    inline bool check_minimiser_pos(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    inline bool check_minimiser_pos2(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    inline bool check_overlap(uint64_t, uint64_t, uint64_t &, uint64_t &);
    inline void fill_buffer(uint64_t*, uint64_t*, size_t, size_t, const uint64_t);
    inline bool lookup_buffer(uint64_t*, uint64_t*, const size_t, const uint64_t, const uint64_t, uint64_t &, const size_t, const size_t, bool &, uint64_t &, uint64_t &);
    inline bool extend_in_text(uint64_t&, uint64_t, uint64_t, bool, const uint64_t, const uint64_t, const uint64_t);
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
    void build(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>&);
    uint64_t streaming_query(const seqan3::bitpacked_sequence<seqan3::dna4>&, uint64_t&);
    void save(const std::filesystem::path&);
    void load(const std::filesystem::path&);
    void print_info();
};


class RSHash2
{
private:
    uint64_t k, m1, m_thres1, m2, m_thres2;
    uint64_t span1, span2;
    uint64_t kmermask, mmermask1, mmermask2;
    mixer_64 m_hasher1, m_hasher2;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r1, r2;
    bit_vector s1, s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select, s2_select;
    bits::compact_vector offsets1, offsets2;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    template<int level>
    inline uint64_t find_minimiser(const uint64_t, const uint64_t, size_t &, size_t &);
    template<int level>
    inline void update_minimiser(const uint64_t, const uint64_t, uint64_t&, size_t &, size_t &);
    template<int level>
    inline bool check(const uint64_t, const uint64_t, uint64_t*, const size_t, const size_t, const size_t, const size_t);
    template<int level>
    inline void fill_buffer(uint64_t *, uint64_t *, size_t, size_t, const uint64_t, const uint64_t);
    template<int level>
    inline bool check_minimiser_pos(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    template<int level>
    inline bool check_minimiser_pos2(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    template<int level>
    inline bool check_overlap(uint64_t, uint64_t, uint64_t &, uint64_t &);
    template<int level>
    inline bool lookup_buffer(uint64_t *, uint64_t *, const size_t, const uint64_t,  const uint64_t, uint64_t &, const size_t, const size_t, bool &, uint64_t &, uint64_t &);
    inline bool extend_in_text(uint64_t&, uint64_t, uint64_t, bool, const uint64_t, const uint64_t);
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
    int build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>&);
    uint64_t streaming_query(const seqan3::bitpacked_sequence<seqan3::dna4>&, uint64_t&);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSHash3
{
private:
    uint64_t k, m1, m_thres1, m2, m_thres2, m3, m_thres3;
    uint64_t span1, span2, span3;
    uint64_t kmermask, mmermask1, mmermask2, mmermask3;
    mixer_64 m_hasher1, m_hasher2, m_hasher3;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r1, r2, r3;
    bit_vector s1, s2, s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select, s2_select, s3_select;
    bits::compact_vector offsets1, offsets2, offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    template<int level>
    inline uint64_t find_minimiser(const uint64_t, const uint64_t, size_t &, size_t &);
    template<int level>
    inline void update_minimiser(const uint64_t, const uint64_t, uint64_t&, size_t &, size_t &);
    template<int level>
    inline bool check(const uint64_t, const uint64_t, uint64_t*, const size_t, const size_t, const size_t, const size_t);
    template<int level>
    inline void fill_buffer(uint64_t *, uint64_t *, size_t, size_t, const uint64_t, const uint64_t);
    template<int level>
    inline bool check_overlap(uint64_t, uint64_t, uint64_t &, uint64_t &);
    template<int level>
    inline bool check_minimiser_pos(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    template<int level>
    inline bool check_minimiser_pos2(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    template<int level>
    inline bool lookup_buffer(uint64_t *, uint64_t *, const size_t, const uint64_t,  const uint64_t, uint64_t &, const size_t, const size_t, bool &, uint64_t &, uint64_t &);
    inline bool extend_in_text(uint64_t&, uint64_t, uint64_t, bool, const uint64_t, const uint64_t);
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
    int build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>&);
    uint64_t streaming_query(const seqan3::bitpacked_sequence<seqan3::dna4>&, uint64_t&);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};