#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/Rank9Sel.hpp>
#include <sux/bits/EliasFano.hpp>
#include <gtl/phmap.hpp>

#include "compact_vector.hpp"

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;


// const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed1 = 1;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed3 = 0xE5'9A'38'5F'03'76'C9'F6;


class RSIndex
{
private:
    uint8_t k, m1, m2, m3, m_thres1, m_thres2;
    uint16_t m_thres3;
    size_t span;
    bit_vector r1;
    rank_support_v<1> r1_rank;
    bit_vector r2;
    rank_support_v<1> r2_rank;
    bit_vector r3;
    rank_support_v<1> r3_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    // sux::bits::Rank9Sel<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    // sux::bits::Rank9Sel<sux::util::AllocType::MALLOC> s2_select;
    bit_vector s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s3_select;
    // sux::bits::Rank9Sel<sux::util::AllocType::MALLOC> s3_select;
    pthash::compact_vector offsets1;
    pthash::compact_vector offsets2;
    pthash::compact_vector offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);


public:
    RSIndex();
    RSIndex(uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m3,
        uint8_t const m_thres1, uint8_t const m_thres2, uint16_t const m_thres3, size_t const span);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    uint64_t access(const uint64_t, const size_t);
    std::vector<uint64_t> rand_text_kmers(const uint64_t);
    uint64_t lookup(const std::vector<uint64_t>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
    void stats();
};


class RSIndexComp
{
private:
    uint8_t k, m1, m2, m3, m_thres1, m_thres2;
    uint16_t m_thres3;
    size_t span;
    bit_vector r1;
    rank_support_v<1> r1_rank;
    sd_vector<> r2;
    rank_support_sd<> r2_rank;
    sd_vector<> r3;
    rank_support_sd<> r3_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    bit_vector s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s3_select;
    pthash::compact_vector offsets1;
    pthash::compact_vector offsets2;
    pthash::compact_vector offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t);
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);


public:
    RSIndexComp();
    RSIndexComp(uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m3,
        uint8_t const m_thres1, uint8_t const m_thres2, uint16_t const m_thres3, size_t const span);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    uint64_t access(const uint64_t, const size_t);
    std::vector<uint64_t> rand_text_kmers(const uint64_t);
    uint64_t lookup(const std::vector<uint64_t>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSHash1
{
private:
    uint8_t k, m1, m_thres1;
    size_t span;
    sd_vector<> r1;
    rank_support_sd<> r1_rank;
    // rrr_vector<15> r1;
    // rank_support_rrr<1, 15> r1_rank;
    // bit_vector r1;
    // rank_support_v<1> r1_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    pthash::compact_vector offsets1;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t);
    inline void fill_buffer(std::vector<uint64_t>&, std::vector<size_t> &, std::vector<size_t> &, std::vector<size_t> &, std::vector<size_t> &, size_t, size_t);
    inline bool lookup_serial(std::vector<uint64_t> &, std::vector<size_t> &, std::vector<size_t> &, std::vector<size_t> &, std::vector<size_t> &, const uint64_t,  const uint64_t, size_t &, size_t &, bool &, size_t &, size_t &, size_t &);
    inline bool extend_in_text(size_t&, size_t, size_t, bool, const uint64_t, const uint64_t, size_t&, size_t&);


public:
    RSHash1();
    RSHash1(uint8_t const k, uint8_t const m1, uint8_t const m_thres1, size_t const span);
    uint8_t getk() { return k; }
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    std::vector<uint64_t> rand_text_kmers(const uint64_t);
    uint64_t access(const uint64_t, const size_t);
    uint64_t lookup(const std::vector<uint64_t>&, bool verbose);
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&, size_t &, size_t &);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSHash2
{
private:
    uint8_t k, m1, m_thres1, m2, m_thres2;
    size_t span;
    sd_vector<> r1;
    rank_support_sd<> r1_rank;
    sd_vector<> r2;
    rank_support_sd<> r2_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    pthash::compact_vector offsets1;
    pthash::compact_vector offsets2;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t);
    template<int level>
    inline void fill_buffer(std::vector<uint64_t>&, size_t, size_t);
    // inline void fill_text_buffer(std::vector<uint64_t>&, size_t, size_t, size_t&, bool, size_t&, size_t&);
    // inline void fill_minimiser_buffer(std::vector<uint64_t>&, std::vector<size_t>&, size_t, size_t);


public:
    RSHash2();
    RSHash2(uint8_t const k, uint8_t const m1, uint8_t const m_thres1,
                 uint8_t const m2, uint8_t const m_thres2, size_t const span);
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
    uint8_t k, m1, m2, m3, m_thres1, m_thres2;
    uint16_t m_thres3;
    size_t span;
    sd_vector<> r1;
    rank_support_sd<> r1_rank;
    sd_vector<> r2;
    rank_support_sd<> r2_rank;
    sd_vector<> r3;
    rank_support_sd<> r3_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    bit_vector s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s3_select;
    // int_vector<0> offsets1;
    // int_vector<0> offsets2;
    // int_vector<0> offsets3;
    pthash::compact_vector offsets1;
    pthash::compact_vector offsets2;
    pthash::compact_vector offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t, double &, double &, double &);
    template<int level>
    inline bool check(const size_t, const size_t, const uint64_t, const uint64_t, const uint64_t);
    template<int level>
    inline void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);


public:
    RSHash3();
    RSHash3(uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m3,
        uint8_t const m_thres1, uint8_t const m_thres2, uint16_t const m_thres3, size_t const span);
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


