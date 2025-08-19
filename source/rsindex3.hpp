#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>
#include <gtl/phmap.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

// namespace sdsl = seqan3::contrib::sdsl;

const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed3 = 0xE5'9A'38'5F'03'76'C9'F6;

const size_t span = 31;


class RSIndex
{
private:
    uint8_t k, m1, m2, m3, m_thres1, m_thres2;
    uint16_t m_thres3;
    bit_vector r1;
    rank_support_v<1> r1_rank;
    bit_vector r2;
    rank_support_v<1> r2_rank;
    bit_vector r3;
    rank_support_v<1> r3_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    bit_vector s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s3_select;
    int_vector<0> offsets1;
    int_vector<0> offsets2;
    int_vector<0> offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sd_vector<> endpoints;
    rank_support_sd<> endpoints_rank;
    select_support_sd<> endpoints_select;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);


public:
    RSIndex();
    RSIndex(uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m3,
        uint8_t const m_thres1, uint8_t const m_thres2, uint16_t const m_thres3);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSIndexComp
{
private:
    uint8_t k, m1, m2, m3, m_thres1, m_thres2;
    uint16_t m_thres3;
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
    int_vector<0> offsets1;
    int_vector<0> offsets2;
    int_vector<0> offsets3;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    // sd_vector<> sequences;
    // rank_support_sd<> endpoints_rank;
    // select_support_sd<> endpoints_select;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);


public:
    RSIndexComp();
    RSIndexComp(uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m3,
        uint8_t const m_thres1, uint8_t const m_thres2, uint16_t const m_thres3);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


