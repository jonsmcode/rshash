#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>
#include <gtl/phmap.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;


const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;


class RSIndex
{
private:
    uint8_t k, m1, m2, m_thres1, m_thres2;
    size_t span;
    bit_vector r1;
    rank_support_v<1> r1_rank;
    bit_vector r2;
    rank_support_v<1> r2_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    int_vector<0> offsets1;
    int_vector<0> offsets2;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);

public:
    RSIndex();
    RSIndex(uint8_t const k, uint8_t const m1, uint8_t const m2,
        const uint8_t m_thres1, const uint8_t m_thres2, size_t const span);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSIndexComp
{
private:
    uint8_t k, m1, m2, m_thres1, m_thres2;
    size_t span;
    bit_vector r1;
    rank_support_v<1> r1_rank;
    sd_vector<> r2;
    rank_support_sd<> r2_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    int_vector<0> offsets1;
    int_vector<0> offsets2;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);

public:
    RSIndexComp();
    RSIndexComp(uint8_t const k, uint8_t const m1, uint8_t const m2,
        const uint8_t m_thres1, const uint8_t m_thres2, size_t const span);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};



class RSIndexComp2
{
private:
    uint8_t k, m1, m2, m_thres1, m_thres2;
    size_t span;
    sd_vector<> r1;
    rank_support_sd<> r1_rank;
    sd_vector<> r2;
    rank_support_sd<> r2_rank;
    bit_vector s1;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    bit_vector s2;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    int_vector<0> offsets1;
    int_vector<0> offsets2;
    gtl::flat_hash_set<uint64_t> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    bit_vector sequences;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);

public:
    RSIndexComp2();
    RSIndexComp2(uint8_t const k, uint8_t const m1, uint8_t const m2,
        const uint8_t m_thres1, const uint8_t m_thres2, size_t const span);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, uint64_t &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};

