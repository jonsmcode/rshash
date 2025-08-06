#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>
#include <gtl/phmap.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

// namespace sdsl = seqan3::contrib::sdsl;
// const uint8_t m_thres1 = 20;
// const uint8_t m_thres2 = 20;

const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;

const size_t span = 31;


class RSIndex
{
private:
    uint8_t k, m1, m2, m_thres;
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
    sd_vector<> endpoints;
    rank_support_sd<> endpoints_rank;
    select_support_sd<> endpoints_select;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);

public:
    RSIndex();
    RSIndex(uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t m_thres);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class RSIndexComp
{
private:
    uint8_t k, m1, m2, m_thres;
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
    sd_vector<> endpoints;
    rank_support_sd<> endpoints_rank;
    select_support_sd<> endpoints_select;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    template<int level>
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);

public:
    RSIndexComp();
    RSIndexComp(uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t m_thres);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


// class RSIndex {
// protected:
//     uint8_t k, m1, m2;
//     bit_vector r1;
//     rank_support_v<1> r1_rank;
//     bit_vector s1;
//     sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
//     bit_vector s2;
//     sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
//     int_vector<0> offsets1;
//     int_vector<0> offsets2;
//     gtl::flat_hash_set<uint64_t> hashmap;
//     sd_vector<> endpoints;
//     rank_support_sd<> endpoints_rank;
//     select_support_sd<> endpoints_select;
//     seqan3::bitpacked_sequence<seqan3::dna4> text;
//     template<int level>
//     void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);

// public:
//     RSIndex() = default;
//     RSIndex(uint8_t k_, uint8_t m1_, uint8_t m2_)
//         : k(k_), m1(m1_), m2(m2_) {}
//     virtual ~RSIndex() = default;
//     uint8_t getk() { return k; }
//     virtual int build(const std::vector<std::vector<seqan3::dna4>>& data) = 0;
//     virtual int save(const std::filesystem::path& path) = 0;
//     virtual int load(const std::filesystem::path& path) = 0;
//     virtual uint64_t streaming_query(const std::vector<seqan3::dna4>& input) = 0;
// };


// class RSIndexComp : public RSIndex
// {
//     using RSIndex::RSIndex;

//     private:
//         sd_vector<> r2;
//         rank_support_sd<> r2_rank;

//     public:
//         int build(const std::vector<std::vector<seqan3::dna4>>&) override;
//         int save(const std::filesystem::path&) override;
//         int load(const std::filesystem::path&) override;
//         uint64_t streaming_query(const std::vector<seqan3::dna4>& input) override;
// };


// class RSIndexUncomp : public RSIndex
// {
//     using RSIndex::RSIndex;

//     private:
//         bit_vector r2;
//         rank_support_v<1> r2_rank;

//     public:
//         int build(const std::vector<std::vector<seqan3::dna4>>&) override;
//         int save(const std::filesystem::path&) override;
//         int load(const std::filesystem::path&) override;
//         uint64_t streaming_query(const std::vector<seqan3::dna4>& input) override;
// };
