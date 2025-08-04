#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>
#include <gtl/phmap.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

// namespace sdsl = seqan3::contrib::sdsl;

const size_t span = 100;

class RSIndex
{
private:
    uint8_t k, m;
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
    void fill_buffer1(std::vector<uint64_t>&, const uint64_t, size_t, size_t);
    void fill_buffer2(std::vector<uint64_t>&, const uint64_t, size_t, size_t);


public:
    RSIndex();
    RSIndex(uint8_t const k, uint8_t const m);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


