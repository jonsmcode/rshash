#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>
#include <gtl/phmap.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

// namespace sdsl = seqan3::contrib::sdsl;

const size_t span = 31;

class RSIndex
{
private:
    uint8_t k, m;
    bit_vector r;
    rank_support_v<1> r_rank;
    bit_vector s;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s_select;
    int_vector<0> offsets;
    gtl::flat_hash_set<uint64_t> hashmap;
    sd_vector<> endpoints;
    rank_support_sd<> endpoints_rank;
    select_support_sd<> endpoints_select;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    void fill_buffer(std::vector<uint64_t>&, const uint64_t, size_t, size_t);


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


