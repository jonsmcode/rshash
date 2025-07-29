#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

namespace sdsl = seqan3::contrib::sdsl;


class UnitigsDictionaryHash2
{
private:
    uint8_t k, m;
    bit_vector r1;
    bit_vector r2;
    bit_vector s1;
    bit_vector s2;
    seqan3::contrib::sdsl::sd_vector<> endpoints;
    seqan3::contrib::sdsl::rank_support_sd<> endpoints_rank;
    seqan3::contrib::sdsl::select_support_sd<> endpoints_select;
    rank_support_v<1> r1_rank;
    rank_support_v<1> r2_rank;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
    int_vector<0> offsets1;
    int_vector<0> offsets2;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    void fill_buffer1(std::vector<uint64_t> &, const uint64_t, size_t, size_t);
    void fill_buffer2(std::vector<uint64_t> &, const uint64_t, size_t, size_t);


public:
    UnitigsDictionaryHash2();
    UnitigsDictionaryHash2(uint8_t const k, uint8_t const m);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


