#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

class Dictionary
{
private:
    uint8_t k, m;
    bit_vector r;
    bit_vector s;
    rank_support_v<1> r_rank;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> simple_select;
    int_vector<0> offset;
    int_vector<0> span;
    int offset_width;
    int span_width;

public:
    Dictionary();
    Dictionary(uint8_t const k, uint8_t const m);
    int build(const std::vector<seqan3::dna4>&);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&,
                        std::vector<uint64_t> &);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
    int save2(const std::filesystem::path&);
    int load2(const std::filesystem::path&);
};
