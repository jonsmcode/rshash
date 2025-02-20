#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

class Dictionary
{
    uint8_t k, m;
    bit_vector r;
    bit_vector s;
    sd_vector<> sdr;
    sd_vector<> sds;
    sd_vector<>::rank_1_type r_rank;
    sd_vector<>::select_1_type s_select;
    int_vector<0> offset;
    int_vector<0> span;
    int offset_width;
    int span_width;

    public:
        Dictionary(uint8_t const k, uint8_t const m);
        int build(std::vector<seqan3::dna4>&);
        int streaming_query(std::vector<seqan3::dna4>&,
            std::vector<seqan3::dna4>&,
            std::vector<uint64_t> &);
};
