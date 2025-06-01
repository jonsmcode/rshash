#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

class Dictionary
{
private:
    uint8_t k, m;
    bit_vector r;
    bit_vector s;
    // rank_support_v<1> r_rank;
    // sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> simple_select;
    int_vector<0> offset;
    // int_vector<0> span;
    int offset_width;
    // int span_width;

public:
    Dictionary();
    Dictionary(uint8_t const k, uint8_t const m);
    uint8_t getk() { return k; }
    int build(const std::vector<seqan3::dna4>&);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&,
                        std::vector<uint64_t> &);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
    // int save2(const std::filesystem::path&);
    // int load2(const std::filesystem::path&);
};


// class LookupDictionary : public Dictionary
class LookupDictionary
{
private:
    uint8_t k, m;
    bit_vector r;
    bit_vector s;
    rank_support_v<1> r_rank;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> simple_select;
    int_vector<0> offset;
    int offset_width;

public:
    LookupDictionary();
    LookupDictionary(uint8_t const k, uint8_t const m);
    uint8_t getk() { return k; }
    int build(const std::vector<seqan3::dna4>&);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&,
                        std::vector<std::pair<uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
//     int load_comp(const std::filesystem::path &);
//     int save_comp(const std::filesystem::path &);
//     void print_statistics();
};


class UnitigsDictionary
{
private:
    uint8_t k, m;
    bit_vector r;
    bit_vector s;
    std::unordered_set<uint64_t> cbk;
    seqan3::contrib::sdsl::sd_vector<> endpoints;
    seqan3::contrib::sdsl::rank_support_sd<> endpoints_rank;
    seqan3::contrib::sdsl::select_support_sd<> endpoints_select;
    // int_vector<0> endpoints;
    rank_support_v<1> r_rank;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s_select;
    int_vector<0> offsets;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    void fill_buffer(std::vector<uint64_t> &, const uint64_t, size_t, size_t);
    void fill_buffer(std::vector<uint64_t>&, std::vector<uint64_t>&, std::vector<uint64_t>&, const uint64_t, size_t, size_t);


public:
    UnitigsDictionary();
    UnitigsDictionary(uint8_t const k, uint8_t const m);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class CompUnitigsDictionary
{
private:
    uint8_t k, m;
    seqan3::contrib::sdsl::sd_vector<> r;
    bit_vector s;
    std::unordered_set<uint64_t> cbk;
    sd_vector<> endpoints;
    rank_support_sd<> endpoints_rank;
    select_support_sd<> endpoints_select;
    rank_support_sd<> r_rank;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s_select;
    int_vector<0> offsets;
    seqan3::bitpacked_sequence<seqan3::dna4> text;
    void fill_buffer(std::vector<uint64_t> &, const uint64_t, size_t, size_t);
    void fill_buffer(std::vector<uint64_t>&, std::vector<uint64_t>&, std::vector<uint64_t>&, const uint64_t, size_t, size_t);

public:
    CompUnitigsDictionary();
    CompUnitigsDictionary(uint8_t const k, uint8_t const m);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


class LocateDictionary
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
    LocateDictionary();
    LocateDictionary(uint8_t const k, uint8_t const m);
    int build(const std::vector<seqan3::dna4>&);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&,
                        std::vector<uint64_t> &);
    int streaming_query(const std::vector<seqan3::dna4>&,
                        const std::vector<seqan3::dna4>&,
                        std::vector<std::pair<uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};
