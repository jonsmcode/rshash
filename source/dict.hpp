#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <sux/bits/EliasFano.hpp>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;

// namespace sdsl = seqan3::contrib::sdsl;

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


class UnitigsDictionaryHash
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
    UnitigsDictionaryHash();
    UnitigsDictionaryHash(uint8_t const k, uint8_t const m);
    uint8_t getk() { return k; }
    int build(const std::vector<std::vector<seqan3::dna4>>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&);
    uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
};


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


// class UnitigsDictionaryHash3
// {
// private:
//     uint8_t k, m;
//     uint64_t M;
//     bit_vector r1;
//     bit_vector r2;
//     bit_vector s1;
//     bit_vector s2;
//     seqan3::contrib::sdsl::sd_vector<> endpoints;
//     seqan3::contrib::sdsl::rank_support_sd<> endpoints_rank;
//     seqan3::contrib::sdsl::select_support_sd<> endpoints_select;
//     rank_support_v<1> r1_rank;
//     rank_support_v<1> r2_rank;
//     sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select;
//     sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s2_select;
//     int_vector<0> offsets1;
//     int_vector<0> offsets2;
//     seqan3::bitpacked_sequence<seqan3::dna4> text;
//     size_t get_endpoints(const std::vector<std::vector<seqan3::dna4>> &input);
//     template <typename ViewType>
//     void compute_level1(const std::vector<std::vector<seqan3::dna4>> &input,
//         const ViewType& view, const uint64_t text_length);
//     template <typename ViewType>
//     void mark_minimisers(const std::vector<std::vector<seqan3::dna4>> &input,
//         const ViewType& view, sdsl::bit_vector &r);
//     template <typename ViewType>
//     void mark_unfreq_minimisers(const std::vector<std::vector<seqan3::dna4>> &input,
//         const ViewType& view, sdsl::bit_vector &r, uint8_t* count, sdsl::rank_support_v<1> &r_rank, uint8_t const m_thres);
//     template <typename ViewType>
//     uint8_t* count_minimisers(const std::vector<std::vector<seqan3::dna4>> &input,
//         const ViewType& view, sdsl::rank_support_v<1> &r_rank, const uint8_t m_thres);
//     template <typename ViewType>
//     uint8_t* update_counts(const std::vector<std::vector<seqan3::dna4>> &input,
//         const ViewType& view, uint8_t* count, sdsl::bit_vector &r, sdsl::bit_vector &rnew,
//         sdsl::rank_support_v<1> &r_rank, sdsl::rank_support_v<1> &rnew_rank);
//     uint64_t build_S1(uint8_t* count, sdsl::rank_support_v<1> &r_rank);
//     template <typename ViewType>
//     void fill_offsets(const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view,
//     uint8_t* count, sdsl::bit_vector &r, sdsl::rank_support_v<1> &r_rank,
//     sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s_select, int_vector<0> &offsets,
//     uint64_t n, uint64_t N);
//     template <typename ViewType>
//     std::vector<std::vector<seqan3::dna4>> get_frequent_skmers(
//         const std::vector<std::vector<seqan3::dna4>> &input,
//         const ViewType& view, sdsl::bit_vector &r);
//     template <int> void fill_buffer(std::vector<uint64_t> &, const uint64_t, size_t, size_t);


// public:
//     UnitigsDictionaryHash3();
//     UnitigsDictionaryHash3(uint8_t const k, uint8_t const m);
//     uint8_t getk() { return k; }
//     int build(const std::vector<std::vector<seqan3::dna4>>&);
//     uint64_t streaming_query(const std::vector<seqan3::dna4>&);
//     uint64_t streaming_query(const std::vector<seqan3::dna4>&, std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &);
//     int save(const std::filesystem::path&);
//     int load(const std::filesystem::path&);
// };


class UnitigsDictionarySIMD
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
    UnitigsDictionarySIMD();
    UnitigsDictionarySIMD(uint8_t const k, uint8_t const m);
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
