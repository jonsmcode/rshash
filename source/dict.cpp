#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
// #include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>


using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;


class Dictionary
{
    uint8_t k, m;
    bit_vector r;
    bit_vector s;
    std::vector<bool> offset;
    std::vector<bool> span;
    int offset_width;
    int span_width;

    public:
        Dictionary(uint8_t const, uint8_t const);
        int build(std::vector<seqan3::dna4> const);
        int query(std::vector<seqan3::dna4> const);
        int streaming_query(std::vector<seqan3::dna4> const);
};

Dictionary::Dictionary(uint8_t const k, uint8_t const m) {
    this->k = k;
    this->m = m;
    r = bit_vector(1 << (2*m)); // bitvector for 4^m minimizers
}


int Dictionary::build(std::vector<seqan3::dna4> const text)
{
    // auto minimisers = text | seqan3::views::minimiser_hash(
    //     seqan3::shape{seqan3::ungapped{m}}, seqan3::window_size{k}, seqan3::seed{0});
    // uint64_t const seed = 0x8F'3F'73'B5'CF'1C'9A'DE;
    uint64_t const seed = 0;
    auto minimiser_hashs = text | seqan3::views::kmer_hash(seqan3::ungapped{m})
                                | std::views::transform([](uint64_t i) {return i ^ seed;})
                                | seqan3::views::minimiser(k);
    auto minimisers = minimiser_hashs | std::views::transform([seed](uint64_t i){ return i ^ seed;});
    // auto minimiser_hashs = text | seqan3::views::minimiser_hash(
    //     seqan3::shape{seqan3::ungapped{m}}, seqan3::window_size{k});
    // uint64_t seed = 0x8F'3F'73'B5'CF'1C'9A'DE;
    // auto minimizers = minimizer_hashs | std::views::transform([seed](uint64_t i){ return i ^ seed;});
    // auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{m}) | seqan3::views::minimiser(k);

    for(auto minimiser : minimisers) {
        r[minimiser] = 1;
        seqan3::debug_stream << minimiser << "\n";
    }

    sd_vector<> sdr(r);
    sd_vector<>::rank_1_type sdr_rank(&sdr);

    for (size_t i=0; i<r.size(); i++)
        seqan3::debug_stream  << r[i] << " ";
    std::cout << "\n";

    const int c = sdr_rank((1 << (2*m))-1);

    // uint8_t* count = new uint8_t[c] {0};
    uint8_t count[c];
    std::fill(count, count + c, 0);

    // todo: hashtable CB for more than 255 minimisers
    int n = 0;
    for(auto minimiser : minimisers) {
        count[sdr_rank(minimiser)]++; // + CB[minimiser] (if > 255)
        n++;
    }
    std::cout  << n <<  "\n";
    for (size_t i=0; i<c; i++)
        std::cout  << +count[i] <<  " ";
    std::cout << "\n";

    s = bit_vector(n, 0);
    s[0] = 1;
    int j = 0;
    for (size_t i=0; i<c-1; i++) {
        j += count[i]; // + CB
        s[j] = 1;
    }
    for (size_t i=0; i<n; i++)
        std::cout  << s[i] <<  " ";
    std::cout << "\n";

    offset_width = std::bit_width(text.size());
    span_width = std::bit_width(k);
    // const int width = offset_width + span_width;
    // offset = new std::vector<bool>[offset_width*n];
    offset.reserve(offset_width*n);
    span.reserve(span_width*n);

    // std::cout << offset.size() << "\n";
    // std::cout << span.size() << "\n";
    std::fill(count, count + c, 0);
    // zero CB
    // auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{m})
                           // | seqan3::views::minimiser(k);
                           // | std::views::transform([](uint64_t i) {return i ^ seed;});
    // auto forward_strand = text | seqan3::views::kmer_hash(seqan3::ungapped{m});
    auto forward_strand = seqan3::views::kmer_hash(seqan3::ungapped{m});
    seqan3::debug_stream << forward_strand << "\n";
    seqan3::detail::minimiser_view view = seqan3::detail::minimiser_view(forward_strand, k - m + 1);
    // seqan3::debug_stream << view.next_minimiser() << "\n";
    // for(auto kmer_hash : kmer_hashs) {
    // for(int i=0; i<n; i++) {
    //     while(!view.next_minimiser()) {}
    //     // seqan3::debug_stream << minimisers->minimiser_position_offset << "\n";
    // }
    // for(auto minimiser : minimisers) {
    //     seqan3::debug_stream << minimiser << "\n";
    //     // offset[offset_width*minimiser] = 0bi;
    //     // span[span_width*minimiser] = 0bj;
    //     i++;
    // }
    // auto minimisers = minimiser_hashs | std::views::transform([seed](uint64_t i){ return i ^ seed;});


    return 0;
}


// int query(std::vector<seqan3::dna4> const &query) {
//     return 0;
// }


// int streaming_query(std::vector<seqan3::dna4> const &query) {
//     return 0;
// }


int main()
{
    // const int m = 2; // minimiser
    // const int k = 5; // k-mer

    // auto reference_stream = seqan3::sequence_file_input{reference_file};

    // // read into memory
    // std::vector<std::vector<seqan3::dna4>> text;
    // for (auto& record : reference_stream) {
    //     text.push_back(record.sequence());
    // }

    // std::vector<seqan3::dna4> text{"TCATCAGTAGCTACATTACG"_dna4};
    std::vector<seqan3::dna4> text{"TCATCAGTAGCTACATTACG"_dna4};
    // std::vector<seqan3::dna4> text{"CATTAC"_dna4};

    Dictionary dict(5, 2);
    dict.build(text);

    // std::vector<seqan3::dna4> query{"GTAGCTA"_dna4};

    
    return 0;
}