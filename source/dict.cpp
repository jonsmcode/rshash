#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/minimiser_hash_and_positions.hpp>
#include <seqan3/search/views/minimiser_and_window_hash.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>

// #include <boost/dynamic_bitset.hpp>
// #include <bitset>

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;


// class Hash
// {

//     auto urange_size = std::ranges::distance(it_start, it_end);
//         auto step = (shape_.size() > urange_size + 1) ? 0 : urange_size - shape_.size() + 1;
//         text_left = std::ranges::next(it_start, step, it_end);

//         // shape size = 3
//         // Text:      1 2 3 4 5 6 7 8 9
//         // text_left: ^
//         // text_right:    ^
//         // distance(text_left, text_right) = 2
//         if (shape_.size() <= std::ranges::distance(text_left, it_end) + 1)
//         {
//             roll_factor = pow(sigma, static_cast<size_t>(std::ranges::size(shape_) - 1));
//             hash_full();
//         }

//         text_right = it_end;

//     value_type operator*() const noexcept
//     {
//         return hash_value + to_rank(*text_right);
//     }

//     void hash_forward()
//     {
//         if (shape_.all())
//         {
//             hash_roll_forward();
//         }
//         else
//         {
//             std::ranges::advance(text_left, 1);
//             hash_full();
//         }
//     }

//     void hash_full() {

//     }

//     void hash_roll_forward()
//     {
//         hash_value -= to_rank(*(text_left)) * roll_factor;
//         hash_value += to_rank(*(text_right));
//         hash_value *= sigma;

//         std::ranges::advance(text_left, 1);
//         std::ranges::advance(text_right, 1);
//     }
// }



seqan3::dna4_vector kmer_to_string(uint64_t kmer, size_t const kmer_size)
{
    seqan3::dna4_vector result(kmer_size);
    for (size_t i = 0; i < kmer_size; ++i)
    {
        result[kmer_size - 1 - i].assign_rank(kmer & 0b11);
        kmer >>= 2;
    }
    return result;
}


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
        Dictionary(uint8_t const, uint8_t const);
        int build(std::vector<seqan3::dna4>&);
        int streaming_query(std::vector<seqan3::dna4>&, std::vector<seqan3::dna4>&);
};

Dictionary::Dictionary(uint8_t const k, uint8_t const m) {
    this->k = k;
    this->m = m;
    r = bit_vector(1 << (2*m)); // bitvector for 4^m minimizers
}


int Dictionary::build(std::vector<seqan3::dna4> &text)
{
    auto view = bsc::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k});

    for(auto && minimiser : text | view) {
        r[minimiser.minimiser_value] = 1;
    }

    sdr = sd_vector<>(r);
    r_rank = sd_vector<>::rank_1_type(&sdr);

    for (size_t i=0; i<r.size(); i++)
        seqan3::debug_stream  << r[i] << " ";
    std::cout << "\n";

    const int c = r_rank((1 << (2*m))-1);

    uint8_t count[c];
    std::fill(count, count + c, 0);
    // zero CB

    // todo: hashtable CB for more than 255 minimisers
    // todo: implement minmizer s.t. minmizer range at most window length
    int n = 0;
    for(auto && minimiser : text | view) {
        int r = r_rank(minimiser.minimiser_value);
        int mo = minimiser.occurrences;
        while(mo > k) {
            mo -= k;
            count[r]++; // + CB[minimiser] (if > 255)
            n++;
        }
        count[r]++; // + CB[minimiser] (if > 255)
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

    sds = sd_vector<>(s);
    s_select = sd_vector<>::select_1_type(&sds);

    size_t offset_width = std::bit_width(text.size());
    size_t span_width = std::bit_width(k);
    // offset.width(offset_width + span_width);
    offset.width(offset_width);
    offset.resize(n);
    span.width(span_width);
    span.resize(n);
    // std::cout << "offset width " << offset_width << ", span width " << span_width << std::endl;
    // std::cout << "total width " << (int) offset.width() << ", total bits " << offset.bit_size() << std::endl;

    std::fill(count, count + c, 0);
    // zero CB

    for (auto && minimiser : text | view) {
        int r = r_rank(minimiser.minimiser_value);
        int s = s_select(minimiser.minimiser_value);
        // seqan3::debug_stream << minimiser.minimiser_value << kmer_to_string(minimiser.minimiser_value, m) << ',' << minimiser.range_position
                             // << ',' << minimiser.occurrences << '\n';
        int mo = minimiser.occurrences;
        int j = 0;
        while(mo > k) {
            mo -= k;
            // offset[s + count[r]] = (minimiser.range_position + j*k << offset_width) | k;
            offset[s + count[r]] = minimiser.range_position + j*k;
            span[s + count[r]] = k;
            j++;
            count[r]++;
        }
        // offset[s + count[r]] = (minimiser.range_position + j*k << offset_width) | mo;
        offset[s + count[r]] = minimiser.range_position + j*k;
        span[s + count[r]] = mo;
        count[r]++;
    }
    // uint64_t mask = (1 << span_width) - 1;
    // for (size_t i=0; i<n; i++)
    //     std::cout << (offset[i] >> offset_width) << " " << (offset[i] & mask) << "\n";
    // std::cout << "\n";
    for (size_t i=0; i<n; i++)
        std::cout << offset[i] << "\n";
    std::cout << "\n";
    for (size_t i=0; i<n; i++)
        std::cout << span[i] << "\n";
    std::cout << "\n";

    return 0;
}


int Dictionary::streaming_query(
    std::vector<seqan3::dna4> &text,
    std::vector<seqan3::dna4> &query)
{
    auto view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});
    // uint64_t mask = static_cast<uint64_t>(1 << span_width - 1);

    // todo: parallelize queries
    for(auto && minimiser : query | view) {
        size_t r = r_rank(minimiser.minimiser_value)+1;
        // todo: fast select bitvector implementation
        // todo: test compressed size array
        size_t p = s_select(r);
        size_t q = s_select(r+1);

        size_t b = q - p;
        // seqan3::debug_stream << kmer_to_string(minimiser.minimiser_value, m) << ','
        //                      << kmer_to_string(minimiser.window_value, k) << '\n';
        // uint64_t o = (static_cast<uint64_t>(offset[p]) >> offset_width);
        // uint64_t s = static_cast<uint64_t>(offset[p] & mask);
        // std::cout << std::bitset<10>(offset[p]) << "\n";
        // std::cout << std::bitset<6>(o) << "\n";
        // std::cout << std::bitset<4>(s) << "\n";
        // std::cout << o << " " << s << "\n";

        std::vector<size_t> const buffer;
        for(int i = 0; i < b; i++) {
            // std::cout << offset[p+i] << "\n";
            // std::cout << span[p+i] << "\n";
            // for(int j = 0; j < span[p+i]; j++) {
            //     buffer.push_back(hash());
            // }
            // for(auto &&minimiser = text[offset[p]] | view; minimiser != text[offset[p] + span[p+i]] | view; ++minimiser) {
            //     seqan3::debug_stream << kmer_to_string(minimiser.minimiser_value, m) << ','
            //                  << kmer_to_string(minimiser.window_value, k) << '\n';
            // }
            for (auto && minimiser : text | view) {
                seqan3::debug_stream << kmer_to_string(minimiser.minimiser_value, m) << ','
                             << kmer_to_string(minimiser.window_value, k) << '\n';
            }
        }
        

    }

    return 0;
}

// int Dictionary::hash_buf(std::vector<size_t> const &buffer, std::vector<seqan3::dna4> const &query)
// {
//     for(int i = 0; i < b; i++) {

//     }
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
    std::vector<seqan3::dna4> text{"TCATCAGTAGCTACATTACG"_dna4};

    Dictionary dict(5, 2);
    dict.build(text);

    std::vector<seqan3::dna4> query{"GTAGCTAGCTACA"_dna4};
    dict.streaming_query(text, query);

    
    return 0;
}