#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>

#include "../source/minimiser_rev_xor_views.hpp"
#include "../source/minimiser_rev_hash_views.hpp"


const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 1;

const uint8_t m_thres1 = 10;
const uint8_t m_thres2 = 10;

const uint8_t k = 31;
const uint8_t m = 16;

const uint64_t M = 1ULL << (m+m);

const size_t span = 100;

namespace sdsl = seqan3::contrib::sdsl;


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


uint64_t load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    uint64_t N = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
        if(N >= 1000000000)
            break;
    }
    return N;
}


template <typename ViewType>
inline uint8_t* count_minimisers(const ViewType& view, const uint8_t m_thres, const std::vector<std::vector<seqan3::dna4>> &input,
                        sdsl::bit_vector &r, sdsl::rank_support_v<1> &r_rank)
{
    for(auto & record : input) {
        for(auto && minimiser : record | view)
            r[minimiser.minimiser_value] = 1;
    }

    r_rank = sdsl::rank_support_v<1>(&r);
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    auto update_count = [&](uint64_t minimiser_value, size_t occurrences) {
        size_t i = r_rank(minimiser_value);
        count[i] += occurrences/span + 1;
        if(count[i] > m_thres)
            count[i] = m_thres;
    };

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view)
            update_count(minimiser.minimiser_value, minimiser.occurrences);
    }

    return count;
}


template <typename ViewType>
inline uint8_t* build_r_and_count(const ViewType& view, const uint8_t m_thres, const std::vector<std::vector<seqan3::dna4>> &input,
                            sdsl::bit_vector &rnew, sdsl::rank_support_v<1> &rnew_rank,
                            uint8_t* count, sdsl::rank_support_v<1> &r_rank)
{
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            if(count[i] < m_thres)
                rnew[minimiser.minimiser_value] = 1;
        }
    }

    rnew_rank = sdsl::rank_support_v<1>(&rnew);

    size_t c = rnew_rank(M);
    uint8_t* new_count = new uint8_t[c];
    std::memset(new_count, 0, c*sizeof(uint8_t));

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            if(rnew[minimiser.minimiser_value])
                new_count[rnew_rank(minimiser.minimiser_value)] = count[r_rank(minimiser.minimiser_value)];
        }
    }

    return new_count;
}


template <typename ViewType>
inline uint8_t* compute_level(const ViewType& view, const uint8_t m_thres, const std::vector<std::vector<seqan3::dna4>> &input,
                        sdsl::bit_vector &r,
                        sdsl::rank_support_v<1> &r_rank)
{
    sdsl::bit_vector rtmp = sdsl::bit_vector(M, 0);
    sdsl::rank_support_v<1> rtmp_rank;

    std::cout << "count minimisers and fill temporary bitvector R and count array...\n";
    uint8_t* count_tmp = count_minimisers(view, m_thres, input, rtmp, rtmp_rank);
    std::cout << "build final bitvector and count array...\n";
    uint8_t* count = build_r_and_count(view, m_thres, input, r, r_rank, count_tmp, rtmp_rank);

    delete[] count_tmp;
    // delete rtmp

    return count;
}


template <typename ViewType>
inline std::vector<std::vector<seqan3::dna4>> get_frequent_skmers(
    const ViewType& view, const std::vector<std::vector<seqan3::dna4>> &input, sdsl::bit_vector &r)
{
    std::vector<std::vector<seqan3::dna4>> freq_sequences;
    for(auto & sequence : input) {
        size_t start_position;
        bool freq;
        bool current_freq;
        for(auto && minimiser : sequence | view) {
            freq = !r[minimiser.minimiser_value];
            if(freq)
                start_position = 0;
            break;
        }
        for(auto && minimiser : sequence | view) {
            current_freq = !r[minimiser.minimiser_value];

            if(!freq && current_freq)
                start_position = minimiser.range_position;
            if(freq && !current_freq) {
                std::vector<seqan3::dna4> freq_sequence;
                for(size_t i=start_position; i < minimiser.range_position+k-1; i++)
                    freq_sequence.push_back(sequence[i]);
                freq_sequences.push_back(freq_sequence);
            }

            freq = current_freq;
        }
    }

    return freq_sequences;
}


void stats(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m, .window_size = k, .seed=seed1});
    auto view2 = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed2});

    std::cout << "counting (super-)kmers...\n";
    size_t N = 0;
    uint64_t skmers = 0;
    uint64_t kmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            size_t o = minimiser.occurrences;
            const uint8_t w = o/span + 1;
            kmers += o;
            skmers += w;
        }
        N += sequence.size();
    }

    std::cout << "compute level 1...\n";
    sdsl::bit_vector r1 = sdsl::bit_vector(M, 0);
    sdsl::rank_support_v<1> r1_rank;
    uint8_t* count1 = compute_level(view1, m_thres1, input, r1, r1_rank);

    std::cout << "extract uncovered sequence parts...\n";
    std::vector<std::vector<seqan3::dna4>> remaining_sequences1;
    remaining_sequences1 = get_frequent_skmers(view1, input, r1);

    size_t len_rem_seqs = 0;
    for(auto & sequence : remaining_sequences1) {
        len_rem_seqs += sequence.size();
    }
    std::cout << "remaining superkmers " << remaining_sequences1.size() << " (" << (double)remaining_sequences1.size()/skmers*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/N*100 << "%)\n";

    std::cout << "compute level 2...\n";
    sdsl::bit_vector r2 = sdsl::bit_vector(M, 0);
    sdsl::rank_support_v<1> r2_rank;
    uint8_t* count2 = compute_level(view2, m_thres2, remaining_sequences1, r2, r2_rank);

    std::cout << "extract uncovered sequence parts...\n";
    std::vector<std::vector<seqan3::dna4>> remaining_sequences2;
    remaining_sequences2 = get_frequent_skmers(view2, remaining_sequences1, r2);

    len_rem_seqs = 0;
    for(auto & sequence : remaining_sequences2) {
        len_rem_seqs += sequence.size();
    }
    std::cout << "remaining superkmers " << remaining_sequences2.size() << " (" << (double)remaining_sequences2.size()/skmers*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/N*100 << "%)\n";

    std::cout << "filling HT...\n";
    std::unordered_set<uint64_t> freq_kmers;
    std::unordered_set<uint64_t> freq_minimzer;

    auto view3 = srindex::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k, .seed=seed2});
    for(auto & sequence : remaining_sequences2) {
        for(auto && minimiser : sequence | view3) {
            freq_minimzer.insert(minimiser.minimiser_value);
            freq_kmers.insert(minimiser.window_value);
        }
    }

    std::cout << "computing distributions...\n";

    size_t c1 = r1_rank(M);
    uint64_t counter1[m_thres1] = {0};
    for(size_t i=0; i < c1; i++) {
        for(size_t j=0; j < m_thres1; j++) {
            if(count1[i] == j) {
                counter1[j]++;
                break;
            }
        }
    }
    uint64_t n1 = 0;
    for(size_t i=0; i < c1; i++)
        n1 += count1[i];

    size_t c2 = r2_rank(M);
    uint64_t counter2[m_thres2] = {0};
    for(size_t i=0; i < c2; i++) {
        for(size_t j=0; j < m_thres2; j++) {
            if(count2[i] == j) {
                counter2[j]++;
                break;
            }
        }
    }
    uint64_t n2 = 0;
    for(size_t i=0; i < c2; i++)
        n2 += count2[i];
    

    uint64_t cum = 0;
    uint64_t cum_skmers = 0;
    std::cout << "\nminimiser 1 distribution:\n";
    for(uint64_t j=0; j < m_thres1; j++) {
        cum += counter1[j];
        cum_skmers += (j)*counter1[j];
        std::cout << "occurrences " << (j) << ": " << counter1[j] << " " << (double) counter1[j]/c1*100 << "%  cum: " << cum << " " << (double) cum/c1*100 << "% covering " << cum_skmers << ' ' << (double) cum_skmers/skmers*100 << "% superkmers\n";
    }
    std::cout << "avg superkmers1: " << (double) n1/c1 <<  '\n';
    // std::cout << "minimisers going to level 2: " << c-c1 << "  " << (double) (c-c1)/c*100 << "% to cover " << (double) (n-n1)/n*100 << "% superkmers\n";
    
    cum = 0;
    std::cout << "\nminimiser 2 distribution:\n";
    for(int j=0; j < m_thres2; j++) {
        cum += counter2[j];
        cum_skmers += (j)*counter2[j];
        std::cout << "occurrences " << (j) << ": " << counter2[j] << " " << (double) counter2[j]/c2*100 << "%  cum: " << cum << " " << (double) cum/c2*100 << "% covering " << cum_skmers << ' ' << (double) cum_skmers/skmers*100 << "% superkmers\n";
    }
    std::cout << "avg superkmers2: " << (double) n2/c2 <<  '\n';
    // std::cout << "minimisers going to HT: " << c1-c2 << " " << (double) freq_minimzer.size()/(n2)*100 << "%\n";

    // std::cout << "superkmers to cover by HT: " << 100.0 - (double) cum_skmers/(n1+n2+n3)*100 << "%\n";
    // std::cout << "freq minimisers left: " << freq_minimzer.size() << " " << (double) freq_minimzer.size()/n*100 << "%\n";
    // std::cout << "kmers in HT: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";

    std::cout << "\n====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no superkmers: " << skmers << "\n";
    std::cout << "no minimiser1: " << n1 << "\n";
    std::cout << "no distinct minimiser1: " << c1 << "\n";
    std::cout << "no minimiser2: " << n2 << "\n";
    std::cout << "no distinct minimiser2: " << c2 << "\n";
    std::cout << "no minimiser HT: " << freq_minimzer.size() << " " << (double)freq_minimzer.size()/skmers*100 << "%\n";
    std::cout << "no kmers HT: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";
    std::cout << "density r1: " << (double) c1/M*100 << "%\n";
    std::cout << "density r2: " << (double) c2/M*100 << "%\n";

    delete[] count1;
    delete[] count2;
    // delete r1
    // delete r2
}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}
