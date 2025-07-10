#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "../source/minimiser_rev_hash_views.hpp"
#include "../source/minimiser_rev_xor_views.hpp"
#include "../source/minimiser_rev_hash_views2.hpp"


const uint8_t m_thres = 10;
const uint8_t m_thres2 = 10;
const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
// const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed2 = 1;


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


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    uint64_t N = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
        // if(N >= 1000000000)
        //     break;
        // if(N >= 100000000)
        //     break;
    }
}


void stats(const std::vector<std::vector<seqan3::dna4>> &input)
{
    const uint8_t k = 31;
    size_t N = 0;
    for(auto & record : input) {
        N += record.size();
    }
    // uint8_t m = std::log(N)/std::log(4);
    // if(m > 16)
    //     m = 16;
    // uint8_t m = 15;
    uint8_t m = 16;
    std::cout << "m: " << +m << '\n';

    // const size_t span = k-m;
    const size_t span = 100;

    // auto view1 = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed1});
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m, .window_size = k, .seed=seed1});
    auto view2 = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed2});
    auto view12 = srindex::views::two_minimisers_and_window_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});

    const uint64_t M = 1ULL << (m+m);

    std::cout << "scanning minimizers...\n";
    seqan3::contrib::sdsl::bit_vector r = seqan3::contrib::sdsl::bit_vector(M, 0);
    
    for(auto & record : input) {
        for(auto && minimiser : record | view1)
            r[minimiser.minimiser_value] = 1;
    }
    seqan3::contrib::sdsl::rank_support_v<1> r_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r);

    size_t c = r_rank(M);
    std::cout << "found " << c << " distinct minimisers \n";
    std::cout << "counting minimisers...\n";
    
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    uint64_t kmers = 0;
    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;

            const uint64_t w = o/span + 1;
            if(count[i] + w >= m_thres)
                count[i] = m_thres;
            else
                count[i] += w;
            kmers += o;
            n += w;
        }
    }

    std::cout << "filling R1...\n";
    seqan3::contrib::sdsl::bit_vector r1 = seqan3::contrib::sdsl::bit_vector(M, 0);
    seqan3::contrib::sdsl::bit_vector r2_ = seqan3::contrib::sdsl::bit_vector(M, 0);

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view12) {
            size_t i = r_rank(minimiser.minimiser1_value);
            if(count[i] < m_thres)
                r1[minimiser.minimiser1_value] = 1;
            else
                r2_[minimiser.minimiser2_value] = 1;
        }
    }

    delete[] count;

    seqan3::contrib::sdsl::rank_support_v<1> r1_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r1);
    seqan3::contrib::sdsl::rank_support_v<1> r2__rank = seqan3::contrib::sdsl::rank_support_v<1>(&r2_);

    std::cout << "count minimisers2...\n";
    size_t c2_ = r2__rank(M);
    uint8_t* count2_ = new uint8_t[c2_];
    std::memset(count2_, 0, c2_*sizeof(uint8_t));

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view2) {
            if(r2_[minimiser.minimiser_value]) {
                size_t i = r2__rank(minimiser.minimiser_value);
                size_t o = minimiser.occurrences;

                const uint64_t w = o/span + 1;
                if(count2_[i] + w >= m_thres2)
                    count2_[i] = m_thres2;
                else
                    count2_[i] += w;
            }
        }
    }

    std::cout << "filling R2 and HT...\n";

    seqan3::contrib::sdsl::bit_vector r2 = seqan3::contrib::sdsl::bit_vector(M, 0);
    std::unordered_set<uint64_t> freq_kmers;
    std::unordered_set<uint64_t> freq_minimzer;

    auto view3 = srindex::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k, .seed=seed2});

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view3) {
            if(r2_[minimiser.minimiser_value]) {
                size_t i = r2__rank(minimiser.minimiser_value);
                if(count2_[i] < m_thres2)
                    r2[minimiser.minimiser_value] = 1;
                else {
                    freq_kmers.insert(minimiser.window_value);
                    freq_minimzer.insert(minimiser.minimiser_value);
                }
            }
        }
    }

    // for (const auto& kmer: freq_kmers) {
    // // for (const auto& kmer: freq_minimzer) {
    //     seqan3::debug_stream << kmer_to_string(kmer, k);
    //     // seqan3::debug_stream << kmer_to_string(kmer, m) << ' ';
    // }
    // std::cout << '\n';
    std::cout << "no freq minimiser: " << freq_minimzer.size() << '\n';
    std::cout << "no freq kmers: " << freq_kmers.size() << '\n';

    delete[] count2_;

    seqan3::contrib::sdsl::rank_support_v<1> r2_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r2);


    std::cout << "count found minimisers1...\n";
    size_t c1 = r1_rank(M);
    uint8_t* count1 = new uint8_t[c1];
    std::memset(count1, 0, c1*sizeof(uint8_t));
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                count1[i] += minimiser.occurrences/span+1;
            }
        }
    }
    uint64_t n1 = 0;
    for(uint64_t i=0; i < c1; i++)
        n1 += count1[i];
    std::cout << "distinct minimiser1: " << c1 << " total number: " << n1 << '\n';

    std::cout << "count found minimisers2...\n";
    size_t c2 = r2_rank(M);
    uint8_t* count2 = new uint8_t[c2];
    std::memset(count2, 0, c2*sizeof(uint8_t));

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view2) {
            if(r2[minimiser.minimiser_value]) {
                size_t i = r2_rank(minimiser.minimiser_value);
                count2[i] += minimiser.occurrences/span+1;
            }
        }
    }
    uint64_t n2 = 0;
    for(uint64_t i=0; i < c2; i++)
        n2 += count2[i];

    std::cout << "distinct minimiser2: " << c2 << " total number: " << n2 << '\n';
    std::cout << "freq kmers: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";

    uint64_t thresholds1[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    uint64_t* counter1 = new uint64_t[9];
    std::memset(counter1, 0, 9*sizeof(uint64_t));

    for(uint64_t i=0; i < c1; i++) {
        for(int j=8; j >= 0; j--) {
            if(count1[i] >= thresholds1[j]) {
                counter1[j]++;
                break;
            }
        }
    }

    uint64_t thresholds2[18] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100, 1000, 10000, 100000};
    uint64_t* counter2 = new uint64_t[18];
    std::memset(counter2, 0, 18*sizeof(uint64_t));
    for(uint64_t i=0; i < c2; i++) {
        for(int j=17; j >= 0; j--) {
            if(count2[i] >= thresholds2[j]) {
                counter2[j]++;
                break;
            }
        }
    }

    uint64_t cum = 0;
    uint64_t cum_skmers = 0;
    std::cout << "\nminimiser 1 distribution:\n";
    for(int j=0; j < 9; j++) {
        cum += counter1[j];
        cum_skmers += (j+1)*counter1[j];
        std::cout << "occurrences " << thresholds1[j] << ": " << counter1[j] << " " << (double) counter1[j]/c1*100 << "%  cum: " << cum << " " << (double) cum/c1*100 << "% covering " << (double) cum_skmers/n*100 << "% superkmers\n";
    }

    std::cout << "superkmers covered by minimisers 1: " << (double) n1/n*100 << "%\n";
    std::cout << "avg superkmers1: " << (double) n1/c1 <<  '\n';

    std::cout << "minimiser going to level 2: " << c-c1 << "  " << (double) (c-c1)/c*100 << "% to cover " << (double) (n-n1)/n*100 << "% superkmers\n";
    
    cum = 0;
    uint64_t cum_skmers2 = 0;
    std::cout << "\nminimiser 2 distribution:\n";
    for(int j=0; j < 9; j++) {
        cum += counter2[j];
        cum_skmers2 += (j+1)*counter2[j];
        std::cout << "occurrences " << thresholds2[j] << ": " << counter2[j] << " " << (double) counter2[j]/c2*100 << "%  cum: " << cum << " " << (double) cum/c2*100 << "% covering " << (double) cum_skmers2/n*(n-n1)/n*100 << "% (" << (double) cum_skmers2/n*(n-n1)/n*100 + (double) cum_skmers/n*100 << "%) superkmers\n";
    }
    for(int j=9; j < 17; j++) {
        cum += counter2[j];
        cum_skmers2 += (j+1)*counter2[j];
        std::cout << thresholds2[j] << "<= occurrences < " << thresholds2[j+1] << ": " << counter2[j] << " " << (double) counter2[j]/c2*100 << "%  cum: " << cum << " " << (double) cum/c2*100 << "% covering " << (double) cum_skmers2/n*(n-n1)/n*100 << "% (" << (double) cum_skmers2/n*(n-n1)/n*100 + (double) cum_skmers/n*100 << "%) superkmers\n";
    }
    cum += counter2[17];
    std::cout << "occurrences <= " << thresholds2[17] << ": " << counter2[17] << " " << (double) counter2[17]/c2*100 << "%  cum: " << cum << " " << (double) cum/c2*100 << "%\n";

    std::cout << "avg superkmers2: " << (double) n2/c2 <<  '\n';
    std::cout << "superkmers to cover by HT: " << 100.0 - (double) cum_skmers/n*100 - (double) cum_skmers2/n*(n-n1)/n*100 << "%\n";
    std::cout << "kmers in HT: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";

    std::cout << "\n====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no minimiser: " << n << "\n";
    std::cout << "no distinct minimiser: " << c << "\n";
    std::cout << "no minimiser1: " << n1 << "\n";
    std::cout << "no distinct minimiser1: " << c1 << "\n";
    std::cout << "no minimiser2: " << n2 << "\n";
    std::cout << "no distinct minimiser2: " << c2 << "\n";
    std::cout << "freq kmers: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";
    std::cout << "density r1: " << (double) c1/M*100 << "%\n";
    std::cout << "density r2: " << (double) c2/M*100 << "%\n";

    delete[] count1;
    delete[] counter1;
    delete[] count2;
    delete[] counter2;

}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/bigdata/ag_abi/jonas/datasets/human.k31.unitigs.fa.ust.fa.gz";
    // std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}