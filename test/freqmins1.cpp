#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <seqan3/core/debug_stream.hpp>


#include "../source/minimiser_rev_xor_hash_views.hpp"


const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 1;

const uint8_t m_thres1 = 20;
const uint8_t m_thres2 = 20;


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

static inline constexpr uint64_t compute_mask(uint64_t const size) {
        assert(size > 0u);
        assert(size <= 64u);

        if (size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }


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


void stats(const std::vector<std::vector<seqan3::dna4>> &input)
{
    const uint8_t k = 31;
    uint8_t m = 16;
    std::cout << "m: " << +m << '\n';

    const size_t span = 100;

    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m, .window_size = k, .seed=seed1});
    auto view2 = srindex::views::two_minimisers_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});
    auto view3 = srindex::views::two_minimisers_and_occurence_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});
    auto view4 = srindex::views::two_minimisers_and_window_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});

    const uint64_t M = 1ULL << (m+m);

    std::cout << "scanning minimizers...\n";
    seqan3::contrib::sdsl::bit_vector r = seqan3::contrib::sdsl::bit_vector(M, 0);
    
    for(auto & record : input) {
        for(auto && minimiser : record | view1)
            r[minimiser.minimiser_value] = 1;
    }

    std::cout << "counting minimisers...\n";
    seqan3::contrib::sdsl::rank_support_v<1> r_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r);
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    size_t N = 0;
    uint64_t n = 0;
    uint64_t kmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;
            const uint8_t w = o/span + 1;
            count[i] += w;
            if(count[i] > m_thres1)
                count[i] = m_thres1;
            kmers += o;
            n += w;
        }
        N += sequence.size();
    }

    std::cout << "mark (un)frequent minimisers...\n";
    seqan3::contrib::sdsl::bit_vector r1 = seqan3::contrib::sdsl::bit_vector(M, 0);
    seqan3::contrib::sdsl::bit_vector r2 = seqan3::contrib::sdsl::bit_vector(M, 0);
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view2) {
            size_t i = r_rank(minimiser.minimiser1_value);
            if(count[i] < m_thres1)
                r1[minimiser.minimiser1_value] = 1;
            else
                r2[minimiser.minimiser2_value] = 1;
        }
    }

    std::cout << "counting minimisers1...\n";
    seqan3::contrib::sdsl::rank_support_v<1> r1_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r1);
    size_t c1 = r1_rank(M);
    uint8_t* count1 = new uint8_t[c1];
    std::memset(count1, 0, c1*sizeof(uint8_t));

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value])
                count1[r1_rank(minimiser.minimiser_value)] = count[r_rank(minimiser.minimiser_value)];
        }
    }
    uint64_t n1 = 0;
    for(uint64_t i=0; i < c1; i++)
        n1 += count1[i];

    delete[] count;
    // delete r

    std::cout << "counting minimisers2...\n";
    seqan3::contrib::sdsl::rank_support_v<1> r2_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r2);
    size_t c2 = r2_rank(M);
    uint8_t* count2 = new uint8_t[c2];
    std::memset(count2, 0, c2*sizeof(uint8_t));

    auto update_count2 = [&](uint64_t minimiser_value, size_t occurrences) {
        size_t i = r2_rank(minimiser_value);
        count2[i] += occurrences/span + 1;
        if(count2[i] > m_thres2)
            count2[i] = m_thres2;
    };

    for(auto & sequence : input) {
        size_t delta;
        for(auto && minimiser : sequence | view3)
        {
            bool level1 = r1[minimiser.minimiser1_value];
            bool next_level1 = r1[minimiser.new_minimiser1_value];
            bool m2_change = minimiser.minimiser2_value != minimiser.new_minimiser2_value;

            // in level 2 update count2 when going to level 1 or minimiser 2 changes
            if(!level1) {
                if(next_level1 || m2_change)
                    // update_count2(minimiser.minimiser2_value, minimiser.occurrences2);
                    update_count2(minimiser.minimiser2_value, minimiser.occurrences2 - delta);
            }
            // // when changing level update occurences i.e. delta of minimiser in level 2 respectively
            if(level1 != next_level1)
                delta = minimiser.occurrences2;
            // if minimiser 2 changes set delta to zero
            delta *= !m2_change;
        }
    }


    std::cout << "mark unfrequent minimisers2, filling HT...\n";
    seqan3::contrib::sdsl::bit_vector r3 = seqan3::contrib::sdsl::bit_vector(M, 0);
    std::unordered_set<uint64_t> freq_kmers;
    std::unordered_set<uint64_t> freq_minimzer;

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view3) {
            if(!r1[minimiser.minimiser1_value]) {
                size_t i = r2_rank(minimiser.minimiser2_value);
                if(count2[i] < m_thres2)
                    r3[minimiser.minimiser2_value] = 1;
                else {
                    // freq_kmers.insert(minimiser.window_value);
                    freq_minimzer.insert(minimiser.minimiser2_value);
                }
            }
        }
    }

    std::cout << "counting unfrequent minimisers2...\n";
    seqan3::contrib::sdsl::rank_support_v<1> r3_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r3);
    size_t c3 = r3_rank(M);
    uint8_t* count3 = new uint8_t[c3];
    std::memset(count3, 0, c3*sizeof(uint8_t));

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view3) {
            if(r3[minimiser.minimiser2_value])
                count3[r3_rank(minimiser.minimiser2_value)] = count2[r2_rank(minimiser.minimiser2_value)];
        }
    }
    uint64_t n2 = 0;
    for(uint64_t i=0; i < c3; i++)
        n2 += count3[i];

    // delete r2
    delete[] count2;

    std::cout << "computing distributions...\n";
    uint64_t* counter1 = new uint64_t[m_thres1];
    std::memset(counter1, 0, m_thres1*sizeof(uint64_t));
    for(uint64_t i=0; i < c1; i++) {
        for(int j=0; j < m_thres1; j++) {
            if(count1[i] == j+1) {
                counter1[j]++;
                break;
            }
        }
    }
    uint64_t* counter2 = new uint64_t[m_thres2];
    std::memset(counter2, 0, m_thres2*sizeof(uint64_t));
    for(uint64_t i=0; i < c3; i++) {
        for(int j=0; j < m_thres2; j++) {
            if(count3[i] == j+1) {
                counter2[j]++;
                break;
            }
        }
    }

    uint64_t cum = 0;
    uint64_t cum_skmers = 0;
    std::cout << "\nminimiser 1 distribution:\n";
    for(int j=0; j < m_thres1; j++) {
        cum += counter1[j];
        cum_skmers += (j+1)*counter1[j];
        std::cout << "occurrences " << (j+1) << ": " << counter1[j] << " " << (double) counter1[j]/c1*100 << "%  cum: " << cum << " " << (double) cum/c1*100 << "% covering " << cum_skmers << ' ' << (double) cum_skmers/n*100 << "% superkmers\n";
    }
    std::cout << "avg superkmers1: " << (double) n1/c1 <<  '\n';
    std::cout << "minimisers going to level 2: " << c-c1 << "  " << (double) (c-c1)/c*100 << "% to cover " << (double) (n-n1)/n*100 << "% superkmers\n";
    
    cum = 0;
    std::cout << "\nminimiser 2 distribution:\n";
    for(int j=0; j < m_thres2; j++) {
        cum += counter2[j];
        cum_skmers += (j+1)*counter2[j];
        std::cout << "occurrences " << (j+1) << ": " << counter2[j] << " " << (double) counter2[j]/c3*100 << "%  cum: " << cum << " " << (double) cum/c3*100 << "% covering " << cum_skmers << ' ' << (double) cum_skmers/n*100 << "% superkmers\n";
    }
    std::cout << "avg superkmers2: " << (double) n2/c3 <<  '\n';
    // std::cout << "minimisers going to HT: " << c1-c2 << " " << (double) freq_minimzer.size()/(n2)*100 << "%\n";

    // std::cout << "superkmers to cover by HT: " << 100.0 - (double) cum_skmers/(n1+n2+n3)*100 << "%\n";
    // std::cout << "freq minimisers left: " << freq_minimzer.size() << " " << (double) freq_minimzer.size()/n*100 << "%\n";
    // std::cout << "kmers in HT: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";

    std::cout << "\n====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no minimiser: " << n << "\n";
    // std::cout << "no distinct minimiser: " << c << "\n";
    std::cout << "no minimiser1: " << n1 << "\n";
    std::cout << "no distinct minimiser1: " << c1 << "\n";
    std::cout << "no minimiser2: " << n2 << "\n";
    std::cout << "no distinct minimiser2: " << c3 << "\n";
    std::cout << "no minimiser HT: " << freq_minimzer.size() << "\n";
    std::cout << "no kmers HT: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";
    std::cout << "density r1: " << (double) c1/M*100 << "%\n";
    std::cout << "density r2: " << (double) c3/M*100 << "%\n";

    delete[] count1;
    delete[] counter1;
    delete[] count3;
    delete[] counter2;

}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}