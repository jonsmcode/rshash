#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

#include "../source/minimiser_rev_xor_views3.hpp"

const uint8_t k = 31;
const uint8_t m1 = 15;
const uint8_t m2 = 16;
const uint8_t m3 = 18;
const uint8_t m_thres = 20;
const uint8_t span = 31;


const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed3 = 0xE5'9A'38'5F'03'76'C9'F6;


namespace sdsl = seqan3::contrib::sdsl;


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
        if(N >= 1000000000)
            break;
        // if(N >= 100000000)
        //     break;
    }
}


void stats(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});

    std::cout << +m1 << " " << +m2 << " " << +m3 << "\n";

    const uint64_t M1 = 1ULL << (m1+m1);
    const uint64_t M2 = 1ULL << (m2+m2);
    const uint64_t M3 = 1ULL << (m3+m3);

    std::cout << "find minimizers...\n";
    sdsl::bit_vector r1tmp = sdsl::bit_vector(M1, 0);

    size_t N = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        for(auto && minimiser : record | view1) {
            r1tmp[minimiser.minimiser_value] = 1;
        }
        N += record.size();
        no_sequences++;
    }
    sdsl::rank_support_v<1> r1tmp_rank = sdsl::rank_support_v<1>(&r1tmp);

    std::cout << "count minimizers...\n";
    size_t c1tmp = r1tmp_rank(M1);
    uint8_t* count1tmp = new uint8_t[c1tmp];
    std::memset(count1tmp, 0, c1tmp*sizeof(uint8_t));

    uint64_t kmers = 0;
    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            size_t i = r1tmp_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;
            size_t w = o/span + 1;
            count1tmp[i] += w;
            if(count1tmp[i] > m_thres)
                count1tmp[i] = m_thres;
            kmers += o;
            n += w;
        }
    }

    std::cout << "filling R_1...\n";
    sdsl::bit_vector r1 = sdsl::bit_vector(M1, 0);
    for(auto & sequence : input) {
        for(auto && minimisers : sequence | view1) {
            if(count1tmp[r1tmp_rank(minimisers.minimiser_value)] < m_thres)
                r1[minimisers.minimiser_value] = 1;
        }
    }
    sdsl::rank_support_v<1> r1_rank = sdsl::rank_support_v<1>(&r1);

    delete[] count1tmp;

    std::cout << "count minimizers1 again...\n";
    size_t c1 = r1_rank(M1);
    uint8_t* count1 = new uint8_t[c1];
    std::memset(count1, 0, c1*sizeof(uint8_t));
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                count1[i] += minimiser.occurrences/span + 1;
            }
        }
    }
    uint64_t n1 = 0;
    for(size_t i=0; i < c1; i++)
        n1 += count1[i];

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers1;
    std::vector<size_t> skmer_positions;
    size_t length = 0;
    for(auto & sequence : input) {
        size_t start_position;
        bool level_up;
        bool current_level_up;

        for(auto && minimiser : sequence | view1) {
            level_up = r1[minimiser.minimiser_value];
            break;
        }
        if(!level_up)
            start_position = 0;
        for(auto && minimiser : sequence | view1) {
            current_level_up = r1[minimiser.minimiser_value];

            if(level_up && !current_level_up)
                start_position = minimiser.range_position;
            if(!level_up && current_level_up) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers1.push_back(skmer);
                skmer_positions.push_back(length + start_position);
            }

            level_up = current_level_up;
        }
        if(!current_level_up) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers1.push_back(skmer);
            skmer_positions.push_back(length + start_position);
        }

        length += sequence.size();
    }
    
    size_t len_rem_seqs1 = 0;
    for(auto & sequence : freq_skmers1)
        len_rem_seqs1 += sequence.size();


    std::cout << "level 2...\n";
    auto view2 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
    sdsl::bit_vector r2tmp = sdsl::bit_vector(M2, 0);

    for(auto & record : freq_skmers1) {
        for(auto && minimiser : record | view2)
            r2tmp[minimiser.minimiser_value] = 1;
    }
    sdsl::rank_support_v<1> r2tmp_rank = sdsl::rank_support_v<1>(&r2tmp);

    std::cout << "count minimizers2...\n";
    size_t c2tmp = r2tmp_rank(M2);
    uint8_t* count2tmp = new uint8_t[c2tmp];
    std::memset(count2tmp, 0, c2tmp*sizeof(uint8_t));

    uint64_t n2tmp = 0;
    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view2) {
            size_t i = r2tmp_rank(minimiser.minimiser_value);
            count2tmp[i] += minimiser.occurrences/span + 1;
            n2tmp += minimiser.occurrences/span + 1;
            if(count2tmp[i] > m_thres)
                count2tmp[i] = m_thres;
        }
    }

    std::cout << "fill R2...\n";
    sdsl::bit_vector r2 = sdsl::bit_vector(M2, 0);
    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view2) {
            if(count2tmp[r2tmp_rank(minimiser.minimiser_value)] < m_thres)
                r2[minimiser.minimiser_value] = 1;
        }
    }
    sdsl::rank_support_v<1> r2_rank = sdsl::rank_support_v<1>(&r2);

    std::cout << "fill count 2...\n";
    size_t c2 = r2_rank(M2);
    uint8_t* count2 = new uint8_t[c2];
    std::memset(count2, 0, c2*sizeof(uint8_t));

    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view2) {
            count2[r2_rank(minimiser.minimiser_value)] = count2tmp[r2tmp_rank(minimiser.minimiser_value)];
        }
    }

    delete[] count2tmp;

    uint64_t n2 = 0;
    for(uint64_t i=0; i < c2; i++) {
        n2 += count2[i];
    }


    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers2;
    for(auto & sequence : freq_skmers1) {
        size_t start_position;
        bool level_up;
        bool current_level_up;

        for(auto && minimiser : sequence | view2) {
            level_up = r2[minimiser.minimiser_value];
            break;
        }
        if(!level_up)
            start_position = 0;
        for(auto && minimiser : sequence | view2) {
            current_level_up = r2[minimiser.minimiser_value];

            if(level_up && !current_level_up)
                start_position = minimiser.range_position;
            if(!level_up && current_level_up) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers2.push_back(skmer);
            }

            level_up = current_level_up;
        }
        if(!current_level_up) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers2.push_back(skmer);
        }

    }
    
    size_t len_rem_seqs2 = 0;
    for(auto & sequence : freq_skmers2)
        len_rem_seqs2 += sequence.size();
    std::cout << "remaining superkmers " << freq_skmers2.size() << " (" << (double) freq_skmers2.size()/n*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs2 << " (" << (double) len_rem_seqs2/N*100 << "%)\n";


    std::cout << "level 3...\n";

    auto view3 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m3, .window_size = k, .seed=seed3});
    sdsl::bit_vector r3tmp = sdsl::bit_vector(M3, 0);

    for(auto & record : freq_skmers2) {
        for(auto && minimiser : record | view3)
            r3tmp[minimiser.minimiser_value] = 1;
    }
    sdsl::rank_support_v<1> r3tmp_rank = sdsl::rank_support_v<1>(&r3tmp);

    std::cout << "count minimizers3...\n";
    size_t c3tmp = r3tmp_rank(M3);
    uint8_t* count3tmp = new uint8_t[c3tmp];
    std::memset(count3tmp, 0, c3tmp*sizeof(uint8_t));

    uint64_t n3tmp = 0;
    for(auto & sequence : freq_skmers2) {
        for(auto && minimiser : sequence | view3) {
            size_t i = r3tmp_rank(minimiser.minimiser_value);
            count3tmp[i] += minimiser.occurrences/span + 1;
            n3tmp += minimiser.occurrences/span + 1;
            if(count3tmp[i] > m_thres)
                count3tmp[i] = m_thres;
        }
    }

    std::cout << "fill R3...\n";
    sdsl::bit_vector r3 = sdsl::bit_vector(M3, 0);
    for(auto & sequence : freq_skmers2) {
        for(auto && minimiser : sequence | view3) {
            if(count3tmp[r3tmp_rank(minimiser.minimiser_value)] < m_thres)
                r3[minimiser.minimiser_value] = 1;
        }
    }
    sdsl::rank_support_v<1> r3_rank = sdsl::rank_support_v<1>(&r3);

    std::cout << "fill count 3...\n";
    size_t c3 = r3_rank(M3);
    uint8_t* count3 = new uint8_t[c3];
    std::memset(count3, 0, c3*sizeof(uint8_t));

    for(auto & sequence : freq_skmers2) {
        for(auto && minimiser : sequence | view3) {
            count3[r3_rank(minimiser.minimiser_value)] = count3tmp[r3tmp_rank(minimiser.minimiser_value)];
        }
    }

    delete[] count3tmp;

    uint64_t n3 = 0;
    for(size_t i=0; i < c3; i++) {
        n3 += count3[i];
    }

    std::cout << "build level 4, HT...\n";

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers3;
    for(auto & sequence : freq_skmers2) {
        size_t start_position;
        bool level_up;
        bool current_level_up;

        for(auto && minimiser : sequence | view3) {
            level_up = r3[minimiser.minimiser_value];
            break;
        }
        if(!level_up)
            start_position = 0;
        for(auto && minimiser : sequence | view3) {
            current_level_up = r3[minimiser.minimiser_value];

            if(level_up && !current_level_up)
                start_position = minimiser.range_position;
            if(!level_up && current_level_up) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers3.push_back(skmer);
            }

            level_up = current_level_up;
        }
        if(!current_level_up) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers3.push_back(skmer);
        }
    }

    size_t len_rem_seqs3 = 0;
    for(auto & sequence : freq_skmers3)
        len_rem_seqs3 += sequence.size();
    std::cout << "remaining superkmers " << freq_skmers3.size() << " (" << (double) freq_skmers3.size()/n*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs3 << " (" << (double) len_rem_seqs3/N*100 << "%)\n";

    std::cout << "filling HT...\n";
    std::unordered_set<uint64_t> freq_kmers;
    auto view4 = srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1});
    for(auto & sequence : freq_skmers3) {
        for(auto && minimiser : sequence | view4) {
            freq_kmers.insert(std::min<uint64_t>(minimiser.window_value, minimiser.window_value_rev));
        }
    }

    uint64_t counter1[m_thres] = {0};
    for(uint64_t i=0; i < c1; i++) {
        for(int j=m_thres-1; j >= 0; j--) {
            if(count1[i] >= j+1) {
                counter1[j]++;
                break;
            }
        }
    }

    uint64_t counter2[m_thres] = {0};
    for(uint64_t i=0; i < c2; i++) {
        for(int j=m_thres-1; j >= 0; j--) {
            if(count2[i] >= j+1) {
                counter2[j]++;
                break;
            }
        }
    }

    uint64_t counter3[m_thres] = {0};
    for(uint64_t i=0; i < c3; i++) {
        for(int j=m_thres-1; j >= 0; j--) {
            if(count3[i] >= j+1) {
                counter3[j]++;
                break;
            }
        }
    }

    std::cout << "\n====== report ======\n";
    std::cout << "text:\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "#kmers: " << kmers <<  '\n';
    std::cout << "#ukmers: todo\n";
    std::cout << "#kmers/#ukmers: todo\n\n";

    std::cout << "Level1: m1=" << +m1 << ", threshold=" << +m_thres << "\n";
    std::cout << "#minis1: " << n << "\n";
    std::cout << "#uminis1: " << c1tmp << "\n";
    std::cout << "#skmers1: todo\n";
    std::cout << "#unfreq minis1: " << n1 << " (" << (double) n1/n*100 << "%)\n";
    std::cout << "#unfreq uminis1: " << c1 << " ("<< (double) c1/c1tmp*100 << "%)\n\n";
    uint64_t cum = 0;
    uint64_t cum_skmers = 0;
    for(int j=0; j < m_thres; j++) {
        cum += counter1[j];
        cum_skmers += (j+1)*counter1[j];
        std::cout << "occurrences " << (j+1) << ": " << counter1[j] << " " << (double) counter1[j]/c1tmp*100 << "%  cum: " << cum << " " << (double) cum/c1tmp*100 << "% covering " << cum_skmers << ' ' << (double) cum_skmers/n*100 << "% windows\n";
    }
    std::cout << "avg windows/unfreq minis1: " << (double) n1/c1 <<  '\n';

    std::cout << "uminis1 going to level 2: " << c1tmp-c1 << '\n';
    std::cout << "to cover " << freq_skmers1.size() << " remaining skmers1";
    std::cout << " (average length: " << (double) len_rem_seqs1/freq_skmers1.size()-span << ") in " << n-n1 << " many windows\n";
    std::cout << "split into: " << n2tmp << " minis2\n\n";

    std::cout << "Level2: m2=" << +m2 << ", threshold=" << +m_thres << "\n";
    std::cout << "#minis2: " << n2tmp << "\n";
    std::cout << "#uminis2: " << c2tmp << "\n";
    std::cout << "#skmers2: todo\n";
    std::cout << "#unfreq minis2: " << n2 << " (" << (double) n2/n2tmp*100 << "%)\n";
    std::cout << "#unfreq uminis2: " << c2 << " ("<< (double) c2/c2tmp*100 << "%)\n";
    cum = 0;
    uint64_t cum_skmers2 = 0;
    for(int j=0; j < m_thres; j++) {
        cum += counter2[j];
        cum_skmers2 += (j+1)*counter2[j];
        std::cout << "occurrences " << (j+1) << ": " << counter2[j] << " " << (double) counter2[j]/c2tmp*100 << "%  cum: " << cum << " " << (double) cum/c2tmp*100 << "% covering " << cum_skmers2 << ' ' << (double) cum_skmers2/n2tmp*100 << "% windows\n";
    }
    std::cout << "avg windows/unfreq minis2: " << (double) n2/c2 <<  '\n';

    std::cout << "uminis2 going to level 3: " << c2tmp-c2 << '\n';
    std::cout << "to cover " << freq_skmers2.size() << " remaining skmers2";
    std::cout << " (average length: " << (double) len_rem_seqs2/freq_skmers2.size()-span << ") in " << n2tmp-n2 << " many windows\n";
    std::cout << "split into " << n3tmp << " minis2\n\n";

    std::cout << "Level3: m3=" << +m3 << ", threshold=" << +m_thres << "\n";
    std::cout << "#minis3: " << n3tmp << "\n";
    std::cout << "#uminis3: " << c3tmp << "\n";
    std::cout << "#skmers3: todo\n";
    std::cout << "#unfreq minis3: " << n3 << " (" << (double) n3/n3tmp*100 << "%)\n";
    std::cout << "#unfreq uminis3: " << c3 << " ("<< (double) c3/c3tmp*100 << "%)\n\n";
    cum = 0;
    uint64_t cum_skmers3 = 0;
    std::cout << "\nminimiser 3 distribution:\n";
    for(int j=0; j < m_thres; j++) {
        cum += counter3[j];
        cum_skmers3 += (j+1)*counter3[j];
        std::cout << "occurrences " << (j+1) << ": " << counter3[j] << " " << (double) counter3[j]/c3tmp*100 << "%  cum: " << cum << " " << (double) cum/c3tmp*100 << "% covering " << cum_skmers3 << ' ' << (double) cum_skmers3/n3tmp*100 << "% windows\n";
    }
    std::cout << "avg skmers: " << (double) n3/c3 <<  '\n';
    std::cout << "remaining skmers " << freq_skmers3.size() << " (avgerage length: " << (double) len_rem_seqs3/freq_skmers3.size()-span << ") in " << n3tmp-n3 << " many windows\n\n";

    std::cout << "kmers in HT: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";
    

    delete[] count1;
    delete[] count2;
    delete[] count3;
}



int main(int argc, char** argv)
{
    // std::filesystem::path path = "/bigdata/ag_abi/jonas/datasets/human.k31.unitigs.fa.ust.fa.gz";
    std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}