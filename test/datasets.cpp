#include <filesystem>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

using namespace seqan3::literals;

uint64_t const seed = 0;


void random_dna(std::vector<seqan3::dna4> &text, const uint64_t length) {
    for(int i=0; i < length; i++) {
        seqan3::dna4 d;
        char r = rand() % 4;
        text.push_back(d.assign_rank(r));
    }
}

void different_dna(std::vector<seqan3::dna4> &text, std::vector<seqan3::dna4> &difftext) {
    char r;
    seqan3::dna4 d;
    for(int i=0; i < text.size(); i++) {
        r = (text[i].to_rank() + 1) % 4;
        difftext.push_back(d.assign_rank(r));
    }
}

std::vector<seqan3::dna4> sample_dna(std::vector<seqan3::dna4> &text, const uint64_t position, const uint64_t length) {
    std::vector<seqan3::dna4> sample;
    for(uint64_t i=0; i < length; i++) {
        sample.push_back(text[position+i]);
    }
    return sample;
}

// void guard_kmer(std::vector<seqan3::dna4> &text, std::vector<seqan3::dna4> &guarded_kmer, std::vector<seqan3::dna4> &kmer, const uint8_t k) {
//     char r;
//     seqan3::dna4 d;
//     for(int i=0; i < k; i++) {
//         r = (text[i].to_rank() + 1) % 4;
//         guarded_kmer.push_back(d.assign_rank(r));
//     }
//     for(int i=0; i < k; i++)
//         guarded_kmer.push_back(kmer[i]);
//     r = (kmer[k-1].to_rank() + 1) % 4;
//     for(int i=0; i < k; i++)
//         guarded_kmer.push_back(d.assign_rank(r));
// }


// void random_positions(std::vector<uint64_t> &positions, const uint64_t n, const uint8_t k, const uint64_t range) {
//     if(n == 0)
//         return;
//     if(n > range/k) {
//         std::cout << "n too big";
//         return;
//     }
//     uint64_t p = rand() % range;
//     positions.push_back(p);
//     int i = 1;
//     while(i < n) {
//         p = rand() % range;
//         auto it = std::upper_bound(positions.begin(), positions.end(), p) - positions.begin();
//         if(it == 0) {
//             if(labs(p - positions[it]) > k) {
//                 positions.insert(positions.begin()+it, p);
//                 i++;
//             }
//         }
//         else {
//             if(labs(p - positions[it-1]) > k && labs(p - positions[it]) > k) {
//                 positions.insert(positions.begin()+it, p);
//                 i++;
//             }
//         }
//     }
// }

// void random_positions(std::vector<uint64_t> &positions, const uint64_t n, const uint8_t k, const uint64_t range) {
//     if(n == 0)
//         return;
//     if(n > range/k) {
//         std::cout << "n too big";
//         return;
//     }
//     uint64_t p = rand() % range;
//     positions.push_back(p);
//     int i = 1;
//     while(i < n) {
//         p = rand() % range;
//         auto it = std::upper_bound(positions.begin(), positions.end(), p) - positions.begin();
//         if(it == 0) {
//             if(labs(p - positions[it]) > k) {
//                 positions.insert(positions.begin()+it, p);
//                 i++;
//             }
//         }
//         else {
//             if(labs(p - positions[it-1]) > k && labs(p - positions[it]) > k) {
//                 positions.insert(positions.begin()+it, p);
//                 i++;
//             }
//         }
//     }
// }


void insert_kmer(std::vector<seqan3::dna4> &text, const std::vector<seqan3::dna4> &kmer, const uint64_t position)
{
    for(uint64_t i=0; i < kmer.size(); i++)
        text[position+i] = kmer[i];
}

void insert_kmer(std::vector<seqan3::dna4> &text, const std::vector<seqan3::dna4> &kmer,
    const std::vector<uint64_t> &positions)
{
    for(uint64_t p : positions)
        for(size_t i=0; i < kmer.size(); i++)
            text[p+i] = kmer[i];
}

void bf(const std::vector<seqan3::dna4> &text,
        const std::vector<std::vector<seqan3::dna4>> &queries, const uint8_t k,
        std::vector<std::vector<uint64_t>> &positions)
{
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k})
                      | std::views::transform([](uint64_t i) {return i ^ seed;});
    for(auto query : queries) {
        std::vector<uint64_t> pos;
        for (auto && qkmer : query | kmer_view) {
            uint64_t i = 0;
            for (auto && tkmer : text | kmer_view) {
                if(qkmer == tkmer)
                    pos.push_back(i);
                i++;
            }
        }
        positions.push_back(pos);
    }
}

void save_dataset(const std::string file, const std::vector<seqan3::dna4> &text,
                  const std::vector<std::vector<seqan3::dna4>> &queries,
                  const std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &positions)
{
    using types = seqan3::type_list<std::vector<seqan3::dna4>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;

    std::ofstream textout(file + "_text.fasta");
    textout << ">text\n";
    for (auto nucleotide : text)
        textout << seqan3::to_char(nucleotide);
    textout << '\n';
    // seqan3::sequence_file_output tfout{textout, seqan3::format_fasta{}};
    // sequence_record_type recordt{std::move(text), std::move("text")};
    // tfout.push_back(recordt);

    std::ofstream queryout(file + "_query.fasta");
    seqan3::sequence_file_output qfout{queryout, seqan3::format_fasta{}};
    int i = 0;
    for(auto query : queries) {
        queryout << ">query-" + std::to_string(i) + '\n';
        for (auto nucleotide : query)
            queryout << seqan3::to_char(nucleotide);
        queryout << '\n';
        // sequence_record_type recordq{std::move(query), std::move("query-"+std::to_string(i))};
        // qfout.push_back(recordq);
        i++;
    }

    std::ofstream positionsfile;
    positionsfile.open(file + ".positions");
    for (auto const& [query, querypos, sequence, textposition] : positions) {
        positionsfile << '(' << query << ',' << querypos << ',' << sequence << ',' << textposition << "), ";
    }
    positionsfile.close();
}


void build_dataset(const uint64_t textlength, const uint64_t querylength, const uint64_t no_queries, const uint8_t k,
                   const double hitrate, const std::string path)
{
    std::vector<seqan3::dna4> text;
    std::vector<std::vector<seqan3::dna4>> queries;

    random_dna(text, textlength);

    // sample snippets of length sample_size from text that are sample_dist apart
    const uint64_t sample_size = 500;
    const uint64_t sample_dist = textlength/5000;

    std::vector<std::vector<seqan3::dna4>> samples;
    for(uint64_t i = 0; i < textlength-sample_size; i += sample_dist+sample_size) {
        samples.push_back(sample_dna(text, i, sample_size));
    }

    // insert text samples into queries with query_dist characters apart
    const uint64_t query_dist = (sample_size-k+1)/hitrate - sample_size-k+1;

    uint64_t occs = 0;
    for(uint64_t q = 0; q < no_queries; q++) {
        std::vector<seqan3::dna4> query;
        random_dna(query, querylength);

        uint64_t j = 0;
        for(uint64_t i = 0; j < querylength-sample_size; i++) {
            std::vector<seqan3::dna4> sample = samples[i%samples.size()];
            insert_kmer(query, sample, j);
            j += query_dist;
            occs += sample_size-k+1;
        }
        queries.push_back(query);
    }

    std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> positions;

    std::string dataset = "random_n" + std::to_string(textlength) + "_o" + std::to_string(occs);
    std::string file = path + dataset;

    save_dataset(file, text, queries, positions);
}



int main(int argc, char** argv)
{
    const uint64_t textlengths[] = {10000000, 100000000, 1000000000};
    const uint8_t k = 31;
    const double hitrates[] = {0.001, 0.01, 0.1};

    for(const uint64_t textlength : textlengths) {
        for(const double hitrate : hitrates) {
            std::cout << "building random dataset with text length " << textlength << " with hitrate " << hitrate << "\n";
            const uint64_t numberqueries = 100;
            const uint64_t querylength = textlength/numberqueries;
            build_dataset(textlength, querylength, numberqueries, k, hitrate, "../test/datasets/synthetic/");
        }
    }

    return 0;
}

