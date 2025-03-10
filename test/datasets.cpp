#include <filesystem>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

using namespace seqan3::literals;

uint64_t const seed = 0;


void random_dna(std::vector<seqan3::dna4> &text, const int length) {
    for(int i=0; i < length; i++) {
        seqan3::dna4 d;
        char r = rand() % 4;
        text.push_back(d.assign_rank(r));
    }
}

void random_positions(std::vector<int> &positions, const int n, const int k, const int range) {
    if(n == 0)
        return;
    int p = rand() % range;
    positions.push_back(p);
    int i = 1;
    while(i < n) {
        p = rand() % range;
        auto it = std::upper_bound(positions.begin(), positions.end(), p) - positions.begin();
        if(it == 0) {
            if(abs(p - positions[it]) > k) {
                positions.insert(positions.begin()+it, p);
                i++;
            }
        }
        else {
            if(abs(p - positions[it-1]) > k && abs(p - positions[it]) > k) {
                positions.insert(positions.begin()+it, p);
                i++;
            }
        }
    }
}

void insert_kmer(std::vector<seqan3::dna4> &text, const std::vector<seqan3::dna4> &kmer,
    const std::vector<int> &positions)
{
    for(int p : positions)
        for(int i=0; i < kmer.size(); i++)
            text[p+i] = kmer[i];
}

void bf(const std::vector<seqan3::dna4> &text,
        const std::vector<std::vector<seqan3::dna4>> &queries, const int k,
        std::vector<std::vector<int>> &positions)
{
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k})
                      | std::views::transform([](uint64_t i) {return i ^ seed;});
    for(auto query : queries) {
        std::vector<int> pos;
        for (auto && qkmer : query | kmer_view) {
            int i = 0;
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
                  const std::vector<std::vector<int>> &positions)
{
    using types = seqan3::type_list<std::vector<seqan3::dna4>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;

    std::ofstream textout(file + "_text.fasta");
    seqan3::sequence_file_output tfout{textout, seqan3::format_fasta{}};
    sequence_record_type recordt{std::move(text), std::move("text")};
    tfout.push_back(recordt);

    std::ofstream queryout(file + "_query.fasta");
    seqan3::sequence_file_output qfout{queryout, seqan3::format_fasta{}};
    int i = 0;
    for(auto query : queries) {
        sequence_record_type recordq{std::move(query), std::move("query-"+std::to_string(i))};
        qfout.push_back(recordq);
        i++;
    }

    std::ofstream positionsfile;
    positionsfile.open(file + ".positions");
    for(auto position : positions) {
        for(int p : position)
            positionsfile << p << ' ';
        positionsfile << '\n';
    }
    positionsfile.close();
}


// void build_random_dataset1(const std::string path, const int id,
//     const int k, const int seq_length, const int query_length, const int no_queries, const int n)
// {
//     std::vector<seqan3::dna4> text;
//     std::vector<seqan3::dna4> query;
//     std::vector<std::vector<seqan3::dna4>> queries;
//     std::vector<seqan3::dna4> kmer;
//     std::vector<int> tpositions;
//     std::vector<int> qpositions;

//     random_dna(kmer, k);
//     random_dna(text, seq_length);
//     random_positions(tpositions, n, k, seq_length-k-1);
//     insert_kmer(text, kmer, tpositions);
//     random_dna(query, query_length);
//     queries.push_back(query);
//     random_positions(qpositions, 1, k, query_length-k-1);
//     insert_kmer(query, kmer, qpositions);

//     std::string dataset = "r-" + std::to_string(id) + "_k" + std::to_string(k);
//     std::string file = path + dataset;

//     save_dataset(file, text, queries, tpositions);
// }


void build_random_dataset(const std::string path, const int id,
    const int k, const int seq_length, const int query_length, const int no_queries, const int n)
{
    std::vector<seqan3::dna4> text;
    std::vector<std::vector<seqan3::dna4>> queries;
    std::vector<seqan3::dna4> kmer;
    std::vector<int> tpositions;
    std::vector<std::vector<int>> positions;

    random_dna(kmer, k);
    random_dna(text, seq_length);
    random_positions(tpositions, n, k, seq_length-k-1);
    insert_kmer(text, kmer, tpositions);

    for(int q = 0; q < no_queries; q++) {
        std::vector<seqan3::dna4> query;
        std::vector<int> qpositions;
        random_dna(query, query_length);
        random_positions(qpositions, 1, k, query_length-k-1);
        insert_kmer(query, kmer, qpositions);
        queries.push_back(query);
    }

    bf(text, queries, k, positions);

    std::string dataset = "r-" + std::to_string(id) + "_k" + std::to_string(k);
    std::string file = path + dataset;

    save_dataset(file, text, queries, positions);
}




int main(int argc, char** argv)
{
    int ks[] = {20, 31};
    int seq_lengths[] = {100000, 1000000, 10000000};

    int n = 100;

    int id = 0;
    for(int k : ks) {
        for(int tl : seq_lengths) {
                std::cout << "building random lookup dataset " << id << "\n";
                build_random_dataset("test/datasets/1-lookup-random/", id, k, tl, k, 1, n);
                id++;
            }
    }

    n = 10;
    int ql = 150;

    id = 0;
    for(int k : {31}) {
        for(int qn : {1000, 10000}) {
            for(int tl : seq_lengths) {
                std::cout << "building random stream dataset " << id << "\n";
                build_random_dataset("test/datasets/2-stream-random/", id, k, tl, ql, qn, n);
                id++;
            }
        }
    }

    return 0;
}

