#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include "../source/minimiser_views.hpp"


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path o{};
    uint8_t k{31};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.q, sharg::config{.short_id = 'q', .long_id = "query", .description = "provide query file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw sharg::user_input_error("provide input file.");
    if(!parser.is_option_set('k'))
        throw sharg::user_input_error("specify k");
    if(!parser.is_option_set('q'))
        throw sharg::user_input_error("provide query file.");

    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
}

int main(int argc, char** argv)
{
    sharg::parser parser{"HT", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
        check_arguments(parser, args);
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::cout << "loading text...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
    load_file(args.i, text);

    std::cout << "building hashtable...\n";
    std::unordered_set<uint64_t> ht;

    for(auto & sequence : text) {
        for(auto && window : sequence | rshash::views::kmerview({.window_size = args.k})) {
            uint64_t kmer = std::min<uint64_t>(window.kmer_value, window.kmer_value_rev);
            ht.insert(kmer);
        }
    }
     
    std::cout << "loading queries...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> queries;
    load_file(args.q, queries);

    std::cout << "querying...\n";
    uint64_t kmers = 0;
    uint64_t found = 0;
    
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
    for (auto query : queries) {
        for (auto && window : query | rshash::views::kmerview({.window_size = args.k})) {
            found += ht.contains(std::min<uint64_t>(window.kmer_value, window.kmer_value_rev));
        }
        kmers += query.size() - args.k + 1;
    }
    std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        
    double ns_per_kmer = (double) elapsed.count() / kmers;
        
    std::cout << "==== query report:\n";
    std::cout << "num_kmers = " << kmers << '\n';
    std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers*100 << "%)\n";
    std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
 
    return 0;
}