#include <seqan3/core/debug_stream.hpp>
#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "rsindex3_simd.hpp"
#include "bench.hpp"

struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path d{};
    uint8_t k{};
    uint8_t m1{17};
    uint8_t m2{17};
    uint8_t m3{17};
    uint8_t t1{64};
    uint8_t t2{64};
    uint16_t t3{64};
    size_t s{15};
    bool c{0};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_positional_option(args.cmd, sharg::config{.description = "command options: build, query"});
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.q, sharg::config{.short_id = 'q', .long_id = "query", .description = "provide query file"});
    parser.add_option(args.d, sharg::config{.short_id = 'd', .long_id = "dict", .description = "provide dict file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
    parser.add_option(args.m1, sharg::config{.long_id = "m1", .description = "minimiser1 length"});
    parser.add_option(args.m2, sharg::config{.long_id = "m2", .description = "minimiser2 length"});
    parser.add_option(args.m3, sharg::config{.long_id = "m3", .description = "minimiser3 length"});
    parser.add_option(args.t1, sharg::config{.long_id = "t1", .description = "threshold1"});
    parser.add_option(args.t2, sharg::config{.long_id = "t2", .description = "threshold2"});
    parser.add_option(args.t3, sharg::config{.long_id = "t3", .description = "threshold3"});
    parser.add_option(args.s, sharg::config{.short_id = 's', .long_id = "span", .description = "span"});
    parser.add_flag(args.c, sharg::config{.short_id = 'c', .long_id = "comp", .description = "compress level 2 and 3"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('d'))
        throw sharg::user_input_error("provide index file.");
    if(args.cmd == "build") {
        if(!parser.is_option_set('i'))
            throw sharg::user_input_error("provide input file.");
        if(!parser.is_option_set('k'))
            throw sharg::user_input_error("specify k");
    }
    else if(args.cmd == "query") {
        if(!parser.is_option_set('q'))
            throw sharg::user_input_error("provide query file.");
    }
    else if(args.cmd == "lookup") {
        
    }
    else
        throw sharg::user_input_error("illegal command");
    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    size_t N = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
    }
}


int main(int argc, char** argv)
{
    sharg::parser parser{"rsindex", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
        check_arguments(parser, args);
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    if(args.cmd == "build") {
        std::cout << "loading text...\n";
        std::vector<std::vector<seqan3::dna4>> text;
        load_file(args.i, text);

        std::cout << "building dict...\n";
        if(args.c) {
            RSIndexComp index = RSIndexComp(args.k, args.m1, args.m2, args.m3, args.t1, args.t2, args.t3, args.s);
            index.build(text);
            index.save(args.d);            
        }
        else {
            RSIndex index = RSIndex(args.k, args.m1, args.m2, args.m3, args.t1, args.t2, args.t3, args.s);
            index.build(text);
            index.save(args.d);
        }
    }
    else if(args.cmd == "query") {
        std::cout << "loading queries...\n";
        std::vector<std::vector<seqan3::dna4>> queries;
        load_file(args.q, queries);

        std::cout << "loading dict...\n";
        RSIndexComp index = RSIndexComp();
        index.load(args.d);
        std::cout << "querying...\n";

        uint64_t kmers = 0;
        uint64_t found = 0;
        uint64_t extensions = 0;

        std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
        for (auto query : queries) {
            found += index.streaming_query(query, extensions);
            kmers += query.size() - index.getk() + 1;
        }
        std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        double ns_per_kmer = (double) elapsed.count() / kmers;
        
        std::cout << "==== query report:\n";
        std::cout << "num_kmers = " << kmers << '\n';
        std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers*100 << "%)\n";
        std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
        std::cout << "num extensions = " << extensions << '\n';
    }
    else if(args.cmd == "lookup") {
        uint64_t found = 0;
        double ns_per_kmer;
        const int rounds = 5;
        std::vector<uint64_t> kmers;
        std::chrono::high_resolution_clock::time_point t_start, t_stop;

        std::cout << "loading dict...\n";
        if(args.c) {
            RSIndexComp index = RSIndexComp();
            index.load(args.d);
            kmers = index.rand_text_kmers(1000000);
            std::cout << "bench lookup...\n";

            
        }
        else {
            RSIndex index = RSIndex();
            index.load(args.d);
            index.stats();
            kmers = index.rand_text_kmers(1000000);
            std::cout << "bench lookup...\n";

            t_start = std::chrono::high_resolution_clock::now();
            // for(int r = 0; r < rounds; r++) {
            //     found = index.lookup(kmers);
            // }
            found = index.lookup(kmers);
            t_stop = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
            // ns_per_kmer = (double) elapsed.count() / (kmers.size() * rounds);
            ns_per_kmer = (double) elapsed.count() / kmers.size();
            std::cout << "==== positive lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "pos_time_per_kmer = " << ns_per_kmer << '\n';

            kmers = rand_kmers(1000000, index.getk());
            std::cout << "bench lookup...\n";

            t_start = std::chrono::high_resolution_clock::now();
            // for(int r = 0; r < rounds; r++) {
            //     found = index.lookup(kmers);
            // }
            found = index.lookup(kmers);
            t_stop = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
            // ns_per_kmer = (double) elapsed.count() / (kmers.size() * rounds);
            ns_per_kmer = (double) elapsed.count() / kmers.size();
            std::cout << "==== negative lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_negative_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "neg_time_per_kmer = " << ns_per_kmer << '\n';
        }
    }
 
    return 0;
}

