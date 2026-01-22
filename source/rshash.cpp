#include <seqan3/core/debug_stream.hpp>
#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "rshash.hpp"


std::vector<uint64_t> rand_kmers(const uint64_t n, const uint64_t k) {
    const uint64_t mask = compute_mask(2u * k);

    std::uniform_int_distribution<uint64_t> distr;
    std::mt19937_64 m_rand(1);
    std::vector<uint64_t> kmers;
    kmers.reserve(n);

    for (uint64_t i = 0; i < n; ++i) {
        const uint64_t kmer = distr(m_rand) & mask;
        kmers.push_back(kmer);
    }

    return kmers;
}


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path d{};
    uint8_t k{31};
    uint8_t l{2};
    uint8_t m1{16};
    uint8_t m2{18};
    uint8_t m3{20};
    uint8_t t1{65};
    uint8_t t2{65};
    uint16_t t3{65};
    bool c{false};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_positional_option(args.cmd, sharg::config{.description = "command options: build, query"});
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.q, sharg::config{.short_id = 'q', .long_id = "query", .description = "provide query file"});
    parser.add_option(args.d, sharg::config{.short_id = 'd', .long_id = "dict", .description = "provide dict file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
    parser.add_option(args.l, sharg::config{.short_id = 'l', .long_id = "level", .description = "no level"});
    parser.add_option(args.m1, sharg::config{.long_id = "m1", .description = "minimiser1 length"});
    parser.add_option(args.m2, sharg::config{.long_id = "m2", .description = "minimiser2 length"});
    parser.add_option(args.m3, sharg::config{.long_id = "m3", .description = "minimiser3 length"});
    parser.add_option(args.t1, sharg::config{.long_id = "t1", .description = "threshold1"});
    parser.add_option(args.t2, sharg::config{.long_id = "t2", .description = "threshold2"});
    parser.add_option(args.t3, sharg::config{.long_id = "t3", .description = "threshold3"});
    parser.add_flag(args.c, sharg::config{.short_id = 'c', .long_id = "comp", .description = "compress first level"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('d'))
        throw sharg::user_input_error("provide index file.");
    if(args.cmd == "build") {
        if(!parser.is_option_set('i'))
            throw sharg::user_input_error("provide input file.");
        if(!parser.is_option_set('k'))
            throw sharg::user_input_error("specify k");
        if(!parser.is_option_set('l'))
            throw sharg::user_input_error("specify level");
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
        // if(N > 1'000'000'000)
        //     break;
    }
}

uint64_t kmer_to_int(std::vector<seqan3::dna4> &kmerdna4, const uint8_t k) {
    uint64_t kmer = 0;
    for (uint8_t j=0; j < k; j++) {
        uint64_t const new_rank = seqan3::to_rank(kmerdna4[j]);
        kmer = (kmer >> 2) | (new_rank << 2*(k-1));
    }
    return kmer;
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
        if(args.l == 1) {
            RSHash1 index = RSHash1(args.k, args.m1, args.t1);
            index.build(text);
            index.save(args.d);
        }
        else if(args.l == 2) {
            RSHash2 index = RSHash2(args.k, args.m1, args.m2, args.t1, args.t2);
            index.build(text);
            index.save(args.d);
        }
        else if(args.l == 3) {
            if(args.c) {
                RSHash3C index = RSHash3C(args.k, args.m1, args.m2, args.m3, args.t1, args.t2, args.t3);
                index.build(text);
                index.save(args.d);
                return 0;
            }
            else {
                RSHash3 index = RSHash3(args.k, args.m1, args.m2, args.m3, args.t1, args.t2, args.t3);
                index.build(text);
                index.save(args.d);
                return 0;
            }
        }
    }
    else if(args.cmd == "query") {
        std::cout << "loading queries...\n";
        std::vector<std::vector<seqan3::dna4>> queries;
        load_file(args.q, queries);

        uint64_t kmers = 0;
        uint64_t found = 0;
        uint64_t extensions = 0;
        std::chrono::nanoseconds elapsed;
        if(args.l == 1) {
            std::cout << "loading dict...\n";
            RSHash1 index = RSHash1();
            index.load(args.d);
            std::cout << "querying...\n";

            std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
            for (auto query : queries) {
                found += index.streaming_query(query, extensions);
                kmers += query.size() - index.getk() + 1;
            }
            std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        
        }
        else if(args.l == 2) {
            std::cout << "loading dict...\n";
            RSHash2 index = RSHash2();
            index.load(args.d);
            std::cout << "querying...\n";

            std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
            for (auto query : queries) {
                found += index.streaming_query(query, extensions);
                kmers += query.size() - index.getk() + 1;
            }
            std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        }
        else if(args.l == 3) {
            if(args.c) {
                std::cout << "loading dict...\n";
                RSHash3C index = RSHash3C();
                index.load(args.d);
                std::cout << "querying...\n";

                std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
                for (auto query : queries) {
                    found += index.streaming_query(query, extensions);
                    kmers += query.size() - index.getk() + 1;
                }
                std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
            }
            else {
                std::cout << "loading dict...\n";
                RSHash3 index = RSHash3();
                index.load(args.d);
                std::cout << "querying...\n";

                std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
                for (auto query : queries) {
                    found += index.streaming_query(query, extensions);
                    kmers += query.size() - index.getk() + 1;
                }
                std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
            }
        }
        double ns_per_kmer = (double) elapsed.count() / kmers;
        
        std::cout << "==== query report:\n";
        std::cout << "num_kmers = " << kmers << '\n';
        std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers*100 << "%)\n";
        std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
        std::cout << "extensions = " << extensions << '\n';
    }
    else if(args.cmd == "lookup") {
        std::cout << "loading dict...\n";
        uint64_t found = 0;
        double ns_per_kmer;
        const int rounds = 5;
        const bool verbose = false;

        if(args.l == 1) {
            RSHash1 index = RSHash1();
            index.load(args.d);
            std::vector<uint64_t> kmers;

            std::cout << "bench lookup...\n";
            double error = 0.0;
            double lookup_time_sum = 0.0;
            int round = 0;
            std::chrono::high_resolution_clock::time_point t_start, t_stop;
            std::chrono::nanoseconds elapsed;

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = index.rand_text_kmers(1000000);
                t_start = std::chrono::high_resolution_clock::now();
                found = index.lookup(kmers, verbose);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                std::cout << "round " << round << " found " << found << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== positive lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "pos_time_per_kmer = " << lookup_time_sum/round << '\n';

            std::cout << "bench lookup...\n";

            round = 0;
            error = 0.0;
            lookup_time_sum = 0.0;

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = rand_kmers(1000000, index.getk());
                t_start = std::chrono::high_resolution_clock::now();
                found = index.lookup(kmers, verbose);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                std::cout << "round " << round << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== negative lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_negative_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "neg_time_per_kmer = " << ns_per_kmer << '\n';
        }
        else if(args.l == 2) {
            RSHash2 index = RSHash2();
            index.load(args.d);
            std::vector<uint64_t> kmers;

            std::cout << "bench lookup...\n";
            double error = 0.0;
            double lookup_time_sum = 0.0;
            int round = 0;
            std::chrono::high_resolution_clock::time_point t_start, t_stop;
            std::chrono::nanoseconds elapsed;

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = index.rand_text_kmers(1000000);
                t_start = std::chrono::high_resolution_clock::now();
                found = index.lookup(kmers, verbose);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                std::cout << "round " << round << " found " << found << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== positive lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "pos_time_per_kmer = " << lookup_time_sum/round << '\n';

            std::cout << "bench lookup...\n";

            round = 0;
            error = 0.0;
            lookup_time_sum = 0.0;

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = rand_kmers(1000000, index.getk());
                t_start = std::chrono::high_resolution_clock::now();
                found = index.lookup(kmers, verbose);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                std::cout << "round " << round << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== negative lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_negative_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "neg_time_per_kmer = " << ns_per_kmer << '\n';
        }
        else if(args.l == 3) {
            RSHash3C index = RSHash3C();
            index.load(args.d);
            std::vector<uint64_t> kmers;

            std::cout << "bench lookup...\n";
            double error = 0.0;
            double lookup_time_sum = 0.0;
            int round = 0;
            std::chrono::high_resolution_clock::time_point t_start, t_stop;
            std::chrono::nanoseconds elapsed;

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = index.rand_text_kmers(1000000);
                t_start = std::chrono::high_resolution_clock::now();
                found = index.lookup(kmers, verbose);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                std::cout << "round " << round << " found " << found << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== positive lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "pos_time_per_kmer = " << lookup_time_sum/round << '\n';

            std::cout << "bench lookup...\n";

            round = 0;
            error = 0.0;
            lookup_time_sum = 0.0;

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = rand_kmers(1000000, index.getk());
                t_start = std::chrono::high_resolution_clock::now();
                found = index.lookup(kmers, verbose);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                std::cout << "round " << round << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== negative lookup:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_negative_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "neg_time_per_kmer = " << ns_per_kmer << '\n';

        }
    }
 
    return 0;
}

