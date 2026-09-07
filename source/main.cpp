#include <sharg/all.hpp>
#include "rshash.hpp"


std::vector<uint64_t> rand_kmers(const uint64_t n, const uint64_t k)
{
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
    uint8_t level{2};
    uint8_t m1{18};
    uint8_t m2{21};
    uint8_t m3{23};
    uint16_t t1{64};
    uint16_t t2{64};
    uint16_t t3{64};
    uint16_t t{0};
    bool loc{false};
    bool ht{false};
    std::vector<uint32_t> shapes{std::numeric_limits<uint32_t>::max()};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args)
{
    parser.add_positional_option(args.cmd, sharg::config{.description = "command options: build, lookup, locate, bench"});
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.q, sharg::config{.short_id = 'q', .long_id = "query", .description = "provide query file"});
    parser.add_option(args.d, sharg::config{.short_id = 'd', .long_id = "dict", .description = "provide dict file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
    parser.add_option(args.level, sharg::config{.long_id = "level", .description = "no level"});
    parser.add_option(args.m1, sharg::config{.long_id = "m1", .description = "minimiser1 length"});
    parser.add_option(args.m2, sharg::config{.long_id = "m2", .description = "minimiser2 length"});
    parser.add_option(args.m3, sharg::config{.long_id = "m3", .description = "minimiser3 length"});
    parser.add_option(args.t1, sharg::config{.long_id = "t1", .description = "threshold1"});
    parser.add_option(args.t2, sharg::config{.long_id = "t2", .description = "threshold2"});
    parser.add_option(args.t3, sharg::config{.long_id = "t3", .description = "threshold3"});
    parser.add_option(args.t, sharg::config{.short_id = 't', .description = "max k-mer/shape frequency threshold"});
    parser.add_flag(args.loc, sharg::config{.long_id = "loc", .description = "enable locate"});
    parser.add_flag(args.ht, sharg::config{.long_id = "ht", .description = "do not use hashtable on last level"});
    parser.add_option(args.shapes, sharg::config{.long_id = "shapes", .description = "list of shape values"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('d'))
        throw sharg::user_input_error("provide index file.");
    if(args.cmd == "build") {
        if(!parser.is_option_set('i'))
            throw sharg::user_input_error("provide input file.");
    }
    else if(args.cmd == "lookup") {
        if(!parser.is_option_set('q'))
            throw sharg::user_input_error("provide query file.");
    }
    else if(args.cmd == "locate") {
        if(!parser.is_option_set('q'))
            throw sharg::user_input_error("provide query file.");
    }
    else if(args.cmd != "bench")
        throw sharg::user_input_error("illegal command");
    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


size_t load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    size_t max_length = 0;
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
        max_length = std::max(max_length, record.sequence().size());
    }
    return max_length;
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
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }

    if(args.cmd == "build")
    {
        std::cout << "loading text...\n";
        std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
        load_file(args.i, text);

        std::cout << "building dict...\n";
        RSHash index;
        if(args.shapes[0] == std::numeric_limits<uint32_t>::max())
            index = RSHash(args.k, args.level, args.m1, args.m2, args.m3, args.t1, args.t2, args.t3, args.t, args.loc, !args.ht);
        else {
            const Shapes32 shapes = shape32_create(args.shapes);
            index = RSHash(shapes, args.level, args.m1, args.m2, args.m3, args.t1, args.t2, args.t3, args.t, args.loc, !args.ht);
        }
        index.build(text);
        index.save(args.d);
    }
    else if(args.cmd == "lookup")
    {
        std::cout << "loading queries...\n";
        std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> queries;
        load_file(args.q, queries);

        uint64_t kmers = 0;
        uint64_t found = 0;
        uint64_t extensions = 0;
        std::chrono::nanoseconds elapsed;
        RSHash index = RSHash();
        index.load(args.d);
        std::cout << "querying...\n";

        std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
        for (auto query : queries) {
            found += index.streaming_lookup(query, extensions);
            if(index.getshapes().shapes[0].value != std::numeric_limits<uint32_t>::max())
                kmers += query.size() - index.getshapes().length + 1;
            else
                kmers += query.size() - index.getk() + 1;
        }
        std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        
        double ns_per_kmer = (double) elapsed.count() / kmers;
        
        std::cout << "==== query report:\n";
        std::cout << "num_kmers = " << kmers << '\n';
        std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers*100 << "%)\n";
        std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
        std::cout << "extensions = " << extensions << '\n';
    }
    else if(args.cmd == "locate")
    {
        std::cout << "loading queries...\n";
        std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> queries;
        const size_t max_query_length = load_file(args.q, queries);

        std::chrono::nanoseconds elapsed;
        size_t kmers = 0;

        std::cout << "loading dict...\n";
        RSHash index = RSHash();
        index.load(args.d);

        if(!index.has_locate()) {
            std::cerr << "index does not support locate\n";
            return -1;
        }

        std::cout << "locating...\n";
        uint64_t found_positions = 0;
        uint64_t found_kmers = 0;
        // const size_t max_kmer_freq = 1000000;
        std::vector<uint64_t> positions;
        // positions.reserve(max_query_length*max_kmer_freq);
        std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
        for (auto query : queries) {
            found_kmers += index.streaming_locate(query, positions, found_positions);
            found_positions += positions.size();
            if(index.getshapes().shapes[0].value != std::numeric_limits<uint32_t>::max())
                kmers += query.size() - index.getshapes().length + 1;
            else
                kmers += query.size() - index.getk() + 1;
            // positions.clear();
        }
        std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        
        double ns_per_kmer = (double) elapsed.count() / kmers;
        
        std::cout << "==== query report:\n";
        std::cout << "num_kmers = " << kmers << '\n';
        std::cout << "num_positive_kmers = " << found_kmers << "\n";
        std::cout << "num_positions = " << found_positions << '\n';
        std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
    }
    else if(args.cmd == "bench")
    {
        std::cout << "loading dict...\n";
        uint64_t found = 0;
        double ns_per_kmer;
        const int rounds = 5;
        std::vector<uint64_t> kmers;

        RSHash index = RSHash();
        index.load(args.d);

        std::cout << "bench pos lookup...\n";
        double error = 0.0;
        std::vector<double> times;
        double lookup_time_sum = 0.0;
        int round = 0;
        std::chrono::high_resolution_clock::time_point t_start, t_stop;
        std::chrono::nanoseconds elapsed;

        while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
            kmers = index.rand_text_kmers(1000000);
            t_start = std::chrono::high_resolution_clock::now();
            found = index.lookup(kmers);
            t_stop = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
            ns_per_kmer = (double) elapsed.count() / kmers.size();
            lookup_time_sum += ns_per_kmer;
            round++;
            error += std::abs((lookup_time_sum/round) - ns_per_kmer);
            times.push_back(ns_per_kmer);
            std::cout << "round " << round << " found " << found << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
        }

        std::cout << "==== positive lookup:\n";
        std::cout << "num_kmers = " << kmers.size() << '\n';
        std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
        // std::cout << "pos_time_per_kmer = " << lookup_time_sum/round << '\n';
        std::cout << "pos_time_per_kmer = " << std::accumulate(times.begin(), times.end(), 0.0) / times.size() << '\n';
        std::cout << "pos_time_per_kmer_variance = " << std::sqrt(std::accumulate(times.begin(), times.end(), 0.0, [&](double acc, double x) { return acc + (x - (lookup_time_sum/round)) * (x - (lookup_time_sum/round)); }) / times.size()) << '\n';

        std::cout << "bench neg lookup...\n";
        round = 0;
        error = 0.0;
        lookup_time_sum = 0.0;
        times.clear();

        while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
            kmers = rand_kmers(1000000, index.getk());
            t_start = std::chrono::high_resolution_clock::now();
            found = index.lookup(kmers);
            t_stop = std::chrono::high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
            ns_per_kmer = (double) elapsed.count() / kmers.size();
            lookup_time_sum += ns_per_kmer;
            round++;
            error += std::abs((lookup_time_sum/round) - ns_per_kmer);
            times.push_back(ns_per_kmer);
            std::cout << "round " << round << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
        }

        std::cout << "==== negative lookup:\n";
        std::cout << "num_kmers = " << kmers.size() << '\n';
        std::cout << "num_negative_kmers = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
        // std::cout << "neg_time_per_kmer = " << ns_per_kmer << '\n';
        std::cout << "neg_time_per_kmer = " << std::accumulate(times.begin(), times.end(), 0.0) / times.size() << '\n';
        std::cout << "neg_time_per_kmer_variance = " << std::sqrt(std::accumulate(times.begin(), times.end(), 0.0, [&](double acc, double x) { return acc + (x - (lookup_time_sum/round)) * (x - (lookup_time_sum/round)) ; }) / times.size()) << '\n';

        if(index.has_locate()) {
            std::cout << "bench locate...\n";

            std::cout << "bench pos locate...\n";
            round = 0;
            error = 0.0;
            lookup_time_sum = 0.0;
            times.clear();
            std::vector<uint64_t> positions;
            // positions.reserve(100000);

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = index.rand_text_kmers(1000000);
                t_start = std::chrono::high_resolution_clock::now();
                found = index.locate(kmers, positions);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                times.push_back(ns_per_kmer);
                std::cout << "round " << round << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== locate:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_positive_positions = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "pos_locate_time_per_kmer = " << std::accumulate(times.begin(), times.end(), 0.0) / times.size() << '\n';
            std::cout << "pos_locate_time_per_kmer_variance = " << std::sqrt(std::accumulate(times.begin(), times.end(), 0.0, [&](double acc, double x) { return acc + (x - (lookup_time_sum/round)) * (x - (lookup_time_sum/round)); }) / times.size()) << '\n';


            std::cout << "bench neg locate...\n";
            round = 0;
            error = 0.0;
            lookup_time_sum = 0.0;
            times.clear();

            while((round < 10 || error/round > 0.05 * (lookup_time_sum/round)) && round <= 50) {
                kmers = rand_kmers(1000000, index.getk());
                t_start = std::chrono::high_resolution_clock::now();
                found = index.locate(kmers, positions);
                t_stop = std::chrono::high_resolution_clock::now();
                elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
                ns_per_kmer = (double) elapsed.count() / kmers.size();
                lookup_time_sum += ns_per_kmer;
                round++;
                error += std::abs((lookup_time_sum/round) - ns_per_kmer);
                times.push_back(ns_per_kmer);
                std::cout << "round " << round << " time per kmer: " << ns_per_kmer << ", avg: " << (lookup_time_sum/round) << ", error: " << error/round << '\n';
            }

            std::cout << "==== locate:\n";
            std::cout << "num_kmers = " << kmers.size() << '\n';
            std::cout << "num_negative_positions = " << found << " (" << (double) found/kmers.size()*100 << "%)\n";
            std::cout << "neg_locate_time_per_kmer = " << std::accumulate(times.begin(), times.end(), 0.0) / times.size() << '\n';
            std::cout << "neg_locate_time_per_kmer_variance = " << std::sqrt(std::accumulate(times.begin(), times.end(), 0.0, [&](double acc, double x) { return acc + (x - (lookup_time_sum/round)) * (x - (lookup_time_sum/round)); }) / times.size()) << '\n';

        }
    }
 
    return 0;
}

