#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
// #include "../source/minimiser_views.hpp"
#include "../source/kmer_view.hpp"
#include "../source/shape.hpp"
#include "../source/flat_map.hpp"
#include <gtl/phmap.hpp>
#include <fstream>
#include <numeric>
#include <random>


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::vector<std::filesystem::path> q{};
    std::filesystem::path o{};
    uint8_t k{31};
    std::vector<uint64_t> shapes{std::numeric_limits<uint64_t>::max()};
    bool locate{false};
    bool bench{false};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.q, sharg::config{.short_id = 'q', .long_id = "query", .description = "provide query file(s); the table is built once for all of them"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
    parser.add_option(args.shapes, sharg::config{.long_id = "shapes", .description = "list of shape values"});
    parser.add_flag(args.locate, sharg::config{.long_id = "locate", .description = "locate with the FlatMap of the rshash last level (key -> positions) instead of lookup with a hash set"});
    parser.add_flag(args.bench, sharg::config{.long_id = "bench", .description = "benchmark single random text k-mers (pos) and random k-mers (neg) like rshash bench"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw sharg::user_input_error("provide input file.");
    if(!parser.is_option_set('q') && !args.bench)
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

// peak resident set size of the process in kB (like "Maximum resident set size" of /usr/bin/time -v)
uint64_t peak_rss_kb()
{
    std::ifstream status("/proc/self/status");
    std::string line;
    while(std::getline(status, line))
        if(line.rfind("VmHWM:", 0) == 0)
            return std::stoull(line.substr(6));
    return 0;
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
        std::cerr << "[Error] " << ext.what() << '\n';
        return -1;
    }

    std::chrono::high_resolution_clock::time_point t_build_start = std::chrono::high_resolution_clock::now();

    std::cout << "loading text...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
    load_file(args.i, text);

    // start of each sequence in the concatenated text, the last entry is the text length
    std::vector<uint64_t> seq_offsets(text.size() + 1, 0);
    for(size_t s = 0; s < text.size(); s++)
        seq_offsets[s+1] = seq_offsets[s] + text[s].size();
    if(args.locate && seq_offsets.back() > std::numeric_limits<uint32_t>::max()) {
        std::cerr << "text too long for 32 bit positions\n";
        return -1;
    }

    std::cout << "building " << (args.locate ? "hashmap" : "hashtable") << "...\n";
    const bool use_shapes = args.shapes[0] != std::numeric_limits<uint64_t>::max();

    Shapes64 shapes;
    if(use_shapes)
        shapes = shape64_create(args.shapes);
    const size_t no_shapes = use_shapes ? shapes.shapes.size() : 0;
    // fwd and rev structure per shape, one canonical structure for k-mers
    const size_t no_structures = use_shapes ? 2*no_shapes : 1;
    const size_t window_size = use_shapes ? shapes.length : args.k;

    std::vector<gtl::flat_hash_set<uint64_t>> hts(args.locate ? 0 : no_structures);
    std::vector<FlatMap> hms(args.locate ? no_structures : 0);

    // calls f(key, text position) for all text k-mers/shapes of structure j: the forward (even j) or reverse
    // (odd j) values of shape j/2, for k-mers the canonical k-mers (j = 0)
    auto for_each_text_key = [&](const size_t j, auto &&f) {
        if(use_shapes) {
            const Shape64 shape = shapes.shapes[j / 2];
            const bool fwd = j % 2 == 0;
            for(size_t s = 0; s < text.size(); s++) {
                uint64_t pos = seq_offsets[s];
                for(auto && window : text[s] | rshash::views::kmer_view({.window_size = shape.length})) {
                    if(fwd)
                        f(_pext_u64(window.value, shape.mask.lo) | (_pext_u64(window.value_hi, shape.mask.hi) << shape.lo_weight), pos);
                    else
                        f(_pext_u64(window.value, shape.mask_rev.lo) | (_pext_u64(window.value_hi, shape.mask_rev.hi) << shape.rev_lo_weight), pos);
                    pos++;
                }
            }
        }
        else {
            for(size_t s = 0; s < text.size(); s++) {
                uint64_t pos = seq_offsets[s];
                for(auto && window : text[s] | rshash::views::kmer_view({.window_size = args.k})) {
                    f(std::min<uint64_t>(window.value, window.value_rev), pos);
                    pos++;
                }
            }
        }
    };

    for(size_t j = 0; j < no_structures; j++) {
        if(args.locate) {
            // as on the rshash last level: the (key, position) pairs are sorted instead of hashed one by one
            std::vector<KeyValue> pairs;
            for_each_text_key(j, [&](const uint64_t key, const uint64_t pos) { pairs.push_back({key, static_cast<uint32_t>(pos)}); });
            sort_unique(pairs);
            hms[j] = FlatMap::from_sorted(pairs);
        }
        else
            for_each_text_key(j, [&](const uint64_t key, const uint64_t) { hts[j].insert(key); });
    }

    std::chrono::high_resolution_clock::time_point t_build_stop = std::chrono::high_resolution_clock::now();
    const double build_time = std::chrono::duration_cast<std::chrono::duration<double>>(t_build_stop - t_build_start).count();
    const uint64_t build_mem = peak_rss_kb();

    // space report (bits per text k-mer, i.e. per window of the shapes)
    uint64_t no_text_kmers = 0;
    for(auto & sequence : text)
        no_text_kmers += sequence.size() >= window_size ? sequence.size() - window_size + 1 : 0;
    uint64_t space = 0, keys = 0;
    for(auto &ht : hts) {
        space += ht.capacity()*(sizeof(uint64_t) + 1)*8;
        keys += ht.size();
    }
    for(auto &hm : hms) {
        space += hm.memory_bits();
        keys += hm.size();
    }
    std::cout << "textkmers: " << no_text_kmers << '\n';
    std::cout << "no structures: " << no_structures << '\n';
    std::cout << "no keys: " << keys << '\n';
    std::cout << "Last_level: " << (double) space/no_text_kmers << '\n';
    std::cout << "total: " << (double) space/no_text_kmers << '\n';
    std::cout << "build_time = " << build_time << '\n';
    std::cout << "build_mem = " << build_mem << '\n';

    // lookup: whether any shape (fwd or rev) of the window is found
    auto lookup_window = [&](const rshash::kmerview_result &window) {
        if(use_shapes) {
            for(size_t i = 0; i < no_shapes; i++) {
                const Shape64 &shape = shapes.shapes[i];
                const uint64_t kmer_fwd = _pext_u64(window.value, shape.w_mask.lo) | (_pext_u64(window.value_hi, shape.w_mask.hi) << shape.w_lo_weight);
                const uint64_t kmer_rev = _pext_u64(window.value_rev, shape.w_mask_rev.lo) | (_pext_u64(window.value_rev_hi, shape.w_mask_rev.hi) << shape.w_rev_lo_weight);
                if(hts[2*i].contains(kmer_fwd) || hts[2*i+1].contains(kmer_rev))
                    return true;
            }
            return false;
        }
        return hts[0].contains(std::min<uint64_t>(window.value, window.value_rev));
    };

    // locate: number of positions of all shapes (fwd and rev) of the window; like rshash locate_last_level
    // the positions are not read, only counted
    auto locate_window = [&](const rshash::kmerview_result &window) {
        if(use_shapes) {
            uint64_t no_positions = 0;
            for(size_t i = 0; i < no_shapes; i++) {
                const Shape64 &shape = shapes.shapes[i];
                const uint64_t kmer_fwd = _pext_u64(window.value, shape.w_mask.lo) | (_pext_u64(window.value_hi, shape.w_mask.hi) << shape.w_lo_weight);
                const uint64_t kmer_rev = _pext_u64(window.value_rev, shape.w_mask_rev.lo) | (_pext_u64(window.value_rev_hi, shape.w_mask_rev.hi) << shape.w_rev_lo_weight);
                no_positions += hms[2*i].find(kmer_fwd).size();
                no_positions += hms[2*i+1].find(kmer_rev).size();
            }
            return (uint64_t) no_positions;
        }
        return (uint64_t) hms[0].find(std::min<uint64_t>(window.value, window.value_rev)).size();
    };

    if(args.bench) {
        const uint64_t n = 1000000;
        std::mt19937_64 rng(1);

        // random text windows (pos)
        std::vector<rshash::kmerview_result> pos_kmers;
        pos_kmers.reserve(n);
        std::uniform_int_distribution<uint64_t> distr(0, seq_offsets.back() - 1);
        while(pos_kmers.size() < n) {
            const uint64_t offset = distr(rng);
            const size_t s = std::upper_bound(seq_offsets.begin(), seq_offsets.end(), offset) - seq_offsets.begin() - 1;
            const uint64_t p = offset - seq_offsets[s];
            if(p + window_size > text[s].size())
                continue;
            for(auto && window : text[s] | std::views::drop(p) | std::views::take(window_size) | rshash::views::kmer_view({.window_size = window_size}))
                pos_kmers.push_back(window);
        }

        // random k-mers (neg): the windows of a random sequence
        std::vector<rshash::kmerview_result> neg_kmers;
        neg_kmers.reserve(n);
        seqan3::bitpacked_sequence<seqan3::dna4> random_sequence;
        std::uniform_int_distribution<int> base(0, 3);
        for(size_t i = 0; i < n + window_size - 1; i++)
            random_sequence.push_back(seqan3::dna4{}.assign_rank(base(rng)));
        for(auto && window : random_sequence | rshash::views::kmer_view({.window_size = window_size}))
            neg_kmers.push_back(window);

        // like rshash bench: repeat until the average time per k-mer is stable
        auto bench = [&](const std::vector<rshash::kmerview_result> &kmers, auto &&query, uint64_t &found, uint64_t &no_positions) {
            std::vector<double> times;
            double time_sum = 0.0, error = 0.0;
            int round = 0;
            while((round < 10 || error/round > 0.05 * (time_sum/round)) && round <= 50) {
                found = 0;
                no_positions = 0;
                std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
                for(const auto &kmer : kmers) {
                    const uint64_t r = query(kmer);
                    found += r > 0;
                    no_positions += r;
                }
                std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
                const double ns_per_kmer = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start).count() / kmers.size();
                time_sum += ns_per_kmer;
                round++;
                error += std::abs((time_sum/round) - ns_per_kmer);
                times.push_back(ns_per_kmer);
            }
            return std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        };

        uint64_t found, no_positions;
        if(!args.locate) {
            for(const bool positive : {true, false}) {
                const double t = bench(positive ? pos_kmers : neg_kmers, [&](const auto &kmer) { return (uint64_t) lookup_window(kmer); }, found, no_positions);
                std::cout << "==== " << (positive ? "positive" : "negative") << " lookup:\n";
                std::cout << "num_kmers = " << n << '\n';
                std::cout << "num_" << (positive ? "positive" : "negative") << "_kmers = " << found << " (" << (double) found/n*100 << "%)\n";
                std::cout << (positive ? "pos" : "neg") << "_time_per_kmer = " << t << '\n';
            }
        }
        else {
            for(const bool positive : {true, false}) {
                const double t = bench(positive ? pos_kmers : neg_kmers, locate_window, found, no_positions);
                std::cout << "==== " << (positive ? "positive" : "negative") << " locate:\n";
                std::cout << "num_kmers = " << n << '\n';
                std::cout << "num_" << (positive ? "positive" : "negative") << "_kmers = " << found << " (" << (double) found/n*100 << "%)\n";
                std::cout << "num_positions = " << no_positions << '\n';
                std::cout << (positive ? "pos" : "neg") << "_locate_time_per_kmer = " << t << '\n';
            }
        }
    }

    for(const auto &query_file : args.q) {
        std::cout << "loading queries...\n";
        std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> queries;
        load_file(query_file, queries);

        std::cout << "querying...\n";
        uint64_t kmers = 0;
        uint64_t found_kmers = 0;
        uint64_t found_positions = 0;
        uint64_t found_queries = 0;

        std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
        for (auto& query : queries) {
            bool found_query = false;
            for (auto&& window : query | rshash::views::kmer_view({.window_size = window_size})) {
                bool found;
                if(args.locate) {
                    const uint64_t no_positions = locate_window(window);
                    found_positions += no_positions;
                    found = no_positions > 0;
                }
                else
                    found = lookup_window(window);
                found_query |= found;
                found_kmers += found;
                kmers++;
            }
            found_queries += found_query;
        }
        std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);

        double ns_per_kmer = (double) elapsed.count() / kmers;
        double ns_per_read = (double) elapsed.count() / queries.size();

        std::cout << "==== query report: " << query_file.string() << '\n';
        std::cout << "num_kmers = " << kmers << '\n';
        std::cout << "num_reads = " << queries.size() << '\n';
        std::cout << "num_positive_kmers = " << found_kmers << " (" << (double) found_kmers/kmers*100 << "%)\n";
        if(args.locate)
            std::cout << "num_positions = " << found_positions << '\n';
        std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
        std::cout << "time_per_read = " << ns_per_read << '\n';
        std::cout << "num_found_queries = " << found_queries << " (" << (double) found_queries/queries.size()*100 << "%)\n";
        std::cout << "query_mem = " << peak_rss_kb() << '\n';
        std::cout << "==== end query report\n";
    }

    return 0;
}
