#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
// #include "../source/minimiser_views.hpp"
#include "../source/kmer_view.hpp"
#include "../source/shape.hpp"
#include <gtl/phmap.hpp>


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path o{};
    uint8_t k{31};
    std::vector<uint64_t> shapes{std::numeric_limits<uint64_t>::max()};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.q, sharg::config{.short_id = 'q', .long_id = "query", .description = "provide query file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
    parser.add_option(args.shapes, sharg::config{.long_id = "shapes", .description = "list of shape values"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw sharg::user_input_error("provide input file.");
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
    const bool use_shapes = args.shapes[0] != std::numeric_limits<uint64_t>::max();

    std::vector<gtl::flat_hash_set<uint64_t>> hts(1);

    Shapes64 shapes;
    size_t canonical_no_shapes;
    if(use_shapes) {
        shapes = shape64_create(args.shapes);
        canonical_no_shapes = std::count_if(shapes.shapes.begin(), shapes.shapes.end(), [](const auto& x) { return x.is_canonical; });
        const size_t no_hts = 2*(shapes.shapes.size() - canonical_no_shapes) + canonical_no_shapes;
        hts.resize(no_hts);
    }

    if(use_shapes) {
        for(size_t i = 0; i < shapes.shapes.size(); ++i) {
            const Shape64 shape = shapes.shapes[i];
                for(auto & sequence : text) {
                    for(auto && window : sequence | rshash::views::kmer_view({.window_size = shape.length})) {
                        const uint64_t kmer_fwd = _pext_u64(window.value, shape.mask.lo) | (_pext_u64(window.value_hi, shape.mask.hi) << shape.lo_weight);
                        // const uint64_t kmer_rev = _pext_u64(window.value_rev, shape.mask.lo) | (_pext_u64(window.value_rev_hi, shape.mask.hi) << shape.lo_weight);

                        hts[2*i].insert(kmer_fwd);
                        // hts[2*i + 1].insert(kmer_rev);
                    }
                }
        }
        // for(auto & sequence : text) {
        //     for(auto && window : sequence | rshash::views::longkmerview({.window_size = shapes.length})) {
        //         for(size_t i = 0; i < shapes.shapes.size(); ++i) {
        //             const Shape64 shape = shapes.shapes[i];
        //             const uint64_t kmer_fwd = _pext_u64(window.value_lo, shape.w_mask.lo) | (_pext_u64(window.value_hi, shape.w_mask.hi) << shape.w_lo_weight);

        //             hts[2*i].insert(kmer_fwd);
        //         }
        //     }
        // }
    }
    else {
        for(auto & sequence : text) {
            for(auto && window : sequence | rshash::views::kmer_view({.window_size = args.k})) {
                hts[0].insert(std::min<uint64_t>(window.value, window.value_rev));
            }
        }
    }
    
     
    std::cout << "loading queries...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> queries;
    load_file(args.q, queries);

    std::cout << "querying...\n";
    uint64_t kmers = 0;
    uint64_t found_kmers = 0;
    uint64_t found_positions = 0;
    uint64_t found_queries = 0;

    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    if (use_shapes) {
        for (auto& query : queries) {
            bool found = false;
            for(size_t i = 0; i < shapes.shapes.size() && !found; i++) {
                const Shape64 shape = shapes.shapes[i];
                for (auto&& window : query | rshash::views::kmer_view({.window_size = shape.length})) {
                    const uint64_t kmer_fwd = _pext_u64(window.value, shape.mask.lo) | (_pext_u64(window.value_hi, shape.mask.hi) << shape.lo_weight);
                    if(hts[2*i].contains(kmer_fwd)) {
                        found_kmers++;
                        found = true;
                        // break;
                    }
                    kmers++;
                }
                // kmers += query.size() - shape.length + 1;
            }
            found_queries += found;
        }
        // for (auto& query : queries) {
        //     bool found = false;
        //     for (auto&& window : query | rshash::views::longkmerview({.window_size = shapes.length})) {
        //         for(size_t i = 0; i < shapes.shapes.size() && !found; i++) {
        //             const Shape64 shape = shapes.shapes[i];
        //             const uint64_t kmer_fwd = _pext_u64(window.value_lo, shape.w_mask.lo) | (_pext_u64(window.value_hi, shape.w_mask.hi) << shape.w_lo_weight);
        //             if(hts[2*i].contains(kmer_fwd)) {
        //                 found_kmers++;
        //                 found = true;
        //                 break;
        //             }
        //         }
        //         kmers++;
        //         // kmers += query.size() - shape.length + 1;
        //     }
        //     found_queries += found;
        // }
    }
    else {
        for (auto& query : queries) {
            bool query_found = false;
            for (auto&& window : query | rshash::views::kmer_view({.window_size = args.k})) {
                uint64_t kmer_value = std::min<uint64_t>(window.value, window.value_rev);
                bool found = hts[0].contains(kmer_value);
                query_found |= found;
                found_kmers += found;
                kmers++;
            }
            found_queries += query_found;
        }
    }
    
    std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        
    double ns_per_kmer = (double) elapsed.count() / kmers;
    double ns_per_read = (double) elapsed.count() / queries.size();
        
    std::cout << "==== query report:\n";
    std::cout << "num_kmers = " << kmers << '\n';
    std::cout << "num_reads = " << queries.size() << '\n';
    std::cout << "num_positive_kmers = " << found_kmers << " (" << (double) found_kmers/kmers*100 << "%)\n";
    // std::cout << "num_positions = " << found_positions << '\n';
    std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
    std::cout << "time_per_read = " << ns_per_read << '\n';
    std::cout << "num_found_queries = " << found_queries << " (" << (double) found_queries/queries.size()*100 << "%)\n";

    // report ht size
 
    return 0;
}