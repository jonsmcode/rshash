// reference for correctness checks, from the hashtable.cpp of b273cea ("longer kmers"): every shape is looked up
// on its own, the query is scanned once per shape with windows of the shape length. so every shape occurrence of
// the query is checked, independent of the aligned windows and query padding of rshash and hashtable.cpp.
// added to b273cea: the reverse structures (shape mirrored, text forward) and their lookup with the reverse
// complement of the query, one structure per shape and orientation, locate (number of text occurrences) and no
// early exit after the first shape that is found.
#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include "../source/kmer_view.hpp"
#include "../source/shape.hpp"
#include <gtl/phmap.hpp>


struct cmd_arguments {
    std::filesystem::path i{};
    std::filesystem::path q{};
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
    sharg::parser parser{"HT_reference", argc, argv};
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

    std::cout << "loading text...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
    load_file(args.i, text);

    std::cout << "building hashtable...\n";
    const bool use_shapes = args.shapes[0] != std::numeric_limits<uint64_t>::max();

    // number of text occurrences per key: forward (2i) and reverse (2i+1) structure of shape i, for k-mers one
    // structure of the canonical k-mers
    std::vector<gtl::flat_hash_map<uint64_t, uint64_t>> hts(1);

    Shapes64 shapes;
    if(use_shapes) {
        shapes = shape64_create(args.shapes);
        hts.resize(2*shapes.shapes.size());
    }

    if(use_shapes) {
        for(size_t i = 0; i < shapes.shapes.size(); ++i) {
            const Shape64 shape = shapes.shapes[i];
                for(auto & sequence : text) {
                    for(auto && window : sequence | rshash::views::kmer_view({.window_size = shape.length})) {
                        const uint64_t kmer_fwd = _pext_u64(window.value, shape.mask.lo) | (_pext_u64(window.value_hi, shape.mask.hi) << shape.lo_weight);
                        const uint64_t kmer_rev = _pext_u64(window.value, shape.mask_rev.lo) | (_pext_u64(window.value_hi, shape.mask_rev.hi) << shape.rev_lo_weight);

                        hts[2*i][kmer_fwd]++;
                        hts[2*i + 1][kmer_rev]++;
                    }
                }
        }
    }
    else {
        for(auto & sequence : text) {
            for(auto && window : sequence | rshash::views::kmer_view({.window_size = args.k})) {
                hts[0][std::min<uint64_t>(window.value, window.value_rev)]++;
            }
        }
    }

    auto occurrences = [](const gtl::flat_hash_map<uint64_t, uint64_t> &ht, const uint64_t key) -> uint64_t {
        const auto it = ht.find(key);
        return it == ht.end() ? 0 : it->second;
    };


    std::cout << "loading queries...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> queries;
    load_file(args.q, queries);

    std::cout << "querying...\n";
    uint64_t kmers = 0;          // shape occurrences (shape, orientation, query position) or k-mers of the queries
    uint64_t found_kmers = 0;    // of them found
    uint64_t found_positions = 0;
    uint64_t found_queries = 0;

    if (use_shapes) {
        for (auto& query : queries) {
            bool found = false;
            for(size_t i = 0; i < shapes.shapes.size(); i++) {
                const Shape64 shape = shapes.shapes[i];
                for (auto&& window : query | rshash::views::kmer_view({.window_size = shape.length})) {
                    // forward: the shape on the query, reverse: the mirrored shape on its reverse complement
                    const uint64_t kmer_fwd = _pext_u64(window.value, shape.mask.lo) | (_pext_u64(window.value_hi, shape.mask.hi) << shape.lo_weight);
                    const uint64_t kmer_rev = _pext_u64(window.value_rev, shape.mask_rev.lo) | (_pext_u64(window.value_rev_hi, shape.mask_rev.hi) << shape.rev_lo_weight);
                    for(const uint64_t n : {occurrences(hts[2*i], kmer_fwd), occurrences(hts[2*i + 1], kmer_rev)}) {
                        found_positions += n;
                        found_kmers += n > 0;
                        found |= n > 0;
                        kmers++;
                    }
                }
            }
            found_queries += found;
        }
    }
    else {
        for (auto& query : queries) {
            bool query_found = false;
            for (auto&& window : query | rshash::views::kmer_view({.window_size = args.k})) {
                uint64_t kmer_value = std::min<uint64_t>(window.value, window.value_rev);
                const uint64_t n = occurrences(hts[0], kmer_value);
                query_found |= n > 0;
                found_kmers += n > 0;
                found_positions += n;
                kmers++;
            }
            found_queries += query_found;
        }
    }

    std::cout << "==== query report:\n";
    std::cout << "num_reads = " << queries.size() << '\n';
    std::cout << "num_shape_occurrences = " << kmers << '\n';
    std::cout << "num_positive_shape_occurrences = " << found_kmers << " (" << (double) found_kmers/kmers*100 << "%)\n";
    std::cout << "num_positions = " << found_positions << '\n';
    std::cout << "num_found_queries = " << found_queries << " (" << (double) found_queries/queries.size()*100 << "%)\n";

    return 0;
}
