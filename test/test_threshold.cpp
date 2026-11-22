// Correctness test of the k-mer/shape frequency threshold (rshash build -t).
//
// Builds rshash indexes in-process for a grid of configurations (k-mers and shapes, 1-3 levels, last level with and
// without hash tables, with and without locate, several thresholds) on synthetic, highly repetitive data and compares
// every query read and every query window against a reference that counts all k-mers/shapes of the whole text:
//   - streaming lookup:  found k-mers/shapes per read (and first_hit_lookup: whether a read has a hit), t = 0 only
//   - streaming locate:  found k-mers/shapes and number of positions per read
//   - lookup/locate of single windows (rshash bench interface; lookup t = 0 only)
// An index with a threshold only supports locate (build -t needs --loc), so t > 0 is only tested with locate.
// The reference is the one of test/hashtable.cpp (-t): a window is found by shape i if its forward key is in the
// forward structure or its reverse key in the reverse structure (only the placements within the query), and only if
// the frequency of the shape value (forward + reverse occurrences in the text) is less than t. For k-mers, the
// frequency is the number of text positions of the canonical k-mer. Every index is saved and loaded before querying.
//
//   ./test/test_threshold [--seed N] [--quick] [--input file.fa --query file.fq]

#include <sharg/all.hpp>
#include <filesystem>
#include <random>
#include "../source/rshash.hpp"

using sequence_t = seqan3::bitpacked_sequence<seqan3::dna4>;

struct my_traits : seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

static std::vector<sequence_t> load_file(const std::filesystem::path &filepath)
{
    std::vector<sequence_t> output;
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto &record : stream) {
        sequence_t seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
    return output;
}


// ------------------------------------------------------------------------------------------------ synthetic data

struct Generator {
    std::mt19937_64 rng;
    explicit Generator(uint64_t seed) : rng(seed) {}

    uint8_t base() { return std::uniform_int_distribution<int>(0, 3)(rng); }
    size_t uniform(size_t a, size_t b) { return std::uniform_int_distribution<size_t>(a, b)(rng); }
    bool coin(double p) { return std::bernoulli_distribution(p)(rng); }

    std::vector<uint8_t> random(size_t n) {
        std::vector<uint8_t> s(n);
        for (auto &b : s) b = base();
        return s;
    }
    static std::vector<uint8_t> rc(std::vector<uint8_t> s) {
        std::reverse(s.begin(), s.end());
        for (auto &b : s) b = 3 - b;
        return s;
    }
    std::vector<uint8_t> mutate(std::vector<uint8_t> s, double rate) {
        for (auto &b : s)
            if (coin(rate)) b = (b + uniform(1, 3)) & 3;
        return s;
    }
};

static sequence_t to_sequence(const std::vector<uint8_t> &s)
{
    sequence_t seq;
    for (uint8_t b : s) seq.push_back(seqan3::dna4{}.assign_rank(b));
    return seq;
}

// sequences with repeats of many frequencies: repeat families (exact, mutated, reverse complemented copies), tandem
// repeats, homopolymer runs, reverse-complement palindromes, and short sequences
static std::vector<sequence_t> make_text(Generator &g)
{
    std::vector<std::vector<uint8_t>> families;
    for (size_t i = 0; i < 12; i++)
        families.push_back(g.random(g.uniform(40, 400)));

    std::vector<sequence_t> text;
    for (size_t s = 0; s < 40; s++) {
        std::vector<uint8_t> seq;
        const size_t target = (s % 10 == 9) ? g.uniform(5, 60) : g.uniform(500, 12000);
        while (seq.size() < target) {
            std::vector<uint8_t> piece;
            switch (g.uniform(0, 9)) {
            case 0: case 1: case 2:
                piece = g.random(g.uniform(1, 300));
                break;
            case 3: case 4: {
                // copy of a repeat family; the first families occur very often
                const size_t f = std::min(g.uniform(0, families.size() - 1), g.uniform(0, families.size() - 1));
                piece = families[f];
                if (g.coin(0.5)) piece = g.mutate(piece, g.coin(0.5) ? 0.0 : 0.02);
                if (g.coin(0.4)) piece = Generator::rc(piece);
                if (g.coin(0.3)) piece = std::vector<uint8_t>(piece.begin() + g.uniform(0, piece.size() / 2), piece.end());
                break;
            }
            case 5: {
                // tandem repeat
                const auto unit = g.random(g.uniform(1, 12));
                const size_t copies = g.uniform(3, 60);
                for (size_t c = 0; c < copies; c++) piece.insert(piece.end(), unit.begin(), unit.end());
                if (g.coin(0.3)) piece = g.mutate(piece, 0.01);
                break;
            }
            case 6:
                // homopolymer run
                piece = std::vector<uint8_t>(g.uniform(5, 120), g.base());
                break;
            case 7: {
                // reverse-complement palindrome rc(Y) u Y with a palindromic core u
                const auto half = g.random(g.uniform(2, 8));
                auto core = half;
                const auto half_rc = Generator::rc(half);
                core.insert(core.end(), half_rc.begin(), half_rc.end());
                const auto y = g.random(g.uniform(10, 40));
                piece = Generator::rc(y);
                piece.insert(piece.end(), core.begin(), core.end());
                piece.insert(piece.end(), y.begin(), y.end());
                break;
            }
            default: {
                // copy of an earlier part of this sequence
                if (seq.size() > 100) {
                    const size_t len = g.uniform(30, std::min<size_t>(500, seq.size()));
                    const size_t start = g.uniform(0, seq.size() - len);
                    piece.assign(seq.begin() + start, seq.begin() + start + len);
                    if (g.coin(0.5)) piece = Generator::rc(piece);
                }
                break;
            }
            }
            seq.insert(seq.end(), piece.begin(), piece.end());
        }
        text.push_back(to_sequence(seq));
    }
    return text;
}

// reads sampled from the text (both strands, with and without errors, also across sequence ends) and random reads
static std::vector<sequence_t> make_queries(Generator &g, const std::vector<sequence_t> &text, const size_t n)
{
    std::vector<uint8_t> all;
    for (auto &seq : text)
        for (auto b : seq) all.push_back(seqan3::to_rank(b));

    std::vector<sequence_t> queries;
    for (size_t i = 0; i < n; i++) {
        std::vector<uint8_t> read;
        const size_t len = g.uniform(1, 10) == 1 ? g.uniform(1, 40) : g.uniform(40, 200);
        if (g.coin(0.1) || all.size() <= len)
            read = g.random(len);
        else {
            const size_t start = g.uniform(0, all.size() - len);
            read.assign(all.begin() + start, all.begin() + start + len);
            const double r = g.uniform(0, 3) == 0 ? 0.0 : 0.02;
            read = g.mutate(read, r);
            if (g.coin(0.5)) read = Generator::rc(read);
        }
        queries.push_back(to_sequence(read));
    }
    return queries;
}


// ------------------------------------------------------------------------------------------------ reference

struct Reference {
    bool use_shapes;
    Shapes64 shapes;
    size_t no_shapes, window_size;
    ShapeBorders borders;
    uint64_t t;
    // number of text positions of every key of every structure (forward/reverse per shape, canonical k-mers)
    std::vector<gtl::flat_hash_map<uint64_t, uint64_t>> counts;

    Reference(const std::vector<sequence_t> &text, const bool use_shapes, const Shapes64 &shapes, const size_t k) :
        use_shapes(use_shapes), shapes(shapes), no_shapes(use_shapes ? shapes.shapes.size() : 0),
        window_size(use_shapes ? shapes.length : k), t(0)
    {
        if (use_shapes) borders = shape_borders(shapes);
        counts.resize(use_shapes ? 2 * no_shapes : 1);
        for (size_t j = 0; j < counts.size(); j++) {
            for (auto &seq : text) {
                if (use_shapes) {
                    const Shape64 &shape = shapes.shapes[j / 2];
                    for (auto &&w : seq | rshash::views::kmer_view({.window_size = shape.length})) {
                        if (j % 2 == 0)
                            counts[j][_pext_u64(w.value, shape.mask.lo) | (_pext_u64(w.value_hi, shape.mask.hi) << shape.lo_weight)]++;
                        else
                            counts[j][_pext_u64(w.value, shape.mask_rev.lo) | (_pext_u64(w.value_hi, shape.mask_rev.hi) << shape.rev_lo_weight)]++;
                    }
                }
                else
                    for (auto &&w : seq | rshash::views::kmer_view({.window_size = k}))
                        counts[0][std::min<uint64_t>(w.value, w.value_rev)]++;
            }
        }
    }

    uint64_t count_of(const size_t j, const uint64_t key) const {
        const auto it = counts[j].find(key);
        return it == counts[j].end() ? 0 : it->second;
    }

    static uint64_t pext128(const uint64_t lo, const uint64_t hi, const mask128_t &mask) {
        return _pext_u64(lo, mask.lo) | (_pext_u64(hi, mask.hi) << std::popcount(mask.lo));
    }

    // number of positions of the window (active placements, shapes with frequency < t). the frequency of a shape value
    // x is the number of its forward occurrences (forward structure) plus the number of its occurrences on the reverse
    // strand (reverse structure, key: reversed shape of the reverse complement). the forward placement of shape i in
    // the window and the reverse placement (in the reverse complement window) are frequency-checked separately: if
    // the run of the shape is longer than the common kernel, they are different placements with different values
    uint64_t positions(const rshash::kmerview_result &w, const uint64_t active) const {
        const uint64_t limit = t == 0 ? std::numeric_limits<uint64_t>::max() : t;
        if (use_shapes) {
            uint64_t n = 0;
            const unsigned L = window_size;
            for (size_t i = 0; i < no_shapes; i++) {
                const Shape64 &shape = shapes.shapes[i];
                // forward placement: its value, and the reverse structure key of the same value
                const uint64_t kf = pext128(w.value, w.value_hi, shape.w_mask);
                const uint64_t kf_rev = pext128(w.value_rev, w.value_rev_hi, mask128_shl(shape.mask_rev, 2 * (L - shape.w_dist_right - shape.length)));
                // reverse placement: its reverse structure key, and the value (forward structure key)
                const uint64_t kr = pext128(w.value_rev, w.value_rev_hi, shape.w_mask_rev);
                const uint64_t kr_fwd = pext128(w.value, w.value_hi, mask128_shl(shape.mask, 2 * (L - shape.w_rev_dist_right - shape.length)));
                const uint64_t fwd = count_of(2*i, kf), rev = count_of(2*i+1, kr);
                if (((active >> (2*i)) & 1) && fwd + count_of(2*i+1, kf_rev) < limit)
                    n += fwd;
                if (((active >> (2*i+1)) & 1) && count_of(2*i, kr_fwd) + rev < limit)
                    n += rev;
            }
            return n;
        }
        const uint64_t c = count_of(0, std::min<uint64_t>(w.value, w.value_rev));
        return c < limit ? c : 0;
    }

    // calls f(window, active) for the windows of a query as the streaming lookup/locate scans them
    template <typename F>
    void for_each_window(const sequence_t &query, F &&f) const {
        if (use_shapes) {
            int64_t w = -1;
            for (auto &&window : padded_query(query, borders) | rshash::views::kmer_view({.window_size = window_size})) {
                const uint64_t active = borders.active(++w, query.size());
                if (active) f(window, active);
            }
        }
        else
            for (auto &&window : query | rshash::views::kmer_view({.window_size = window_size}))
                f(window, ~uint64_t{0});
    }
};


// ------------------------------------------------------------------------------------------------ test driver

struct Config {
    std::string name;
    uint8_t k;                      // k-mers (no shapes) only
    std::vector<uint64_t> shapes;   // empty: k-mers
    uint8_t m1, m2, m3, t1, t2, t3;
};

struct Failure {
    size_t count = 0;
    void report(const std::string &msg) {
        if (count++ < 20) std::cerr << "  MISMATCH " << msg << '\n';
    }
};

int main(int argc, char **argv)
{
    sharg::parser parser{"test_threshold", argc, argv};
    uint64_t seed = 42;
    bool quick = false;
    std::filesystem::path input{}, query_file{};
    std::vector<uint64_t> thresholds_arg{};
    std::vector<std::string> only{};
    bool debug = false;
    parser.add_option(seed, sharg::config{.long_id = "seed", .description = "seed of the synthetic data"});
    parser.add_flag(quick, sharg::config{.long_id = "quick", .description = "fewer configurations"});
    parser.add_option(input, sharg::config{.short_id = 'i', .long_id = "input", .description = "text instead of synthetic data"});
    parser.add_option(query_file, sharg::config{.short_id = 'q', .long_id = "query", .description = "queries instead of synthetic reads"});
    parser.add_option(thresholds_arg, sharg::config{.short_id = 't', .long_id = "thresholds", .description = "thresholds to test"});
    parser.add_option(only, sharg::config{.long_id = "only", .description = "only the configurations whose name starts with one of these prefixes (e.g. k for k-mers, shape)"});
    parser.add_flag(debug, sharg::config{.long_id = "debug", .description = "k-mers: print text and minimiser positions of the first wrongly located k-mer"});
    try { parser.parse(); }
    catch (sharg::parser_error const &e) { std::cerr << "[PARSER ERROR] " << e.what() << '\n'; return -1; }

    Generator g(seed);
    const std::vector<sequence_t> text = input.empty() ? make_text(g) : load_file(input);
    const std::vector<sequence_t> queries = query_file.empty() ? make_queries(g, text, quick ? 300 : 1500) : load_file(query_file);
    size_t text_length = 0;
    for (auto &s : text) text_length += s.size();
    std::cout << "text: " << text.size() << " sequences, " << text_length << " bases; queries: " << queries.size() << "\n";

    std::vector<Config> configs = {
        {"k31",            31, {},                                          11, 13, 15, 2, 4, 8},
        {"k21",            21, {},                                          9, 11, 13, 3, 3, 6},
        {"k15-m-even",     15, {},                                          8, 10, 12, 2, 4, 4},
        {"k27-m-even",     27, {},                                          8, 10, 12, 4, 8, 16},
        {"k31-large-t",    31, {},                                          10, 12, 14, 16, 32, 64},
        {"k25-m-even-t",   25, {},                                          10, 12, 14, 8, 16, 32},
        {"shape-asym",      0, {2147483641},                                11, 13, 15, 2, 4, 8},
        {"shape-sym",       0, {1610612733},                                12, 14, 16, 3, 4, 8},
        {"shapes-2",        0, {2147483641, 1610612733},                    10, 12, 14, 2, 4, 8},
        {"shapes-3",        0, {2147483641, 1610612733, 2147483643},        10, 12, 14, 2, 4, 8},
        {"shape-long",      0, {471022436351},                              11, 13, 15, 2, 4, 8},
        // odd minimiser lengths (required with a threshold) for the configurations above with even ones
        {"k15-m-odd",      15, {},                                          7, 9, 11, 2, 4, 4},
        {"k27-m-odd",      27, {},                                          9, 11, 13, 4, 8, 16},
        {"k31-large-t-odd", 31, {},                                         11, 13, 15, 16, 32, 64},
        {"k25-m-odd-t",    25, {},                                          11, 13, 15, 8, 16, 32},
        {"shape-sym-odd",   0, {1610612733},                                13, 15, 17, 3, 4, 8},
        {"shapes-2-odd",    0, {2147483641, 1610612733},                    11, 13, 15, 2, 4, 8},
        {"shapes-3-odd",    0, {2147483641, 1610612733, 2147483643},        11, 13, 15, 2, 4, 8},
        {"shapes-2-large-t", 0, {2147483641, 1610612733},                   11, 13, 15, 16, 32, 64},
    };
    std::vector<uint64_t> thresholds = thresholds_arg.empty() ? std::vector<uint64_t>{0, 1, 2, 3, 4, 6, 9, 17, 33, 100, 1000} : thresholds_arg;
    if (quick) {
        const std::vector<std::string> quick_configs{"k31", "k15-m-odd", "shape-asym", "shapes-2-odd", "shape-long"};
        std::erase_if(configs, [&](const Config &c) { return std::ranges::find(quick_configs, c.name) == quick_configs.end(); });
        if (thresholds_arg.empty()) thresholds = {0, 2, 5, 17};
    }

    if (!only.empty())
        std::erase_if(configs, [&](const Config &c) {
            return std::ranges::none_of(only, [&](const std::string &p) { return c.name.starts_with(p); });
        });

    const std::filesystem::path index_file = std::filesystem::temp_directory_path() / ("test_threshold_" + std::to_string(::getpid()) + ".rshash");
    size_t total_failures = 0, total_runs = 0;
    std::streambuf *cout_buf = std::cout.rdbuf();
    std::ostringstream sink;

    for (const Config &config : configs) {
        const bool use_shapes = !config.shapes.empty();
        std::cout.rdbuf(sink.rdbuf());
        const Shapes64 shapes = use_shapes ? shape64_create(config.shapes) : Shapes64{};
        std::cout.rdbuf(cout_buf);
        Reference ref(text, use_shapes, shapes, config.k);

        for (const uint8_t level : {1, 2, 3})
        for (const bool ht : {true, false})
        for (const bool loc : {false, true})
        for (const uint64_t t : thresholds) {
            // the threshold is only supported for locate (rshash build -t needs --loc)
            if (t > 0 && !loc)
                continue;
            ref.t = t;
            std::ostringstream label;
            label << config.name << " level=" << int(level) << " ht=" << ht << " loc=" << loc << " t=" << t;

            // a threshold needs minimisers of odd length on all levels: the build has to reject the others
            const uint8_t ms[] = {config.m1, config.m2, config.m3};
            bool even_m = false;
            for (int i = 0; i < level; i++) even_m |= ms[i] % 2 == 0;
            const bool reject = t > 0 && even_m;

            sink.str("");
            std::cout.rdbuf(sink.rdbuf());
            bool rejected = false;
            {
                RSHash index = use_shapes
                    ? RSHash(shapes, level, config.m1, config.m2, config.m3, config.t1, config.t2, config.t3, t, loc, ht)
                    : RSHash(config.k, level, config.m1, config.m2, config.m3, config.t1, config.t2, config.t3, t, loc, ht);
                try {
                    index.build(text);
                    index.save(index_file);
                }
                catch (const std::invalid_argument &) { rejected = true; }
            }
            if (reject || rejected) {
                std::cout.rdbuf(cout_buf);
                total_runs++;
                total_failures += reject != rejected;
                std::cout << (reject == rejected ? "ok   " : "FAIL ") << label.str() << (rejected ? "  rejected (even minimiser length)" : "  not rejected") << std::endl;
                continue;
            }
            const uint64_t index_size = std::filesystem::file_size(index_file);
            RSHash index;
            index.load(index_file);
            std::cout.rdbuf(cout_buf);

            Failure fail;
            if (index.has_lookup() != (t == 0))
                fail.report(label.str() + ": has_lookup() is " + std::to_string(index.has_lookup()));
            uint64_t ref_found_total = 0, ref_positions_total = 0, windows_total = 0;
            std::vector<uint64_t> positions;

            for (size_t qi = 0; qi < queries.size(); qi++) {
                const auto &query = queries[qi];
                uint64_t ref_found = 0, ref_positions = 0, windows = 0;
                ref.for_each_window(query, [&](const auto &w, const uint64_t active) {
                    const uint64_t p = ref.positions(w, active);
                    ref_found += p > 0;
                    ref_positions += p;
                    windows++;
                });
                ref_found_total += ref_found;
                ref_positions_total += ref_positions;
                windows_total += windows;

                if (index.query_windows(query.size()) != windows)
                    fail.report(label.str() + " read " + std::to_string(qi) + ": windows " + std::to_string(index.query_windows(query.size())) + " != " + std::to_string(windows));

                if (t == 0) {
                    uint64_t found = 0, extensions = 0;
                    const bool hit = index.streaming_lookup(query, found, extensions);
                    if (found != ref_found || hit != (ref_found > 0))
                        fail.report(label.str() + " read " + std::to_string(qi) + ": streaming lookup found " + std::to_string(found) + ", reference " + std::to_string(ref_found));
                    uint64_t found_first = 0;
                    if (index.first_hit_lookup(query, found_first, extensions) != (ref_found > 0))
                        fail.report(label.str() + " read " + std::to_string(qi) + ": first hit lookup differs, reference " + std::to_string(ref_found));
                }

                if (loc) {
                    uint64_t found_kmers = 0, no_positions = 0;
                    index.streaming_locate(query, positions, found_kmers, no_positions);
                    if (found_kmers != ref_found || no_positions != ref_positions)
                        fail.report(label.str() + " read " + std::to_string(qi) + ": streaming locate found " + std::to_string(found_kmers) + "/" + std::to_string(no_positions) +
                                    " positions, reference " + std::to_string(ref_found) + "/" + std::to_string(ref_positions));
                }

                // single windows of the query (all placements, as rshash bench)
                if (query.size() >= ref.window_size) {
                    for (auto &&w : query | rshash::views::kmer_view({.window_size = ref.window_size})) {
                        const uint64_t p = ref.positions(w, ~uint64_t{0});
                        const std::vector<mask128_t> kmer{mask128_t{w.value, w.value_hi}};
                        const uint64_t f = t == 0 ? index.lookup(kmer) : (p > 0);
                        if (f != (p > 0))
                            fail.report(label.str() + " read " + std::to_string(qi) + ": single lookup " + std::to_string(f) + ", reference positions " + std::to_string(p));
                        if (loc) {
                            uint64_t no_positions = 0;
                            const uint64_t fl = index.locate(kmer, positions, no_positions);
                            if ((fl != (p > 0) || no_positions != p) && debug && use_shapes) {
                                debug = false;
                                std::cerr << "window " << std::bitset<64>(w.value_hi) << " " << std::bitset<64>(w.value) << " rshash " << no_positions << "\n";
                                for (size_t i = 0; i < ref.no_shapes; i++) {
                                    const Shape64 &shape = shapes.shapes[i];
                                    const uint64_t kf = _pext_u64(w.value, shape.w_mask.lo) | (_pext_u64(w.value_hi, shape.w_mask.hi) << shape.w_lo_weight);
                                    const uint64_t kr = _pext_u64(w.value_rev, shape.w_mask_rev.lo) | (_pext_u64(w.value_rev_hi, shape.w_mask_rev.hi) << shape.w_rev_lo_weight);
                                    std::cerr << "  shape " << i << " (" << shape.value << ", w_dist_right " << shape.w_dist_right << ", w_rev_dist_right " << shape.w_rev_dist_right
                                              << ", length " << shape.length << "): fwd " << ref.count_of(2*i, kf) << " rev " << ref.count_of(2*i+1, kr) << "\n";
                                }
                                std::cerr << "  shapes length " << shapes.length << " overlap " << shapes.overlap << " kernel " << shapes.kernel_length << "\n";
                            }
                            if ((fl != (p > 0) || no_positions != p) && debug && !use_shapes) { debug = false;
                                const uint64_t c = std::min<uint64_t>(w.value, w.value_rev);
                                std::cerr << "k-mer occurrences (text position, strand):";
                                uint64_t base = 32;
                                for (auto &seq : text) {
                                    uint64_t pos = base;
                                    for (auto &&x : seq | rshash::views::kmer_view({.window_size = config.k})) {
                                        if (std::min<uint64_t>(x.value, x.value_rev) == c)
                                            std::cerr << " (" << pos << "," << (x.value == c ? '+' : '-') << ")";
                                        pos++;
                                    }
                                    base += seq.size();
                                }
                                std::cerr << "\nlevel 1 minimisers emitted in these windows (value@position):";
                                base = 32;
                                for (auto &seq : text) {
                                    std::vector<uint64_t> occ;
                                    uint64_t pos = base;
                                    for (auto &&x : seq | rshash::views::kmer_view({.window_size = config.k})) {
                                        if (std::min<uint64_t>(x.value, x.value_rev) == c) occ.push_back(pos);
                                        pos++;
                                    }
                                    for (auto &&mm : seq | rshash::views::xor_minimiser_and_positions({.minimiser_size = config.m1, .window_size = config.k, .seed = seed1})) {
                                        const uint64_t q = base + mm.range_position;
                                        for (uint64_t o : occ)
                                            if (q >= o && q + config.m1 <= o + config.k) { std::cerr << " " << mm.minimiser_value << "@" << q; break; }
                                    }
                                    base += seq.size();
                                }
                                std::cerr << "\n";
                            }
                            if (fl != (p > 0) || no_positions != p)
                                fail.report(label.str() + " read " + std::to_string(qi) + ": single locate " + std::to_string(fl) + "/" + std::to_string(no_positions) + ", reference " + std::to_string(p));
                        }
                    }
                }
            }

            total_runs++;
            total_failures += fail.count > 0;
            std::cout << (fail.count ? "FAIL " : "ok   ") << label.str() << "  windows=" << windows_total << " found=" << ref_found_total
                      << " positions=" << ref_positions_total << " index=" << index_size << "B" << (fail.count ? "  mismatches=" + std::to_string(fail.count) : "") << std::endl;
        }
    }

    std::filesystem::remove(index_file);
    std::cout << "\n" << (total_runs - total_failures) << "/" << total_runs << " configurations passed\n";
    return total_failures ? 1 : 0;
}
