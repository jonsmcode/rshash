// compares the occurrences of rshash locate --positions file with the positions where the reads were inserted
// (<prefix>.positions.tsv of test/datasets). a found window (k-mer/shape window of a read) is
//   located:   one of its occurrences is at the position of the inserted read
//   unknown:   not located, but found in the last level whose positions were not reported (--no-last-level-positions)
//   elsewhere: only found at other positions
// a read counts as located if at least one of its windows is located, as unknown if none is located but one is
// unknown, and as elsewhere otherwise (if found at all).
// with rshash locate --per-shape-positions, the found shapes (one per window, shape and orientation) are counted the
// same way: a shape is located if one of its positions is where the read was inserted, in the read's orientation
// (forward shapes of forward reads, reverse shapes of reverse complement reads; any orientation for k-mers).
//
//   eval_locate -t reads.positions.tsv -p occurrences.txt [-l 100] [-r per_read.tsv]

#include <charconv>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sharg/all.hpp>
#include <string>
#include <vector>

struct cmd_arguments {
    std::filesystem::path t{};
    std::filesystem::path p{};
    std::filesystem::path r{};
    uint64_t l{100};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.t, sharg::config{.short_id = 't', .long_id = "truth", .description = "insertion positions (<prefix>.positions.tsv)", .required = true});
    parser.add_option(args.p, sharg::config{.short_id = 'p', .long_id = "positions", .description = "occurrences of rshash locate --positions file (- for stdin)", .required = true});
    parser.add_option(args.l, sharg::config{.short_id = 'l', .long_id = "length", .description = "read length"});
    parser.add_option(args.r, sharg::config{.short_id = 'r', .long_id = "per-read", .description = "write the counts per read to this file"});
}

struct Insertion {
    uint64_t sequence;
    uint64_t position;
    bool fwd;
};

static int64_t to_int(std::string_view s) {
    int64_t v = 0;
    std::from_chars(s.data(), s.data() + s.size(), v);
    return v;
}

int main(int argc, char** argv)
{
    sharg::parser parser{"eval_locate", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
    }
    catch (sharg::parser_error const &ext) {
        std::cerr << "[error] " << ext.what() << '\n';
        return -1;
    }

    std::vector<Insertion> truth;
    {
        std::ifstream in(args.t);
        std::string line;
        std::getline(in, line); // header
        while(std::getline(in, line)) {
            std::vector<std::string_view> f;
            std::string_view v(line);
            for(size_t s = 0, e; s <= v.size(); s = e + 1) {
                e = v.find('\t', s);
                if(e == std::string_view::npos) e = v.size();
                f.push_back(v.substr(s, e - s));
            }
            if(f.size() < 4) continue;
            if(uint64_t(to_int(f[0])) != truth.size()) {
                std::cerr << "[error] reads in " << args.t << " are not numbered 0, 1, ...\n";
                return -1;
            }
            truth.push_back({uint64_t(to_int(f[1])), uint64_t(to_int(f[2])), f[3] == "+"});
        }
    }

    std::ifstream file;
    std::istream *in = &std::cin;
    if(args.p != "-") {
        file.open(args.p);
        if(!file) {
            std::cerr << "[error] cannot read " << args.p << '\n';
            return -1;
        }
        in = &file;
    }

    std::string line;
    int64_t window_length = -1, kernel_offset = -1;
    bool per_shape = false;
    int64_t no_shapes = 0;
    std::vector<int64_t> shape_lengths;
    std::getline(*in, line);
    {
        if(line.find("per_shape=1") != std::string::npos) {
            per_shape = true;
            const size_t n = line.find("shapes="), l = line.find("shape_lengths=");
            no_shapes = to_int(std::string_view(line).substr(n + 7));
            std::string_view v = std::string_view(line).substr(l + 14);
            v = v.substr(0, v.find('\t'));
            for(size_t s = 0, e; s < v.size(); s = e + 1) {
                e = v.find(',', s);
                if(e == std::string_view::npos) e = v.size();
                shape_lengths.push_back(to_int(v.substr(s, e - s)));
            }
        }
        const size_t w = line.find("window_length="), k = line.find("kernel_offset=");
        if(line.empty() || line[0] != '#' || w == std::string::npos || k == std::string::npos) {
            std::cerr << "[error] missing header line #window_length=..\\tkernel_offset=..\n";
            return -1;
        }
        window_length = to_int(std::string_view(line).substr(w + 14));
        kernel_offset = to_int(std::string_view(line).substr(k + 14));
    }
    const int64_t L = args.l;

    std::ofstream per_read;
    if(!args.r.empty()) {
        per_read.open(args.r);
        per_read << "read\tstrand\tfound_windows\tlocated_windows\tunknown_windows\telsewhere_windows\toccurrences\n";
    }

    uint64_t reads_found = 0, reads_located = 0, reads_unknown = 0, reads_elsewhere = 0;
    uint64_t found_windows = 0, located_windows = 0, unknown_windows = 0, elsewhere_windows = 0;
    uint64_t occurrences = 0, located_occurrences = 0, unknown_occurrences = 0;
    // per-shape positions: per shape i (both orientations)
    const size_t n_shapes = std::max<size_t>(shape_lengths.size(), 1);
    std::vector<uint64_t> found_shapes(n_shapes, 0), located_shapes(n_shapes, 0), unknown_shapes(n_shapes, 0), elsewhere_shapes(n_shapes, 0);
    std::vector<bool> seen(truth.size(), false);

    while(std::getline(*in, line)) {
        if(line.empty())
            continue;
        const size_t tab = line.find('\t');
        const uint64_t r = to_int(std::string_view(line).substr(0, tab));
        if(r >= truth.size()) {
            std::cerr << "[error] read " << r << " is not in " << args.t << '\n';
            return -1;
        }
        seen[r] = true;
        const Insertion &ins = truth[r];

        // the occurrences of a read are sorted by window (and shape): count per window and per found shape
        uint64_t fw = 0, lw = 0, uw = 0, ew = 0, occ = 0;
        int64_t current = INT64_MIN, current_shape = -1;
        bool located = false, unknown = false, shape_located = false, shape_unknown = false;
        auto close_shape = [&]() {
            if(current_shape < 0)
                return;
            const size_t i = no_shapes > 0 ? current_shape / 2 : 0;
            found_shapes[i]++;
            if(shape_located) located_shapes[i]++;
            else if(shape_unknown) unknown_shapes[i]++;
            else elsewhere_shapes[i]++;
        };
        auto close_window = [&]() {
            if(current == INT64_MIN)
                return;
            fw++;
            if(located) lw++;
            else if(unknown) uw++;
            else ew++;
        };
        std::string_view rest = std::string_view(line).substr(tab + 1);
        for(size_t s = 0, e; s < rest.size(); s = e + 1) {
            e = rest.find(' ', s);
            if(e == std::string_view::npos) e = rest.size();
            std::string_view o = rest.substr(s, e - s);
            // per shape: <window>:<shape>:<shape start in read>:..., else <kernel offset in read>:...
            size_t c = o.find(':');
            const int64_t w = to_int(o.substr(0, c));
            int64_t shape = 0, q = w;
            if(per_shape) {
                o = o.substr(c + 1); c = o.find(':');
                shape = to_int(o.substr(0, c));
                o = o.substr(c + 1); c = o.find(':');
                q = to_int(o.substr(0, c));
            }
            if(w != current) {
                close_shape();
                close_window();
                current = w;
                current_shape = shape;
                located = unknown = shape_located = shape_unknown = false;
            }
            else if(per_shape && shape != current_shape) {
                close_shape();
                current_shape = shape;
                shape_located = shape_unknown = false;
            }
            occ++;
            if(o.substr(c + 1) == "*") {
                unknown = shape_unknown = true;
                unknown_occurrences++;
                continue;
            }
            const size_t c2 = o.find(':', c + 1);
            const uint64_t sequence = to_int(o.substr(c + 1, c2 - c - 1));
            const int64_t position = to_int(o.substr(c2 + 1));
            int64_t expected;
            bool orientation = true;
            if(per_shape) {
                // start of the shape in the text: forward read at ins.position + q, reverse complement: mirrored
                const int64_t length = shape_lengths[no_shapes > 0 ? shape / 2 : 0];
                expected = ins.fwd ? int64_t(ins.position) + q : int64_t(ins.position) + L - q - length;
                orientation = no_shapes == 0 || (shape % 2 == 0) == ins.fwd;
            }
            else {
                // kernel of the window in the text: forward read at ins.position + q, reverse complement: the window
                // (starting q - kernel_offset in the read) is mirrored within the read
                expected = ins.fwd ? int64_t(ins.position) + q
                    : int64_t(ins.position) + L - (q - kernel_offset) - window_length + kernel_offset;
            }
            if(orientation && sequence == ins.sequence && position == expected) {
                located = shape_located = true;
                located_occurrences++;
            }
        }
        if(per_shape)
            close_shape();
        close_window();

        found_windows += fw; located_windows += lw; unknown_windows += uw; elsewhere_windows += ew;
        occurrences += occ;
        reads_found += fw > 0;
        if(lw) reads_located++;
        else if(uw) reads_unknown++;
        else if(fw) reads_elsewhere++;
        if(per_read.is_open())
            per_read << r << '\t' << (ins.fwd ? '+' : '-') << '\t' << fw << '\t' << lw << '\t' << uw << '\t' << ew << '\t' << occ << '\n';
    }
    if(per_read.is_open())
        for(size_t r = 0; r < truth.size(); r++)
            if(!seen[r])
                per_read << r << '\t' << (truth[r].fwd ? '+' : '-') << "\t0\t0\t0\t0\t0\n";

    std::cout << "reads = " << truth.size() << '\n';
    std::cout << "reads_found = " << reads_found << '\n';
    std::cout << "reads_located = " << reads_located << '\n';
    std::cout << "reads_unknown = " << reads_unknown << '\n';
    std::cout << "reads_elsewhere = " << reads_elsewhere << '\n';
    std::cout << "found_windows = " << found_windows << '\n';
    std::cout << "located_windows = " << located_windows << '\n';
    std::cout << "unknown_windows = " << unknown_windows << '\n';
    std::cout << "elsewhere_windows = " << elsewhere_windows << '\n';
    std::cout << "occurrences = " << occurrences << '\n';
    std::cout << "located_occurrences = " << located_occurrences << '\n';
    std::cout << "unknown_occurrences = " << unknown_occurrences << '\n';
    if(per_shape) {
        auto sum = [](const std::vector<uint64_t> &v) { return std::accumulate(v.begin(), v.end(), uint64_t{0}); };
        std::cout << "found_shapes = " << sum(found_shapes) << '\n';
        std::cout << "located_shapes = " << sum(located_shapes) << '\n';
        std::cout << "unknown_shapes = " << sum(unknown_shapes) << '\n';
        std::cout << "elsewhere_shapes = " << sum(elsewhere_shapes) << '\n';
        for(size_t i = 0; i < n_shapes; i++) {
            std::cout << "found_shapes_" << i << " = " << found_shapes[i] << '\n';
            std::cout << "located_shapes_" << i << " = " << located_shapes[i] << '\n';
        }
    }
    return 0;
}
