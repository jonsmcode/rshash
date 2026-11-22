// sensitivity datasets from real reads: n reads of length l are sampled from each genome (uniformly over all positions
// where a read fits into one sequence, half of them as reverse complement), and copies with e substitution errors per
// read are written, one file per error count. the foreign reads of a genome mix the reads of all other genomes in equal
// parts (same errors), so querying them against the genome's index measures hits by chance between genomes. random
// reads (not from any genome) are written as a second control.
//
// output in directory -o, <name> is the input file name without .fa/.fasta(.gz):
//   <name>.e<e>.reads.fasta     reads sampled from the genome with e errors per read, for every -e
//   <name>.e<e>.foreign.fasta   reads sampled from the other genomes with e errors per read (needs >= 2 genomes)
//   <name>.positions.tsv        where each read was sampled: read index (as in the read files), sequence index,
//                               position in the sequence, strand (+ the read, - its reverse complement)
//   random.control.fasta        random reads

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>
#include <regex>
#include <sharg/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/output.hpp>

using namespace seqan3::literals;


struct cmd_arguments {
    std::vector<std::filesystem::path> i{};
    std::filesystem::path o{};
    uint64_t l{100};
    std::vector<uint64_t> e{};
    uint64_t n{1000000};
    uint64_t seed{42};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "genome file, multiple can be given"});
    parser.add_option(args.o, sharg::config{.short_id = 'o', .long_id = "output", .description = "output directory"});
    parser.add_option(args.l, sharg::config{.short_id = 'l', .long_id = "length", .description = "read length"});
    parser.add_option(args.n, sharg::config{.short_id = 'n', .long_id = "number", .description = "number of reads per genome"});
    parser.add_option(args.e, sharg::config{.short_id = 'e', .long_id = "errors", .description = "substitution errors per read, multiple values can be given (default 0..10)"});
    parser.add_option(args.seed, sharg::config{.short_id = 's', .long_id = "seed", .description = "random seed"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(args.i.empty())
        throw sharg::user_input_error("provide input files.");
    if(!parser.is_option_set('o'))
        throw sharg::user_input_error("provide output directory.");
    if(args.e.empty())
        for(uint64_t e = 0; e <= 10; e++)
            args.e.push_back(e);
    for(uint64_t e : args.e)
        if(e > args.l)
            throw sharg::user_input_error("more errors than read length.");
    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

using sequence_t = seqan3::bitpacked_sequence<seqan3::dna4>;


std::string genome_name(const std::filesystem::path &filepath)
{
    return std::regex_replace(filepath.filename().string(), std::regex("\\.(fa|fasta)(\\.gz)?$"), "");
}

void load_file(const std::filesystem::path &filepath, std::vector<sequence_t> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        sequence_t seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
}

void print_file(const std::vector<sequence_t> &sequences, const std::filesystem::path &filepath)
{
    seqan3::sequence_file_output fout{filepath};
    size_t i = 0;
    for (auto & record : sequences) {
        std::string id = "seq=" + std::to_string(i++);
        fout.emplace_back(record, id);
    }
}

std::vector<sequence_t> random_dnas(const size_t n, const uint64_t length, std::mt19937_64 &gen)
{
    std::vector<sequence_t> sequences;
    sequences.reserve(n);
    for(size_t j = 0; j < n; j++) {
        sequence_t seq;
        for(uint64_t i = 0; i < length; i++)
            seq.push_back(seqan3::dna4{}.assign_rank(gen() % 4));
        sequences.push_back(std::move(seq));
    }
    return sequences;
}

// n reads of length l at uniform positions of the text (each position where a read fits into one sequence is equally
// likely), every second one as reverse complement
std::vector<sequence_t> sample_reads(const std::vector<sequence_t> &text, const size_t n, const uint64_t l,
    std::mt19937_64 &gen, std::ostream &sampled_at)
{
    sampled_at << "read\tsequence\tposition\tstrand\n";
    // starts[s]: number of read positions in sequences before s
    std::vector<uint64_t> starts(text.size() + 1, 0);
    for(size_t s = 0; s < text.size(); s++)
        starts[s + 1] = starts[s] + (text[s].size() >= l ? text[s].size() - l + 1 : 0);
    const uint64_t positions = starts.back();
    if(positions == 0)
        throw std::runtime_error("no sequence is long enough for a read");

    std::vector<sequence_t> reads;
    reads.reserve(n);
    for(size_t j = 0; j < n; j++) {
        const uint64_t p = gen() % positions;
        const size_t s = std::upper_bound(starts.begin(), starts.end(), p) - starts.begin() - 1;
        const uint64_t pos = p - starts[s];
        sequence_t read;
        if(j % 2 == 0)
            for(uint64_t i = 0; i < l; i++)
                read.push_back(text[s][pos + i]);
        else
            for(uint64_t i = 0; i < l; i++)
                read.push_back(seqan3::dna4{}.assign_rank(text[s][pos + l - 1 - i].to_rank() ^ 0b11));
        reads.push_back(std::move(read));
        sampled_at << j << '\t' << s << '\t' << pos << '\t' << (j % 2 == 0 ? '+' : '-') << '\n';
    }
    return reads;
}

// exactly `errors` substitutions per read at distinct positions, each to one of the three other bases
std::vector<sequence_t> with_errors(std::vector<sequence_t> sequences, const uint64_t errors, std::mt19937_64 &gen)
{
    for (auto &seq : sequences) {
        std::vector<size_t> positions(seq.size());
        for (size_t i = 0; i < positions.size(); i++)
            positions[i] = i;
        // partial Fisher-Yates: the first `errors` positions are a uniform sample
        for (size_t i = 0; i < errors; i++)
            std::swap(positions[i], positions[i + gen() % (positions.size() - i)]);
        for (size_t i = 0; i < errors; i++) {
            const size_t pos = positions[i];
            seq[pos] = seqan3::dna4{}.assign_rank((seq[pos].to_rank() + 1 + gen() % 3) % 4);
        }
    }
    return sequences;
}

int main(int argc, char** argv)
{
    sharg::parser parser{"samplereads", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
        check_arguments(parser, args);
    }
    catch (sharg::parser_error const &ext) {
        std::cerr << "[error] " << ext.what() << '\n';
        return -1;
    }

    std::mt19937_64 gen(args.seed);
    std::filesystem::create_directories(args.o);
    const size_t g = args.i.size();

    // reads[i][x]: reads of genome i with the x-th error count
    std::vector<std::vector<std::vector<sequence_t>>> reads(g);
    for(size_t i = 0; i < g; i++) {
        const std::string name = genome_name(args.i[i]);
        std::cout << name << ": loading text...\n";
        std::vector<sequence_t> text;
        load_file(args.i[i], text);
        std::cout << name << ": sampling " << args.n << " reads...\n";
        std::ofstream positions(args.o / (name + ".positions.tsv"));
        const std::vector<sequence_t> sampled = sample_reads(text, args.n, args.l, gen, positions);
        for(uint64_t e : args.e) {
            reads[i].push_back(with_errors(sampled, e, gen));
            print_file(reads[i].back(), args.o / (name + ".e" + std::to_string(e) + ".reads.fasta"));
        }
    }

    // foreign reads of genome i: n / (g - 1) reads of every other genome, the remainder from the first ones
    if(g >= 2)
        for(size_t i = 0; i < g; i++) {
            const std::string name = genome_name(args.i[i]);
            std::cout << name << ": writing foreign reads...\n";
            for(size_t x = 0; x < args.e.size(); x++) {
                std::vector<sequence_t> foreign;
                foreign.reserve(args.n);
                size_t others = 0;
                for(size_t j = 0; j < g; j++) {
                    if(j == i)
                        continue;
                    const size_t count = args.n / (g - 1) + (others < args.n % (g - 1) ? 1 : 0);
                    foreign.insert(foreign.end(), reads[j][x].begin(), reads[j][x].begin() + count);
                    others++;
                }
                print_file(foreign, args.o / (name + ".e" + std::to_string(args.e[x]) + ".foreign.fasta"));
            }
        }

    std::cout << "writing random control reads...\n";
    print_file(random_dnas(args.n, args.l, gen), args.o / "random.control.fasta");
    return 0;
}
