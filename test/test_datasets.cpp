#include <filesystem>
#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

using namespace seqan3::literals;



struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    uint8_t l{100};
    uint8_t e{5};
    uint64_t n{10000};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.l, sharg::config{.short_id = 'l', .long_id = "length", .description = "sequence length"});
    parser.add_option(args.n, sharg::config{.short_id = 'n', .long_id = "number", .description = "number of sequences"});
    parser.add_option(args.e, sharg::config{.short_id = 'e', .long_id = "error", .description = "error rate"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw sharg::user_input_error("provide input file.");

    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


size_t load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &output)
{
    size_t length = 0;
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
        length += record.sequence().size();
    }
    return length;
}

void print_file(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> & text,
    const std::filesystem::path &filepath)
{
    seqan3::sequence_file_output fout{filepath};
    size_t i = 0;
    for (auto & record : text) {
        std::string id = "seq=" + std::to_string(i++);
        fout.emplace_back(record, id);
    }
}

void random_dnas(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const size_t n, const uint8_t length)
{
    for(size_t j=0; j < n; j++) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        for(uint8_t i=0; i < length; i++) {
            seqan3::dna4 d;
            char r = rand() % 4;
            seq.push_back(d.assign_rank(r));
        }
        sequences.push_back(std::move(seq));
    }
}

void insert_sequence(seqan3::bitpacked_sequence<seqan3::dna4> &text,
    const seqan3::bitpacked_sequence<seqan3::dna4> &seq, const uint64_t position)
{
    for(uint64_t i=0; i < seq.size(); i++)
        text[position+i] = seq[i];
}

void insert_reads(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &text,
    const size_t text_length, const size_t read_length,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &reads)
{
    const uint64_t distance = text_length / reads.size(); // assert >= 1, ignore remainder
    uint64_t seq_no = 0;
    uint64_t position = 0;
    for(uint64_t i=0; i < reads.size(); i++) {
        if(text[seq_no].size() < read_length)
            throw std::runtime_error("text sequence shorter than read length");
        if(position + read_length > text[seq_no].size()) {
            seq_no++;
            position = 0;
        }
        insert_sequence(text[seq_no], reads[i], position);
        position += distance;
    }
}

void insert_errors(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const uint8_t error_rate)
{
    for (auto &seq : sequences) {
        for (size_t i = 0; i < seq.size(); i++) {
            if (rand() % 100 < error_rate) {
                seqan3::dna4 value = seq[i];
                value.assign_rank(value.to_rank() ^ 0b11);
                seq[i] = value;
            }
        }
    }
}

int main(int argc, char** argv)
{
    sharg::parser parser{"datasets", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
        check_arguments(parser, args);
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::cout << "building random sequences...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> reads;
    random_dnas(reads, args.n, args.l);

    std::cout << "loading text...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
    const size_t text_length = load_file(args.i, text);

    std::cout << "inserting random sequences into text...\n";
    insert_reads(text, text_length, args.l, reads);
    print_file(text, args.i.string() + ".test.fasta");
    
    std::cout << "writing errors into random sequences...\n";
    insert_errors(reads, args.e);
    print_file(reads, args.i.string() + ".reads.fasta");

    return 0;
}

