#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "rshash.hpp"
#include "build.hpp"
#include "io.hpp"
#include "minimiser_views.hpp"


RSHash1::RSHash1() :
    endpoints(std::vector<uint64_t>{}, 1),
    r1(std::vector<uint64_t>{}, 1),
    m_hasher(seed1)
{}

RSHash1::RSHash1(
    uint8_t const k, uint8_t const m1, uint8_t const m_thres1)
    : k(k), m1(m1),
      m_thres1(m_thres1), span(k-m1+1),
      endpoints(std::vector<uint64_t>{}, 1),
      r1(std::vector<uint64_t>{}, 1),
      m_hasher(seed1),
      kmermask(compute_mask(2u * k)),
      mmermask(compute_mask(2u * m1))
{}



void RSHash1::mark_minimizer_occurences(const size_t no_skmers, const std::vector<uint8_t> &minimizer_occurences)
{
    s1 = bit_vector(no_skmers+1, 0);
    s1[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < minimizer_occurences.size(); i++) {
        j += minimizer_occurences[i];
        s1[j] = 1;
    }
    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), no_skmers+1, 3);
}


void RSHash1::fill_minimizer_offsets(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    std::vector<uint8_t> &minimizer_occurences,
    const size_t text_length, const size_t no_minimizers, const size_t no_skmers)
{
    const size_t offset_width = std::bit_width(text_length+32);
    bits::compact_vector::builder builder;
    builder.resize(no_skmers, offset_width);

    std::fill(minimizer_occurences.begin(), minimizer_occurences.end(), 0);

    size_t length = 32;
    for(auto & sequence : sequences) {
        for (auto && minimiser : sequence | rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1})) {
            if(uint64_t i = r1.rank(minimiser.minimiser_value); r1.rank(minimiser.minimiser_value+1)-i) {
                size_t s = s1_select.select(i);
                builder.set(s + minimizer_occurences[i], length + minimiser.range_position);
                minimizer_occurences[i]++;
            }
        }
        length += sequence.size();
    }
    builder.build(offsets1);
}


std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> RSHash1::get_frequent_skmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences)
{
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> freq_skmers;
    auto skmerview = rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});

    for(auto & sequence : sequences) {
        size_t start_position = 0;
        bool cur_freq, freq;

        for(auto && minimiser : sequence | skmerview) {
            freq = r1.rank(minimiser.minimiser_value+1)-r1.rank(minimiser.minimiser_value);
            break;
        }
        for(auto && minimiser : sequence | skmerview) {
            cur_freq = r1.rank(minimiser.minimiser_value+1)-r1.rank(minimiser.minimiser_value);
            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                seqan3::bitpacked_sequence<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position-1+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers.emplace_back(skmer);
            }
            freq = cur_freq;
        }
        if(!cur_freq) {
            seqan3::bitpacked_sequence<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers.emplace_back(skmer);
        }
    }

    return freq_skmers;
}


void RSHash1::build(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    no_text_kmers = mark_sequences(input, k, endpoints);
    size_t text_length = endpoints.size();

    std::vector<uint64_t> minimizers;
    std::vector<uint8_t> minimizer_occurences;
    const uint64_t no_skmers = get_unfrequent_minimizers(input, m1, m_thres1, k, seed1, minimizers, minimizer_occurences);
    const size_t no_minimizers = minimizers.size();

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimizers, M1);
    minimizers.clear();

    std::cout << "filling bitvector S_1...\n";
    mark_minimizer_occurences(no_skmers, minimizer_occurences);

    std::cout << "filling offsets_1...\n";
    fill_minimizer_offsets(input, minimizer_occurences, text_length, no_minimizers, no_skmers);
    minimizer_occurences.clear();

    std::cout << "get frequent skmers...\n";
    auto freq_skmers = get_frequent_skmers(input);

    std::cout << "build level 2...\n";
    for(auto & sequence : freq_skmers)
        for(auto && window : sequence | rshash::views::kmerview({.window_size = k}))
            hashmap.insert(std::min<uint64_t>(window.kmer_value, window.kmer_value_rev));

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    print_info();
}


const inline uint64_t RSHash1::get_word64(uint64_t pos) {
    uint64_t block = pos >> 5;
    uint64_t shift = (pos & 31) << 1;
    uint64_t lo = text[block];
    uint64_t hi = text[block + 1];

    uint64_t shift_mask = -(shift != 0);
    return (lo >> shift) | ((hi << (64 - shift)) & shift_mask);
}

const inline uint64_t RSHash1::get_base(uint64_t pos) {
    return (text[pos >> 5] >> ((pos & 31) << 1)) & 3ULL;
}


uint64_t RSHash1::access(const uint64_t unitig_id, const size_t offset) {
    return get_word64(offset) & kmermask;
}


std::vector<uint64_t> RSHash1::rand_text_kmers(const uint64_t n) {
    std::uniform_int_distribution<uint32_t> distr;
    std::mt19937 m_rand(1);
    std::vector<std::uint64_t> kmers;
    kmers.reserve(n);
    const uint64_t l = (text.size()-1)*32;

    const uint64_t no_unitigs = number_unitigs();
    for (uint64_t i = 0; i < n;) {
        const uint64_t offset = distr(m_rand) % l;

        const uint64_t r = endpoints.rank(offset+1);
        const uint64_t next_endpoint = endpoints.select(r);

        if(offset + 64 >= next_endpoint)
            continue;

        const uint64_t kmer = access(0, offset);

        if ((i & 1) == 0)
            kmers.push_back(crc(kmer, k));
        else
            kmers.push_back(kmer);

        i++;
    }

    return kmers;
}


uint64_t RSHash1::lookup(const std::vector<uint64_t> &kmers, bool verbose)
{
    uint64_t occurences = 0;

    uint64_t* offsets = new uint64_t[m_thres1-1];
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        const uint64_t kmer_rc = crc(kmer, k);
        const uint64_t minimiser = find_minimiser(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);

        if(uint64_t minimiser_rank = r1.rank(minimiser); r1.rank(minimiser + 1) - minimiser_rank) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;

            occurences += check(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else
            occurences += hashmap.contains(std::min<uint64_t>(kmer, kmer_rc));
    }

    delete[] offsets;

    return occurences;
}


inline bool RSHash1::check(const uint64_t kmer, const uint64_t kmer_rc,
    uint64_t* offsets, const size_t p, const size_t no_skmers,
    const size_t left_minimiser_position, const size_t right_minimiser_position)
{
    for(size_t i = 0; i < no_skmers; i++)
        offsets[i] = offsets1.access(p+i)-span+1;

    for(size_t i = 0; i < no_skmers; i++) {
        const uint64_t o = offsets[i];

        uint64_t hash_rc = get_word64(o + left_minimiser_position) & kmermask;
        uint64_t hash_fwd = get_word64(o + span-1-left_minimiser_position) & kmermask;

        if(kmer == hash_fwd || kmer_rc == hash_rc)
            return true;

        if(left_minimiser_position != k-m1-right_minimiser_position) {
            hash_fwd = get_word64(o + right_minimiser_position) & kmermask;
            hash_rc = get_word64(o + span-1-right_minimiser_position) & kmermask;

            if(kmer == hash_fwd || kmer_rc == hash_rc)
                return true;
        }

    }

    return false;
}


inline void RSHash1::fill_buffer(uint64_t *offsets, uint64_t *buffer, size_t p, size_t N, const uint64_t shift)
{
    for(uint64_t i = 0; i < N; i++)
        offsets[i] = offsets1.access(p+i)+1-span;

    for(uint64_t i = 0; i < N; i++) {
        const uint64_t o = offsets[i];
        
        uint64_t kmer = get_word64(o) & kmermask;
        uint64_t bits = get_word64(o + k);
        *buffer++ = kmer;
        for(size_t j=0; j < span-1; j++) {
            uint64_t const next_base = bits & 3ULL;
            bits >>= 2;
            kmer = (kmer >> 2) | (next_base << shift);
            *buffer++ = kmer;
        }
    }
}

inline bool RSHash1::check_overlap(uint64_t skmer_pos, uint64_t text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    const uint64_t r = endpoints.rank(skmer_pos+span);
    start_pos = endpoints.select(r-1, &end_pos);

    return text_pos >= start_pos && text_pos+k-1 < end_pos;
}


inline bool RSHash1::check_minimiser_pos(uint64_t *buffer, const uint64_t offset,
    const uint64_t kmer, const uint64_t kmerrc, const size_t s,
    const size_t minimiser_pos, bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    if(buffer[s+minimiser_pos] == kmerrc) {
        forward = false;
        text_pos = offset + minimiser_pos;
        if(check_overlap(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[s+span-1-minimiser_pos] == kmer) {
        forward = true;
        text_pos = offset + span-1-minimiser_pos + k - 1;
        if(check_overlap(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }

    return false;
}

inline bool RSHash1::check_minimiser_pos2(uint64_t *buffer, const uint64_t offset,
    const uint64_t kmer, const uint64_t kmerrc, const size_t s,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    if(buffer[s+left_minimiser_pos] == kmerrc) {
        forward = false;
        text_pos = offset + left_minimiser_pos;
        if(check_overlap(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[s+span-1-left_minimiser_pos] == kmer) {
        forward = true;
        text_pos = offset + span-1-left_minimiser_pos + k - 1;
        if(check_overlap(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[s+right_minimiser_pos] == kmer) {
        forward = true;
        text_pos = offset + right_minimiser_pos + k - 1;
        if(check_overlap(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[s+span-1-right_minimiser_pos] == kmerrc) {
        forward = false;
        text_pos = offset + span-1-right_minimiser_pos;
        if(check_overlap(offset, text_pos, start_pos, end_pos))
            return true;
    }

    return false;
}


inline bool RSHash1::lookup_buffer(uint64_t *buffer, uint64_t *offsets, const size_t no_minimiser,
    const uint64_t query, const uint64_t queryrc,
    uint64_t &text_pos, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &start_pos, uint64_t &end_pos)
{
    size_t s = 0;
    if(left_minimiser_pos != k-m1-right_minimiser_pos) {
        for(size_t i = 0; i < no_minimiser; i++) {
            if(check_minimiser_pos2(buffer, offsets[i], query, queryrc, s, left_minimiser_pos, right_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s += span;
        }
    }
    else {
        for(size_t i = 0; i < no_minimiser; i++) {
            if(check_minimiser_pos(buffer, offsets[i], query, queryrc, s, left_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s += span;
        }
    }
    
    return false;
}


inline bool RSHash1::extend_in_text(uint64_t &text_pos, uint64_t start, uint64_t end,
    bool forward, const uint64_t query, const uint64_t query_rc, const uint64_t shift)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = get_base(text_pos);
            return new_rank == (query >> shift);
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = get_base(text_pos);
            return new_rank == (query_rc & 0b11);
        }
    }
    return false;
}


inline uint64_t RSHash1::find_minimiser(const uint64_t kmer, const uint64_t kmer_rc, size_t &left_minimiser_position, size_t &right_minimiser_position)
{
    uint64_t mmer = kmer >> 2*(k - m1);
    uint64_t mmer_rc = kmer_rc & mmermask;
    uint64_t minimiser_value = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
    left_minimiser_position = k-m1;
    right_minimiser_position = 0;

    for (size_t i = 1; i <= k-m1; ++i) {
        mmer = (kmer >> 2*(k-m1-i)) & mmermask;
        mmer_rc = (kmer_rc >> 2*i) & mmermask;
        const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
        if(mmerhash < minimiser_value) {
            minimiser_value = mmerhash;
            left_minimiser_position = k-m1-i;
            right_minimiser_position = i;
        }
        else if(mmerhash == minimiser_value)
            left_minimiser_position = k-m1-i;
    }

    return minimiser_value;
}

inline void RSHash1::update_minimiser(const uint64_t kmer, const uint64_t kmer_rc, uint64_t &minimiser, size_t &left_minimiser_position, size_t &right_minimiser_position)
{
    if(left_minimiser_position-- == 0) {
        minimiser = find_minimiser(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
        return;
    }
    const uint64_t mmer = kmer >> 2*(k - m1);
    const uint64_t mmer_rc = kmer_rc & mmermask;
    const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);

    if(mmerhash < minimiser) {
        minimiser = mmerhash;
        left_minimiser_position = k - m1;
        right_minimiser_position = 0;
        return;
    }
    if(mmerhash == minimiser) {
        right_minimiser_position = 0;
        return;
    }

    right_minimiser_position++;
}


uint64_t RSHash1::streaming_query(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = rshash::views::kmerview({.window_size = k});

    uint64_t occurences = 0;
    uint64_t current_pos_minimiser=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser=std::numeric_limits<uint64_t>::max();
    const uint64_t shift = 2*(k-1);
    uint64_t* offsets = new uint64_t[m_thres1];
    uint64_t* kmer_buffer = new uint64_t[(m_thres1) * span];
    size_t no_minimiser;
    uint64_t text_pos, sequence_begin, sequence_end;
    bool forward;
    bool found = false;
    bool rolling = false;
    uint64_t minimiser, minimiser_rank;
    size_t left_minimiser_position, right_minimiser_position;

    for(auto && kmer : query | view)
    {
        if(found && extend_in_text(text_pos, sequence_begin, sequence_end, forward, kmer.kmer_value, kmer.kmer_value_rev, shift)) {
            occurences++;
            extensions++;
            rolling = false;
        }
        else {
            if(rolling)
                update_minimiser(kmer.kmer_value, kmer.kmer_value_rev, minimiser, left_minimiser_position, right_minimiser_position);
            else {
                minimiser = find_minimiser(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position, right_minimiser_position);
                rolling = true;
            }

            if(minimiser == current_pos_minimiser) {
                found = lookup_buffer(kmer_buffer, offsets, no_minimiser, kmer.kmer_value, kmer.kmer_value_rev, text_pos, left_minimiser_position, right_minimiser_position, forward, sequence_begin, sequence_end);
                occurences += found;
            }
            else if(minimiser != current_neg_minimiser && (minimiser_rank = r1.rank(minimiser), r1.rank(minimiser + 1) - minimiser_rank)) {
                const size_t minimiser_position = s1_select.select(minimiser_rank);
                no_minimiser = s1_select.select(minimiser_rank+1) - minimiser_position;

                fill_buffer(offsets, kmer_buffer, minimiser_position, no_minimiser, shift);
                found = lookup_buffer(kmer_buffer, offsets, no_minimiser, kmer.kmer_value, kmer.kmer_value_rev, text_pos, left_minimiser_position, right_minimiser_position, forward, sequence_begin, sequence_end);
                occurences += found;
                current_pos_minimiser = minimiser;
            }
            else {
                occurences += hashmap.contains(std::min<uint64_t>(kmer.kmer_value, kmer.kmer_value_rev));
                found = false;
                current_neg_minimiser = minimiser;
            }
        }

    }

    delete[] kmer_buffer;
    delete[] offsets;
    
    return occurences;
}


void RSHash1::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    cereal::BinaryOutputArchive archive(out);

    archive(k, m1, m_thres1, s1, endpoints, r1, offsets1, text, no_text_kmers, hashmap);

    out.close();
}

void RSHash1::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    cereal::BinaryInputArchive archive(in);

    archive(k, m1, m_thres1, s1, endpoints, r1, offsets1, text, no_text_kmers, hashmap);

    in.close();

    span = k - m1 + 1;
    kmermask = compute_mask(2u * k);
    mmermask = compute_mask(2u * m1);
    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);

    std::cout << "loaded index...\n";
}


void RSHash1::print_info() {
    const size_t N = text.size()*32;
    const size_t offset_width = std::bit_width(N);
    const uint64_t M1 = 1ULL << (m1+m1);
    const uint64_t no_minimizers = r1.rank(M1);
    const uint64_t no_skmers = s1.size();

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "textkmers: " << no_text_kmers <<  '\n';
    
    std::cout << "no minimiser: " << no_minimizers << "\n";
    std::cout << "no minimiser occurences: " << no_skmers << "\n";
    std::cout << "avg superkmers: " << (double) no_skmers/no_minimizers <<  '\n';
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/no_text_kmers*100 << "%\n";

    std::cout << "density r1: " << (double) no_minimizers/M1*100 << "%\n";
    std::cout << "density s1: " << (double) no_minimizers/(no_skmers+1)*100 <<  "%\n";

    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/no_text_kmers << "\n";
    std::cout << "endpoints: " << (double) endpoints.bitCount()/no_text_kmers << "\n";
    std::cout << "offsets1: " << (double) no_skmers*offset_width/no_text_kmers << "\n";
    std::cout << "R_1: " << (double) r1.bitCount()/no_text_kmers << "\n";
    std::cout << "S_1: " << (double) (no_skmers+1 + s1_select.bitCount())/no_text_kmers << "\n";
    std::cout << "Hashtable: " << (double) 65*hashmap.bucket_count()/no_text_kmers << "\n";
    std::cout << "total: " << (double) (no_skmers*offset_width+2*N+r1.bitCount()+no_skmers+1+s1_select.bitCount()+endpoints.bitCount()+65*hashmap.bucket_count())/no_text_kmers << "\n";
}

