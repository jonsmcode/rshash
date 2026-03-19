#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "rshash.hpp"
#include "build.hpp"
#include "io.hpp"
#include "minimiser_views.hpp"


RSHash2::RSHash2() :
    endpoints(std::vector<uint64_t>{}, 1),
    r1(std::vector<uint64_t>{}, 1),
    r2(std::vector<uint64_t>{}, 1),
    m_hasher1(seed1), m_hasher2(seed2)
{}

RSHash2::RSHash2(
    uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m_thres1, uint8_t const m_thres2)
    : k(k), m1(m1), m_thres1(m_thres1), m2(m2),
      m_thres2(m_thres2), span1(k-m1+1), span2(k-m2+1),
      endpoints(std::vector<uint64_t>{}, 1),
      r1(std::vector<uint64_t>{}, 1),
      r2(std::vector<uint64_t>{}, 1),
      m_hasher1(seed1), m_hasher2(seed2),
      kmermask(compute_mask(2u * k)),
      mmermask1(compute_mask(2u * m1)),
      mmermask2(compute_mask(2u * m2))
{}

template<int level>
void RSHash2::mark_minimizer_occurences(const size_t no_skmers, const std::vector<uint8_t> &minimizer_occurences)
{
    auto & s = [&]() -> auto& {
    if constexpr (level == 1) return s1;
    else if constexpr (level == 2) return s2;
    }();

    s = bit_vector(no_skmers+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < minimizer_occurences.size(); i++) {
        j += minimizer_occurences[i];
        s[j] = 1;
    }

    if constexpr (level == 1)
        s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), no_skmers+1, 3);
    if constexpr (level == 2)
        s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), no_skmers+1, 3);
}


template<int level>
void RSHash2::fill_minimizer_offsets(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    std::vector<size_t> &skmer_positions, std::vector<uint8_t> &minimizer_occurences,
    const size_t text_length, const size_t no_skmers)
{
    auto view = [&]() {
    if constexpr (level == 1)
        return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
    else if constexpr (level == 2)
        return rshash::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level == 1) return r1;
    else if constexpr (level == 2) return r2;
    }();

    auto & s = [&]() -> auto& {
    if constexpr (level == 1) return s1_select;
    else if constexpr (level == 2) return s2_select;
    }();


    const size_t offset_width = std::bit_width(text_length);
    bits::compact_vector::builder builder;
    builder.resize(no_skmers, offset_width);

    std::fill(minimizer_occurences.begin(), minimizer_occurences.end(), 0);

    size_t length = 32;
    size_t skmer_idx = 0;
    for(auto & sequence : sequences) {
        for (auto && minimiser : sequence | view) {
            if(uint64_t minimizer_rank = r.rank(minimiser.minimiser_value); r.rank(minimiser.minimiser_value+1)-minimizer_rank) {
                size_t minimizer_idx = s.select(minimizer_rank);
                if constexpr (level == 1)
                    builder.set(minimizer_idx + minimizer_occurences[minimizer_rank], length + minimiser.range_position);
                if constexpr (level == 2)
                    builder.set(minimizer_idx + minimizer_occurences[minimizer_rank], skmer_positions[skmer_idx] + minimiser.range_position);
                minimizer_occurences[minimizer_rank]++;
            }
        }
        skmer_idx++;
        length += sequence.size();
    }

    if constexpr (level == 1)
        builder.build(offsets1);
    if constexpr (level == 2)
        builder.build(offsets2);
}


template<int level>
void RSHash2::get_frequent_skmers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &freq_skmers, std::vector<size_t> &skmer_positions)
{
    auto skmerview = [&]() {
    if constexpr (level == 1)
        return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed = seed1});
    else if constexpr (level == 2)
        return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m2, .window_size = k, .seed = seed2});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level == 1) return r1;
    else if constexpr (level == 2) return r2;
    }();
    
    size_t length = 32;
    for(auto & sequence : sequences) {
        size_t start_position = 0;
        bool cur_freq, freq;

        for(auto && minimiser : sequence | skmerview) {
            freq = r.rank(minimiser.minimiser_value+1)-r.rank(minimiser.minimiser_value);
            break;
        }
        for(auto && minimiser : sequence | skmerview) {
            cur_freq = r.rank(minimiser.minimiser_value+1)-r.rank(minimiser.minimiser_value);
            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                seqan3::bitpacked_sequence<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position-1+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers.emplace_back(skmer);
                skmer_positions.emplace_back(length + start_position);
            }
            freq = cur_freq;
        }
        if(!cur_freq) {
            seqan3::bitpacked_sequence<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers.emplace_back(skmer);
            skmer_positions.push_back(length + start_position);
        }
        length += sequence.size();
    }
}



void RSHash2::build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input)
{
    no_text_kmers = mark_sequences(input, k, endpoints);
    size_t text_length = endpoints.size();

    std::vector<uint64_t> minimizers1;
    std::vector<uint8_t> minimizers1_occurences;
    const uint64_t no_skmers1 = get_unfrequent_minimizers(input, m1, m_thres1, k, seed1, minimizers1, minimizers1_occurences);
    const size_t no_minimizers1 = minimizers1.size();

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimizers1, M1);
    minimizers1.clear();

    std::cout << "mark skmers1...\n";
    mark_minimizer_occurences<1>(no_skmers1, minimizers1_occurences);

    std::cout << "filling offsets_1...\n";
    std::vector<size_t> skmer_positions1;
    fill_minimizer_offsets<1>(input, skmer_positions1, minimizers1_occurences, text_length, no_skmers1);
    minimizers1_occurences.clear();

    std::cout << "get frequent skmers...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> freq_skmers1;
    get_frequent_skmers<1>(input, freq_skmers1, skmer_positions1);

    std::cout << "count minimizers2...\n";
    std::vector<uint64_t> minimizers2;
    std::vector<uint8_t> minimizers2_occurences;
    const uint64_t no_skmers2 = get_unfrequent_minimizers(freq_skmers1, m2, m_thres2, k, seed2, minimizers2, minimizers2_occurences);
    const size_t no_minimizers2 = minimizers2.size();

    std::cout << "build R_2...\n";
    std::sort(minimizers2.begin(), minimizers2.end());
    const uint64_t M2 = 1ULL << (m2+m2);
    r2 = sux::bits::EliasFano(minimizers2, M2);

    std::cout << "mark skmers1...\n";
    mark_minimizer_occurences<2>(no_skmers2, minimizers2_occurences);

    std::cout << "filling offsets_2...\n";
    fill_minimizer_offsets<2>(freq_skmers1, skmer_positions1, minimizers2_occurences, text_length, no_skmers2);
    minimizers2_occurences.clear();
    skmer_positions1.clear();

    std::cout << "get frequent skmers...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> freq_skmers2;
    std::vector<size_t> skmer_positions2;
    get_frequent_skmers<2>(freq_skmers1, freq_skmers2, skmer_positions2);

    std::cout << "build level 3, HT...\n";
    for(auto & sequence : freq_skmers2) {
        for(auto && kmer : sequence | rshash::views::kmerview({.window_size = k})) {
            hashmap.insert(std::min<uint64_t>(kmer.kmer_value, kmer.kmer_value_rev));
        }
    }

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    print_info();
}


std::vector<uint64_t> RSHash2::rand_text_kmers(const uint64_t n) {
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

const inline uint64_t RSHash2::get_word64(uint64_t pos) {
    uint64_t block = pos >> 5;
    uint64_t shift = (pos & 31) << 1;
    uint64_t lo = text[block];
    uint64_t hi = text[block + 1];

    uint64_t shift_mask = -(shift != 0);
    return (lo >> shift) | ((hi << (64 - shift)) & shift_mask);
}

const inline uint64_t RSHash2::get_base(uint64_t pos) {
    return (text[pos >> 5] >> ((pos & 31) << 1)) & 3ULL;
}


uint64_t RSHash2::access(const uint64_t unitig_id, const size_t offset) {
    return get_word64(offset) & kmermask;
}


uint64_t RSHash2::lookup(const std::vector<uint64_t> &kmers, bool verbose)
{
    uint64_t occurences = 0;

    uint64_t* offsets = new uint64_t[m_thres2-1];
    uint64_t minimiser, minimiser_rank;
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        const uint64_t kmer_rc = crc(kmer, k);
        minimiser = find_minimiser<1>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);

        if(minimiser_rank = r1.rank(minimiser); r1.rank(minimiser + 1) - minimiser_rank) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;

            occurences += check<1>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            minimiser = find_minimiser<2>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
            if(minimiser_rank = r2.rank(minimiser); r2.rank(minimiser + 1) - minimiser_rank) {
                size_t p = s2_select.select(minimiser_rank);
                size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;

                occurences += check<2>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
            }
            else
                occurences += hashmap.contains(std::min<uint64_t>(kmer, kmer_rc));
        }
    }

    delete[] offsets;

    return occurences;
}


template<int level>
inline bool RSHash2::check(const uint64_t kmer, const uint64_t kmer_rc,
    uint64_t* offsets, const size_t p, const size_t no_skmers,
    const size_t left_minimiser_position, const size_t right_minimiser_position)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    
    for(size_t i = 0; i < no_skmers; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i)-span+1;
        else
            offsets[i] = offsets2.access(p+i)-span+1;
    }

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


inline bool RSHash2::extend_in_text(uint64_t &text_pos, uint64_t start, uint64_t end,
    bool forward, const uint64_t query, const uint64_t query_rc)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = get_base(text_pos);
            bool const found = (new_rank == (query >> (2*(k-1))));
            return found;
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = get_base(text_pos);
            bool const found = (new_rank == (query_rc & 0b11));
            return found;
        }
    }
    return false;
}


template<int level>
inline void RSHash2::fill_buffer(uint64_t *offsets, uint64_t *buffer, size_t p, size_t N, const uint64_t mask, const uint64_t shift)
{
    uint64_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    
    for(size_t i = 0; i < N; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i) + 1 - span;
        if constexpr (level == 2)
            offsets[i] = offsets2.access(p+i) + 1 - span;
    }

    for(uint64_t i = 0; i < N; i++) {
        const uint64_t s = offsets[i];

        uint64_t kmer = get_word64(s) & mask;
        uint64_t bits = get_word64(s + k);
        *buffer++ = kmer;
        for(uint64_t j=0; j < span-1; j++) {
            uint64_t const next_base = bits & 3ULL;
            bits >>= 2;
            kmer = (kmer >> 2) | (next_base << shift);
            *buffer++ = kmer;
        }

    }
}


template<int level>
inline bool RSHash2::check_overlap(uint64_t skmer_pos, uint64_t text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    
    const uint64_t r = endpoints.rank(skmer_pos+span);
    start_pos = endpoints.select(r-1, &end_pos);

    return text_pos >= start_pos && text_pos+k-1 < end_pos;
}


template<int level>
inline bool RSHash2::check_minimiser_pos(uint64_t *buffer, const uint64_t offset,
    const uint64_t query, const uint64_t queryrc,
    const size_t s, const size_t e, const size_t minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    if(buffer[s+minimiser_pos] == queryrc) {
        forward = false;
        text_pos = offset + minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[e-1-minimiser_pos] == query) {
        forward = true;
        text_pos = offset + e-1-s-minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }

    return false;
}


template<int level>
inline bool RSHash2::check_minimiser_pos2(uint64_t *buffer, const uint64_t offset,
    const uint64_t query, const uint64_t queryrc,
    const size_t s, const size_t e, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    if(buffer[s+left_minimiser_pos] == queryrc) {
        forward = false;
        text_pos = offset + left_minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[e-1-left_minimiser_pos] == query) {
        forward = true;
        text_pos = offset + e-1-s-left_minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[s+right_minimiser_pos] == query) {
        forward = true;
        text_pos = offset + right_minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[e-1-right_minimiser_pos] == queryrc) {
        forward = false;
        text_pos = offset + e-1-s-right_minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos))
            return true;
    }

    return false;
}


template<int level>
inline bool RSHash2::lookup_buffer(uint64_t* buffer, uint64_t *offsets, const size_t no_skmers,
    const uint64_t query, const uint64_t queryrc,
    uint64_t &text_pos, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &start_pos, uint64_t &end_pos)
{
    size_t span, m;
    if constexpr (level == 1) {
        span = span1;
        m = m1;
    }
    if constexpr (level == 2) {
        span = span2;
        m = m2;
    }

    size_t s = 0, e = 0;
    if(left_minimiser_pos != k-m-right_minimiser_pos) {
        for(size_t i = 0; i < no_skmers; i++) {
            e += span;
            if(check_minimiser_pos2<level>(buffer, offsets[i], query, queryrc, s, e, left_minimiser_pos, right_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s = e;
        }
    }
    else {
        for(size_t i = 0; i < no_skmers; i++) {
            e += span;
            if(check_minimiser_pos<level>(buffer, offsets[i], query, queryrc, s, e, left_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s = e;
        }
    }
    
    return false;
}


template<int level>
inline uint64_t RSHash2::find_minimiser(const uint64_t kmer, const uint64_t kmer_rc, size_t &left_minimiser_position, size_t &right_minimiser_position)
{
    uint64_t m, mmermask;
    mixer_64 m_hasher;
    if constexpr (level == 1) {
        m = m1;
        mmermask = mmermask1;
        m_hasher = m_hasher1;
    }
    if constexpr (level == 2) {
        m = m2;
        mmermask = mmermask2;
        m_hasher = m_hasher2;
    }

    uint64_t mmer = kmer >> 2*(k - m);
    uint64_t mmer_rc = kmer_rc & mmermask;
    uint64_t minimiser = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
    left_minimiser_position = k-m;
    right_minimiser_position = 0;

    for (size_t i = 1; i <= k-m; ++i) {
        mmer = (kmer >> 2*(k-m-i)) & mmermask;
        mmer_rc = (kmer_rc >> 2*i) & mmermask;
        const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
        if(mmerhash < minimiser) {
            minimiser = mmerhash;
            left_minimiser_position = k-m-i;
            right_minimiser_position = i;
        }
        else if(mmerhash == minimiser)
            left_minimiser_position = k-m-i;
    }

    return minimiser;
}

template<int level>
inline void RSHash2::update_minimiser(const uint64_t kmer, const uint64_t kmer_rc, uint64_t &minimiser, size_t &left_minimiser_position, size_t &right_minimiser_position)
{
    uint64_t m, mmermask;
    mixer_64 m_hasher;
    if constexpr (level == 1) {
        m = m1;
        mmermask = mmermask1;
        m_hasher = m_hasher1;
    }
    if constexpr (level == 2) {
        m = m2;
        mmermask = mmermask2;
        m_hasher = m_hasher2;
    }

    if(left_minimiser_position-- == 0) {
        minimiser = find_minimiser<level>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
        return;
    }

    const uint64_t mmer = kmer >> 2*(k - m);
    const uint64_t mmer_rc = kmer_rc & mmermask;
    const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);

    if(mmerhash < minimiser) {
        minimiser = mmerhash;
        left_minimiser_position = k - m;
        right_minimiser_position = 0;
        return;
    }
    if(mmerhash == minimiser) {
        right_minimiser_position = 0;
        return;
    }

    right_minimiser_position++;
}


uint64_t RSHash2::streaming_query(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = rshash::views::kmerview({.window_size = k});

    const uint64_t shift = 2*(k-1);
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* buffer1 = new uint64_t[(m_thres1) * span1];
    uint64_t* buffer2 = new uint64_t[(m_thres2) * span2];
    size_t no_skmers1, no_skmers2;
    uint64_t unitig_begin, unitig_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling1 = false;
    bool rolling2 = false;
    uint64_t occurences = 0;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    size_t left_minimiser2_position, right_minimiser2_position;
    uint64_t minimiser2, minimiser2_rank;

    for(auto && window : query | view)
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, window.kmer_value, window.kmer_value_rev)) {
            occurences++;
            extensions++;
            rolling1 = false;
            rolling2 = false;
        }
        else {
            if(rolling1)
                update_minimiser<1>(window.kmer_value, window.kmer_value_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(window.kmer_value, window.kmer_value_rev, left_minimiser1_position, right_minimiser1_position);
                rolling1 = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, window.kmer_value, window.kmer_value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
                rolling2 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && (minimiser1_rank = r1.rank(minimiser1), r1.rank(minimiser1 + 1) - minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer<1>(offsets1, buffer1, p, no_skmers1, kmermask, shift);
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, window.kmer_value, window.kmer_value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
                current_minimiser1 = minimiser1;
                rolling2 = false;
            }
            else {
                if(rolling2)
                    update_minimiser<2>(window.kmer_value, window.kmer_value_rev, minimiser2, left_minimiser2_position, right_minimiser2_position);
                else {
                    minimiser2 = find_minimiser<2>(window.kmer_value, window.kmer_value_rev, left_minimiser2_position, right_minimiser2_position);
                    rolling2 = true;
                }

                if(minimiser2 == current_minimiser2) {
                    found = lookup_buffer<2>(buffer2, offsets2, no_skmers2, window.kmer_value, window.kmer_value_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, unitig_begin, unitig_end);
                    occurences += found;
                }
                else if(minimiser2 != current_neg_minimiser2 && (minimiser2_rank = r2.rank(minimiser2), r2.rank(minimiser2 + 1) - minimiser2_rank)) {
                    const size_t p = s2_select.select(minimiser2_rank);
                    no_skmers2 = s2_select.select(minimiser2_rank+1) - p;

                    fill_buffer<2>(offsets2, buffer2, p, no_skmers2, kmermask, shift);
                    found = lookup_buffer<2>(buffer2, offsets2, no_skmers2, window.kmer_value, window.kmer_value_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, unitig_begin, unitig_end);
                    occurences += found;
                    current_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                }
                else {
                    occurences += hashmap.contains(std::min<uint64_t>(window.kmer_value, window.kmer_value_rev));
                    found = false;
                    current_neg_minimiser1 = minimiser1;
                    current_neg_minimiser2 = minimiser2;
                }
            }
        }
    }

    delete[] offsets1;
    delete[] offsets2;
    delete[] buffer1;
    delete[] buffer2;
    
    return occurences;
}


int RSHash2::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    cereal::BinaryOutputArchive archive(out);

    archive(k, m1, m2, m_thres1, m_thres2, s1, s2, endpoints, r1, r2, offsets1, offsets2, text, hashmap);

    out.close();
    return 0;
}

int RSHash2::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    cereal::BinaryInputArchive archive(in);

    archive(k, m1, m2, m_thres1, m_thres2, s1, s2, endpoints, r1, r2, offsets1, offsets2, text, hashmap);

    std::cout << "loaded index...\n";
    in.close();

    this->span1 = k - m1 + 1;
    this->span2 = k - m2 + 1;
    this->kmermask = compute_mask(2u * k);
    this->mmermask1 = compute_mask(2u * m1);
    this->mmermask2 = compute_mask(2u * m2);

    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);

    std::cout << "built rank and select ds...\n";
    
    return 0;
}


void RSHash2::print_info() {
    const size_t N = text.size()*32;
    const size_t offset_width = std::bit_width(N);
    const uint64_t M1 = 1ULL << (m1+m1);
    const uint64_t M2 = 1ULL << (m2+m2);
    const uint64_t no_minimizers1 = r1.rank(M1);
    const uint64_t no_skmers1 = s1.size();
    const uint64_t no_minimizers2 = r2.rank(M2);
    const uint64_t no_skmers2 = s2.size();

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "textkmers: " << no_text_kmers <<  '\n';
    
    std::cout << "no minimiser1: " << no_minimizers1 << "\n";
    std::cout << "no distinct minimiser1: " << no_skmers1 << "\n";
    std::cout << "avg superkmers1: " << (double) no_minimizers1/no_skmers1 <<  '\n';
    std::cout << "no minimiser2: " << no_minimizers2 << "\n";
    std::cout << "no distinct minimiser2: " << no_skmers2 << "\n";
    std::cout << "avg superkmers2: " << (double) no_minimizers2/no_skmers2 <<  '\n';
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/no_text_kmers*100 << "%\n";

    std::cout << "density r1: " << (double) no_minimizers1/M1*100 << "%\n";
    std::cout << "density r2: " << (double) no_minimizers2/M2*100 << "%\n";
    std::cout << "density s1: " << (double) s1_select.bitCount()/(no_skmers1+1)*100 <<  "%\n";
    std::cout << "density s2: " << (double) s2_select.bitCount()/(no_skmers2+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/no_text_kmers << "\n";
    std::cout << "endpoints: " << (double) endpoints.bitCount()/no_text_kmers << "\n";
    std::cout << "offsets1: " << (double) no_skmers1*offset_width/no_text_kmers << "\n";
    std::cout << "offsets2: " << (double) no_skmers2*offset_width/no_text_kmers << "\n";
    std::cout << "Hashtable: " << (double) 65*hashmap.bucket_count()/no_text_kmers << "\n";
    std::cout << "R_1: " << (double) r1.bitCount()/no_text_kmers << "\n";
    std::cout << "R_2: " << (double) r2.bitCount()/no_text_kmers << "\n";
    std::cout << "S_1: " << (double) (no_skmers1+1)/no_text_kmers << "\n";
    std::cout << "S_2: " << (double) (no_skmers2+1)/no_text_kmers << "\n";
    
    std::cout << "total: " << (double) (no_skmers1*offset_width+no_skmers2*offset_width+2*N+r1.bitCount()+r2.bitCount()+no_skmers1+1+s1_select.bitCount()+no_skmers2+1+s2_select.bitCount()+endpoints.bitCount()+65*hashmap.bucket_count())/no_text_kmers << "\n";
}