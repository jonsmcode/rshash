#include <filesystem>
#include <bitset>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "rshash.hpp"
#include "io.hpp"
#include "minimiser_views.hpp"


static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if (size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << size) - 1u;
}


RSHash1::RSHash1() : endpoints(std::vector<uint64_t>{}, 1) {}

RSHash1::RSHash1(
    uint8_t const k, uint8_t const m1, uint8_t const m_thres1)
    : k(k), m1(m1),
      m_thres1(m_thres1), span(k-m1+1),
      endpoints(std::vector<uint64_t>{}, 1)
{}


int RSHash1::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto view2 = srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1});
    // auto view1 = srindex::views::xor_minimiser_and_positions2({.minimiser_size = m1, .window_size = k, .seed=seed1});
    // auto view2 = srindex::views::xor_minimiser_and_window2({.minimiser_size = m1, .window_size = k, .seed=seed1});

    std::cout << +m1 << "\n";

    const uint64_t M1 = 1ULL << (m1+m1);

    std::cout << "scan text...\n";
    size_t N = 0;
    uint64_t kmers = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        N += record.size();
        kmers += record.size() - k + 1;
        no_sequences++;
    }

    std::cout << "text length: " << N << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    std::cout << "no sequences: " << no_sequences << "\n";

    std::cout << "get sequences...\n";
    sequences = bit_vector(N+1, 0);
    size_t j = 0;
    sequences[0] = 1;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), N+1);

    std::cout << "count minimizers1...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers1;

    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            minimizers1[minimiser.minimiser_value] += (minimiser.occurrences+span-1)/span;
            if(minimizers1[minimiser.minimiser_value] > m_thres1)
                minimizers1[minimiser.minimiser_value] = m_thres1;
            n += (minimiser.occurrences+span-1)/span;
        }
    }
    
    std::cout << "extract unfrequent minimizers...\n";
    uint64_t n1 = 0;
    std::vector<uint64_t> unfreq_minimizers1;
    for(auto const& [minimizer, count] : minimizers1) {
        if(count < m_thres1) {
            unfreq_minimizers1.push_back(minimizer);
            n1 += count;
        }
    }

    size_t c1 = minimizers1.size();
    size_t c1tmp = unfreq_minimizers1.size();

    std::cout << "unfrequent minimizers: " << unfreq_minimizers1.size() << " (" << (double) unfreq_minimizers1.size()/minimizers1.size()*100 << "%)\n";

    std::cout << "build R_1...\n";
    sd_vector_builder builder(M1, unfreq_minimizers1.size());

    std::sort(unfreq_minimizers1.begin(), unfreq_minimizers1.end());
    for(uint64_t minimizer : unfreq_minimizers1)
        builder.set(minimizer);

    r1 = sd_vector<>(builder);
    r1_rank = rank_support_sd<>(&r1);

    std::cout << "filling bitvector S_1...\n";
    s1 = bit_vector(n1+1, 0);
    s1[0] = 1;
    j = 0;
    for(uint64_t minimizer : unfreq_minimizers1) {
        j += minimizers1[minimizer];
        s1[j] = 1;
    }
    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), n1+1, 3);

    minimizers1.clear();
    unfreq_minimizers1.clear();

    std::cout << "filling offsets_1...\n";
    const size_t offset_width = std::bit_width(N);
    pthash::compact_vector::builder b1;
    b1.resize(n1, offset_width);

    uint8_t* count1 = new uint8_t[c1tmp];
    std::memset(count1, 0, c1tmp*sizeof(uint8_t));

    size_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                size_t s = s1_select.select(i);
                b1.set(s + count1[i], length + minimiser.range_position);
                count1[i]++;
            }
        }
        length += sequence.size();
    }
    b1.build(offsets1);

    delete[] count1;

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers1;
    std::vector<size_t> skmer_positions;
    length = 0;
    for(auto & sequence : input) {
        size_t start_position = 0;
        bool cur_freq, freq;

        for(auto && minimiser : sequence | view1) {
            freq = r1[minimiser.minimiser_value];
            break;
        }
        for(auto && minimiser : sequence | view1) {
            cur_freq = r1[minimiser.minimiser_value];
            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position-1+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers1.push_back(skmer);
                skmer_positions.push_back(length + start_position);
            }
            freq = cur_freq;
        }
        if(!cur_freq) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers1.push_back(skmer);
            skmer_positions.push_back(length + start_position);
        }

        length += sequence.size();
    }

    // uint64_t freq_kmers = 0;
    // for(auto & sequence : input) {
    //     for(auto && minimiser : sequence | view2) {
    //         freq_kmers += !r1[minimiser.minimiser_value];
    //     }
    // }
    // std::cout << "frequent k-mers: " << freq_kmers << " (" << (double) freq_kmers/kmers*100 << "%)\n";
    
    size_t rem_kmers1 = 0;
    for(auto & skmer : freq_skmers1)
        rem_kmers1 += skmer.size() - k + 1;
    // std::cout << "skmers: " << skmers << " minimizer: " << n << " unfrequent minimizers: " << unfreq_minimizers1.size() << " (" << (double) unfreq_minimizers1.size()/minimizers1.size()*100 << "%) \n";
    // std::cout << "remaining superkmers " << freq_skmers1.size() << " (" << (double) freq_skmers1.size()/skmers*100 << "%) ";
    std::cout << "remaining kmers: " << rem_kmers1 << " (" << (double) rem_kmers1/kmers*100 << "%)\n";

    std::cout << "build level 4, HT...\n";

    // todo: simple kmer view
    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view2) {
            hashmap.insert(std::min<uint64_t>(minimiser.window_value, minimiser.window_value_rev));
        }
    }

    std::cout << "copy text...\n";
    for(auto & record : input) {
        std::ranges::move(record, std::back_inserter(text));
    }

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    
    std::cout << "no minimiser: " << n << "\n";
    std::cout << "no distinct minimiser: " << c1tmp << "\n";
    std::cout << "minimiser going to level 2: " << c1-c1tmp << "  " << (double) (c1-c1tmp)/c1*100 << "%\n";
    std::cout << "no minimiser1: " << n1 << "\n";
    std::cout << "no distinct minimiser1: " << c1tmp << "\n";
    std::cout << "avg superkmers1: " << (double) n1/c1tmp <<  '\n';
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/kmers*100 << "%\n";

    std::cout << "density r1: " << (double) c1/M1*100 << "%\n";
    std::cout << "density s1: " << (double) s1_select.bitCount()/(n1+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/kmers << "\n";
    std::cout << "endpoints: " << (double) endpoints.bitCount()/kmers << "\n";
    std::cout << "offsets1: " << (double) n1*offset_width/kmers << "\n";
    // std::cout << "Hashtable: " << (double) 64*hashmap.size()/kmers << "\n";
    std::cout << "Hashtable: " << (double) 65*hashmap.bucket_count()/kmers << "\n";
    std::cout << "R_1: " << (double) (8*size_in_bytes(r1))/kmers << "\n";
    std::cout << "S_1: " << (double) (n1+1 + s1_select.bitCount())/kmers << "\n";


    std::cout << "total: " << (double) (n1*offset_width+2*N+8*size_in_bytes(r1)+n1+1+s1_select.bitCount()+endpoints.bitCount()+65*hashmap.bucket_count())/kmers << "\n";

    return 0;
}



static inline constexpr uint64_t crc(uint64_t x, uint64_t k) {
    // assert(k <= 32);
    uint64_t c = ~x;

    /* swap byte order */
    uint64_t res = __builtin_bswap64(c);

    /* Swap nuc order in bytes */
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;              // ...0000.1111.0000.1111
    const uint64_t c2 = 0x3333333333333333;              // ...0011.0011.0011.0011
    res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4);  // swap 2-nuc order in bytes
    res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2);  // swap nuc order in 2-nuc

    /* Realign to the right */
    res >>= 64 - 2 * k;

    return res;
}


std::vector<uint64_t> RSHash1::rand_text_kmers(const uint64_t n) {
    std::uniform_int_distribution<uint32_t> distr;
    std::mt19937 m_rand(1);
    std::vector<std::uint64_t> kmers;
    kmers.reserve(n);

    const uint64_t no_unitigs = number_unitigs();
    for (uint64_t i = 0; i < n; ++i) {
        const uint64_t unitig_id = distr(m_rand) % no_unitigs;
        const uint64_t offset = distr(m_rand) % unitig_size(unitig_id);
        const uint64_t kmer = access(unitig_id, offset);

        if ((i & 1) == 0)
            kmers.push_back(crc(kmer, k));
        else
            kmers.push_back(kmer);
    }

    return kmers;
}



uint64_t RSHash1::access(const uint64_t unitig_id, const size_t offset)
{
    size_t offset_text = endpoints.select(unitig_id) + offset;

    uint64_t kmer = 0;
    for (size_t i=offset_text; i < offset_text+k; i++) {
        uint64_t const new_rank = seqan3::to_rank(text[i]);
        kmer <<= 2;
        kmer |= new_rank;
    }

    return kmer;
}


inline bool RSHash1::check(const size_t p, const size_t q, const uint64_t mask,
    const uint64_t kmer, const uint64_t kmer_rc,
    double &to, double &th, double &te)
{
    std::chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4;
    for(size_t i = 0; i < q-p; i++) {

        t0 = std::chrono::high_resolution_clock::now();
        size_t o = offsets1.access(p+i);
        t1 = std::chrono::high_resolution_clock::now();
        to += (std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count();

        uint64_t hash = 0;
        for (size_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
        }
        t2 = std::chrono::high_resolution_clock::now();
        th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
        if(hash == kmer || hash == kmer_rc)
            return true;

        uint64_t e = std::min<uint64_t>(endpoints.select(endpoints.rank(o+1)), o+k+span);
        // uint64_t e = o+k+span;

        t3 = std::chrono::high_resolution_clock::now();
        te += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();

        for(size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
            if(hash == kmer || hash == kmer_rc) {
                t4 = std::chrono::high_resolution_clock::now();
                th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count();
                return true;
            }
        }
        t4 = std::chrono::high_resolution_clock::now();
        th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count();
    }

    return false;
}


inline bool RSHash1::check(const size_t p, const size_t q, const uint64_t mask,
    const uint64_t kmer, const uint64_t kmer_rc)
{
    std::chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4;
    for(size_t i = 0; i < q-p; i++) {

        t0 = std::chrono::high_resolution_clock::now();
        size_t o = offsets1.access(p+i);

        uint64_t hash = 0;
        for (size_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
        }
        if(hash == kmer || hash == kmer_rc)
            return true;

        uint64_t e = std::min<uint64_t>(endpoints.select(endpoints.rank(o+1)), o+k+span);
        // uint64_t e = o+k+span;

        for(size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
            if(hash == kmer || hash == kmer_rc) {
                return true;
            }
        }
    }

    return false;
}

uint64_t RSHash1::lookup(const std::vector<uint64_t> &kmers, bool verbose)
{
uint64_t occurences = 0;
if (verbose) {
    const uint64_t mask = compute_mask(2u * k);
    srindex::minimizers::Minimisers_hash2 minimisers = srindex::minimizers::Minimisers_hash2(k, m1, seed1);
    std::chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4, t5, t6 = std::chrono::high_resolution_clock::now();
    double t0_ = 0.0;
    double t1_ = 0.0;
    double t2_ = 0.0;
    double t3_ = 0.0;
    double t4_ = 0.0;
    double t5_ = 0.0;
    double to = 0.0;
    double th = 0.0;
    double te = 0.0;
    size_t skmers1_ = 0;
    uint64_t lookups1 = 0;
    uint64_t ht_lookups = 0;

    for(uint64_t kmer : kmers)
    {
        t5 = std::chrono::high_resolution_clock::now();
        minimisers.compute(kmer);
        t4_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t5)).count();

        t6 = std::chrono::high_resolution_clock::now();
        if(r1[minimisers.minimiser1]) {
            t0 = std::chrono::high_resolution_clock::now();
            uint64_t minimizer_id = r1_rank(minimisers.minimiser1);
            t1 = std::chrono::high_resolution_clock::now();
            size_t p = s1_select.select(minimizer_id);
            size_t q = s1_select.select(minimizer_id+1);
            t2 = std::chrono::high_resolution_clock::now();

            occurences += check(p, q, mask, minimisers.window, minimisers.window_rev, to, th, te);
            // occurences += check<1>(p, q, mask, minimisers.window, minimisers.window_rev);
            t3 = std::chrono::high_resolution_clock::now();

            t0_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count();
            t1_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
            t2_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();
            t5_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t0 - t6)).count();
            skmers1_ += q - p;
            ++lookups1;
        }
        else {
            t4 = std::chrono::high_resolution_clock::now();
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window, minimisers.window_rev));
            t3_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t4)).count();
            t5_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t6)).count();
            ++ht_lookups;
        }
    }
    std::cout << "minimisers: " << t4_/kmers.size() << " ns\n";
    std::cout << "R lookup: " << t5_/kmers.size() << " ns\n";
    std::cout << "r_rank: " << t0_/kmers.size() << " ns\n";
    std::cout << "s_select: " << t1_/kmers.size() << " ns\n";
    std::cout << "check: " << t2_/kmers.size() << " ns\n";
    std::cout << "offsets: " << to/kmers.size() << " ns\n";
    std::cout << "endpoints: " << te/kmers.size() << " ns\n";
    std::cout << "text: " << th/kmers.size() << " ns\n";
    std::cout << "ht: " << t3_/kmers.size() << " ns\n";
    std::cout << "lookups lvl1: " << (double) lookups1/(lookups1+ht_lookups)*100 << "%\n";
    std::cout << "lookups ht: " << (double) ht_lookups/(lookups1+ht_lookups)*100 << "%\n";
    std::cout << "avg skmers1: " << (double) skmers1_/lookups1 << "\n";
}
else {
    const uint64_t mask = compute_mask(2u * k);
    srindex::minimizers::Minimisers_hash2 minimisers = srindex::minimizers::Minimisers_hash2(k, m1, seed1);

    for(uint64_t kmer : kmers)
    {
        minimisers.compute(kmer);

        if(r1[minimisers.minimiser1]) {
            uint64_t minimizer_id = r1_rank(minimisers.minimiser1);
            size_t p = s1_select.select(minimizer_id);
            size_t q = s1_select.select(minimizer_id+1);

            occurences += check(p, q, mask, minimisers.window, minimisers.window_rev);
        }
        else {
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window, minimisers.window_rev));
        }

    }
}

    return occurences;
}




inline bool RSHash1::extend_in_text(size_t &text_pos, size_t start, size_t end,
    bool forward, const uint64_t query, const uint64_t query_rc, uint64_t &fwd_extensions, uint64_t &rev_extensions)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = seqan3::to_rank(text[text_pos]);
            bool const found = (new_rank == (query >> (2*(k-1))));
            fwd_extensions++;
            return found;
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = seqan3::to_rank(text[text_pos]);
            bool const found = (new_rank == (query_rc & 0b11));
            rev_extensions++;
            return found;
        }
    }
    return false;
}


inline void RSHash1::fill_buffer(std::vector<uint64_t> &buffer, std::vector<SkmerInfo> &skmers, size_t p, size_t q)
{
    for(size_t i = 0; i < q-p; i++) {
        size_t o = offsets1.access(p+i);
        uint64_t next_endpoint;
        uint64_t prev_endpoint = endpoints.select(endpoints.rank(o+1)-1, &next_endpoint);
        size_t e = std::min(o+k+span, next_endpoint);

        // todo: get 256 chars around o in one cache line
        // todo: let buffer be bits text[max(o-126,s),...,o,...,min(o+126,e)]
        // check kmer at position p in skmer with reinterpret_cast<uint64_t*>(buffer+2*p)[0] >> (64 - 2*k);
        uint64_t hash = 0;
        for (size_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
        }
        buffer.push_back(hash);
        for (size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
            buffer.push_back(hash);
        }
        skmers.push_back({o, e - o - k + 1, prev_endpoint, next_endpoint});
    }
}


inline bool RSHash1::lookup_serial(std::vector<uint64_t> &buffer, std::vector<SkmerInfo> &skmers,
    const uint64_t query, const uint64_t queryrc,
    size_t &text_pos, size_t &skmer_pos, bool &forward, size_t &kmers_left, size_t &start_pos, size_t &end_pos)
{
    size_t s = 0, e = 0;
    for (size_t i = 0; i < skmers.size(); i++) {
        e += skmers[i].length;
        for(size_t j = s; j < e; j++) {
            if(buffer[j] == query) {
                forward = true;
                text_pos = skmers[i].position + (j - s) + k - 1;
                end_pos = skmers[i].unitig_end;
                return true;
            }
            if(buffer[j] == queryrc) {
                forward = false;
                text_pos = skmers[i].position + (j - s);
                start_pos = skmers[i].unitig_begin;
                return true;
            }
        }
        s = e;
    }
    return false;
}


uint64_t RSHash1::streaming_query(const std::vector<seqan3::dna4> &query,
    uint64_t &buffer_fwd_extensions, uint64_t &buffer_rev_extensions, uint64_t &text_fwd_extensions, uint64_t &text_rev_extensions)
{
    auto view = srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1});

    uint64_t occurences = 0;
    uint64_t current_minimiser=std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> buffer;
    std::vector<SkmerInfo> skmers;
    size_t unitig_begin, unitig_end;
    size_t skmer_pos;
    size_t kmers_left = 0;
    size_t text_pos;
    bool forward;
    bool found = false;

    for(auto && minimisers : query | view)
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, minimisers.window_value, minimisers.window_value_rev, text_fwd_extensions, text_rev_extensions))
            occurences++;
        else if(minimisers.minimiser_value == current_minimiser) {
            found = lookup_serial(buffer, skmers, minimisers.window_value, minimisers.window_value_rev, text_pos, skmer_pos, forward, kmers_left, unitig_begin, unitig_end);
            occurences += found;
        }
        else if(r1[minimisers.minimiser_value]) {
            uint64_t minimizer_id = r1_rank(minimisers.minimiser_value);
            size_t p = s1_select.select(minimizer_id);
            size_t q = s1_select.select(minimizer_id+1);

            buffer.clear();
            skmers.clear();
            fill_buffer(buffer, skmers, p, q);
            found = lookup_serial(buffer, skmers, minimisers.window_value, minimisers.window_value_rev, text_pos, skmer_pos, forward, kmers_left, unitig_begin, unitig_end);
            occurences += found;
            current_minimiser = minimisers.minimiser_value;
        }
        else {
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window_value, minimisers.window_value_rev));
            found = false;
        }
    }
    
    return occurences;
}


int RSHash1::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m1, out);
    seqan3::contrib::sdsl::serialize(r1, out);
    seqan3::contrib::sdsl::serialize(s1, out);
    seqan3::contrib::sdsl::serialize(this->sequences, out);

    cereal::BinaryOutputArchive archive(out);
    archive(this->offsets1);
    archive(this->text);
    archive(this->hashmap);

    out.close();
    return 0;
}

int RSHash1::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m1, in);
    seqan3::contrib::sdsl::load(r1, in);
    seqan3::contrib::sdsl::load(s1, in);
    seqan3::contrib::sdsl::load(this->sequences, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->offsets1);
    archive(this->text);
    archive(this->hashmap);

    this->span = k - m1 + 1;

    std::cout << "loaded index...\n";

    in.close();

    r1_rank = rank_support_sd<>(&r1);
    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);

    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), sequences.size());

    std::cout << "built rank and select ds...\n";
    
    return 0;
}

