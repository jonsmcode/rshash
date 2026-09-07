#include "rshash.hpp"
#include "minimiser_views.hpp"
#include "kxsort.h"

inline uint64_t mark_sequences(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input, const size_t k,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &endpoints)
{
    size_t length = 0;
    uint64_t kmers = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        length += record.size();
        kmers += record.size() - k + 1;
        no_sequences++;
    }

    std::cout << "text length: " << length << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    std::cout << "no sequences: " << no_sequences << "\n";

    std::cout << "mark endpoints BV...\n";
    bit_vector sequences = bit_vector(length+33, 0);
    sequences[0] = 1;
    sequences[32] = 1;
    uint64_t j = 32;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), length+33);
    sequences = bit_vector();

    return kmers;
}

template<int level, typename MinimizerT>
inline void RSHash::filter_freq_minimizers(std::vector<MinimizerT> &minimizers,
    std::vector<uint16_t> &counts, size_t &no_minimizers, size_t &no_skmers)
{
    uint64_t threshold;
    if constexpr (level == 1)
        threshold = m_thres1;
    if constexpr (level == 2)
        threshold = m_thres2;
    if constexpr (level == 3)
        threshold = m_thres3;
    
    uint64_t current_minimizer = minimizers[0].first;
    uint64_t start = 0;
    size_t write_pos_idx = 0, write_min_idx = 0;
    for(size_t i = 1; i < minimizers.size(); i++) {
        if(minimizers[i].first != current_minimizer) {
            size_t occurences = i - start;
            if(occurences <= threshold) {
                minimizers[write_min_idx++].first = current_minimizer;
                for (size_t j = start; j < i; j++)
                    minimizers[write_pos_idx++].second = minimizers[j].second;
                counts.push_back(static_cast<uint16_t>(occurences));
            }
            current_minimizer = minimizers[i].first;
            start = i;
        }
    }
    size_t occurences = minimizers.size() - start;
    if(minimizers.back().first == current_minimizer && occurences <= threshold) {
        minimizers[write_min_idx++].first = current_minimizer;
        for (size_t j = start; j < minimizers.size(); j++)
            minimizers[write_pos_idx++].second = minimizers[j].second;
        counts.push_back(static_cast<uint16_t>(occurences));
    }

    no_minimizers = write_min_idx;
    no_skmers = write_pos_idx;
}


template<int level, typename MinimizerT>
inline uint64_t RSHash::get_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &skmers,
    std::vector<uint64_t> &unfreq_minimizers, std::vector<uint16_t> &counts)
{
    auto view = [&]() {
        if constexpr (level == 1)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
        if constexpr (level == 2)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m3, .window_size = k, .seed=seed3});
    }();

    std::cout << "computing minimizers and positions...\n";
    std::vector<MinimizerT> minimizers;
    if constexpr (level == 1) {
        uint64_t length = 32;
        for(auto & sequence : sequences) {
            for(auto && minimizer : sequence | view)
                minimizers.emplace_back(MinimizerT{minimizer.minimiser_value, length + minimizer.range_position});
            length += sequence.size();
        }
    }
    else {
        for (auto skmerinfo : skmers) {
            auto skmer = sequences[skmerinfo.seq_id] | std::views::drop(skmerinfo.start) | std::views::take(skmerinfo.end - skmerinfo.start);
            size_t sequence_position = endpoints.select(skmerinfo.seq_id+1);    
            for (auto && minimizer : skmer | view)
                minimizers.emplace_back(MinimizerT{minimizer.minimiser_value, sequence_position + skmerinfo.start + minimizer.range_position});
        }
    }

    std::cout << "sorting minimizers...\n";
    // std::sort(minimizers.begin(), minimizers.end(), [](auto const& a, auto const& b) { return a.minimizer_value < b.minimizer_value; });
    // if(MinimizerT().second <= UINT32_MAX)
    //     kx::radix_sort(minimizers.begin(), minimizers.end(), RadixTraitsMinimizer32());
    // else
        kx::radix_sort(minimizers.begin(), minimizers.end(), RadixTraitsMinimizer64());

    std::cout << "filtering frequent minimizers...\n";
    size_t no_minimizers, no_skmers;
    filter_freq_minimizers<level, MinimizerT>(minimizers, counts, no_minimizers, no_skmers);

    std::cout << "filling minimiser offsets...\n";
    bits::compact_vector::builder builder;
    builder.resize(no_skmers, std::bit_width(endpoints.size()));
    for(size_t i = 0; i < no_skmers; i++)
        builder.push_back(minimizers[i].second);
    
    uint64_t* data_ptr = reinterpret_cast<uint64_t*>(minimizers.data());
    for(size_t i = 0; i < no_minimizers; i++)
        data_ptr[i] = minimizers[i].first;
    auto* raw = reinterpret_cast<uint64_t*>(minimizers.data());
    std::vector<uint64_t> tmp(raw, raw + no_minimizers);
    unfreq_minimizers.swap(tmp);

    minimizers.clear();

    if constexpr (level == 1)
        builder.build(offsets1);
    if constexpr (level == 2)  
        builder.build(offsets2);
    if constexpr (level == 3)
        builder.build(offsets3);

    return no_skmers;
}


inline std::vector<uint64_t> pack_dna4_to_uint64(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    std::vector<uint64_t> packed;
    packed.push_back(UINT64_MAX); // padding

    uint64_t word = 0;
    size_t shift = 0;

    for (uint8_t r : input | std::views::join | seqan3::views::to_rank) {
        word |= uint64_t(r) << shift; // pack 2 bits per base
        shift += 2;

        if (shift == 64) {
            packed.push_back(word);
            word = 0;
            shift = 0;
        }
    }

    if (shift != 0)
        packed.push_back(word);

    packed.push_back(UINT64_MAX); // padding

    return packed;
}

template<int level>
void RSHash::mark_occurences(const size_t total_occs, const std::vector<uint16_t> &occurences)
{
    auto & s = [&]() -> auto& {
    if constexpr (level == 1) return s1;
    if constexpr (level == 2) return s2;
    if constexpr (level == 3) return s3;
    if constexpr (level == 4) return s4;
    if constexpr (level == 5) return s5;
    }();

    s = bit_vector(total_occs+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < occurences.size(); i++) {
        j += occurences[i];
        s[j] = 1;
    }

    if constexpr (level == 1)
        s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 2)
        s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 3)
        s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 4)
        s4_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 5)
        s5_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
}


template<int level>
size_t RSHash::get_frequent_skmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &skmers, std::vector<SkmerInfo> &skmers_out)
{
    auto skmerview = [&]() {
        if constexpr (level == 1)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed = seed1});
        if constexpr (level == 2)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m2, .window_size = k, .seed = seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m3, .window_size = k, .seed = seed3});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level == 1) return r1;
    if constexpr (level == 2) return r2;
    if constexpr (level == 3) return r3;
    }();
    
    size_t freq_kmers = 0;
    if constexpr (level == 1)
    {
        for(size_t sequence_id = 0; sequence_id < sequences.size(); ++sequence_id) {
            auto const & sequence = sequences[sequence_id];
            size_t start_position = 0;
            bool cur_freq, freq;
            uint64_t minimiser_rank;

            for(auto && minimiser : sequence | skmerview) {
                freq = r1.contains(minimiser.minimiser_value, minimiser_rank);
                break;
            }
            for(auto && minimiser : sequence | skmerview) {
                cur_freq = r1.contains(minimiser.minimiser_value, minimiser_rank);
                if(freq && !cur_freq)
                    start_position = minimiser.range_position;
                if(!freq && cur_freq) {
                    skmers_out.emplace_back(SkmerInfo{sequence_id, start_position, minimiser.range_position-1+k});
                    freq_kmers += minimiser.range_position - start_position;
                }
                freq = cur_freq;
            }
            if(!cur_freq) {
                skmers_out.emplace_back(SkmerInfo{sequence_id, start_position, sequence.size()});
                freq_kmers += sequence.size() - start_position - k + 1;
            }
        }
    }
    else
    {
        for(auto skmer : skmers) {
            auto skmer_view = sequences[skmer.seq_id] | std::views::drop(skmer.start) | std::views::take(skmer.end - skmer.start);
            bool cur_freq, freq;
            size_t start = 0;
            uint64_t minimiser_rank;

            for(auto && minimiser : skmer_view | skmerview) {
                freq = r.contains(minimiser.minimiser_value, minimiser_rank);
                break;
            }
            for(auto && minimiser : skmer_view | skmerview) {
                cur_freq = r.contains(minimiser.minimiser_value, minimiser_rank);
                if(freq && !cur_freq)
                    start = minimiser.range_position;
                if(!freq && cur_freq) {
                    skmers_out.emplace_back(SkmerInfo{skmer.seq_id, skmer.start + start, skmer.start + minimiser.range_position-1+k});
                    freq_kmers += minimiser.range_position - start;
                }
                freq = cur_freq;
            }
            if(!cur_freq) {
                skmers_out.emplace_back(SkmerInfo{skmer.seq_id, skmer.start + start, skmer.end});
                freq_kmers += skmer.end - skmer.start - start - k + 1;
            }
        }
    }

    return freq_kmers;
}


gtl::flat_hash_map<uint64_t, uint16_t> RSHash::count_kmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    
    gtl::flat_hash_map<uint64_t, uint16_t> kmer_counts;
    for(const auto & skmer_info : freq_skmers) {
        size_t s = (skmer_info.start > shapes.overlap) ? skmer_info.start - shapes.overlap : 0;
        size_t e = std::min(skmer_info.end + shapes.overlap, input[skmer_info.seq_id].size());
        auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
        for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
            if(number_shapes > 0) {
                // todo: multiple shapes
                uint64_t shape_val = _pext_u64(kmer.value, shapes.shapes[0].mask);
                uint64_t shape_val_rev = _pext_u64(kmer.value_rev, shapes.shapes[0].mask);
                uint64_t canonical_shape = std::min<uint64_t>(shape_val, shape_val_rev);
                if(kmer_counts[canonical_shape] < threshold)
                    kmer_counts[canonical_shape]++;
            }
            else {
                uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                if(kmer_counts[canonical_kmer] < threshold)
                    kmer_counts[canonical_kmer]++;
            }
        }
    }
    return kmer_counts;
}


template <typename AddForward>
void RSHash::process_freq_kmers(AddForward&& add,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers,
    gtl::flat_hash_map<uint64_t, uint16_t> &kmer_counts)
{
    for (const auto& skmer_info : freq_skmers) {
        auto skmer = input[skmer_info.seq_id] | std::views::drop(skmer_info.start) | std::views::take(skmer_info.end - skmer_info.start);
        uint32_t pos = endpoints.select(skmer_info.seq_id + 1) + skmer_info.start;
        for (auto&& kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
            uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
            if(threshold == 0 || kmer_counts[canonical_kmer] < threshold)
                add(canonical_kmer, pos);
            ++pos;
        }
    }
}


template <typename AddForward>
void RSHash::process_freq_kmers(AddForward&& add,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers,
    gtl::flat_hash_map<uint64_t, uint16_t> &kmer_counts, Shape32 &shape)
{
    for (const auto& skmer_info : freq_skmers) {
        size_t s = (skmer_info.start > shape.overlap) ? skmer_info.start - shape.overlap : 0;
        size_t e = std::min(skmer_info.end + shape.overlap, input[skmer_info.seq_id].size());
        auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);

        uint32_t pos = endpoints.select(skmer_info.seq_id + 1) + skmer_info.start;

        for (auto&& kmer : skmer | rshash::views::kmerview({.window_size = shape.length})) {

            uint64_t shape_val     = _pext_u64(kmer.value, shape.mask);
            uint64_t shape_val_rc  = _pext_u64(kmer.value_rev, shape.mask);
            uint64_t canonical_shape = std::min<uint64_t>(shape_val, shape_val_rc);

            add(canonical_shape, pos);

            ++pos;
        }
    }
}


template <size_t MarkId, typename EF, typename Offsets>
void RSHash::build_level(gtl::flat_hash_map<uint64_t, std::vector<uint64_t>>& map, EF& ef, Offsets& offsets)
{
    std::vector<decltype(map.cbegin())> order;
    order.reserve(map.size());

    size_t total_positions = 0;
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        order.push_back(it);
        total_positions += it->second.size();
    }

    std::ranges::sort(order, {}, [](auto it) { return it->first; });

    std::vector<uint64_t> keys;
    std::vector<uint16_t> occurrences;
    bits::compact_vector::builder builder;

    keys.reserve(order.size());
    occurrences.reserve(order.size());
    builder.resize(total_positions, std::bit_width(endpoints.size()));

    for (auto it : order) {
        keys.push_back(it->first);
        occurrences.push_back(it->second.size());

        for (uint32_t pos : it->second)
            builder.push_back(pos);
    }

    map.clear();

    ef = sux::bits::EliasFano(keys, 1ULL << (2*window_size));

    mark_occurences<MarkId>(total_positions, occurrences);
    builder.build(offsets);
}


void RSHash::last_level(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    gtl::flat_hash_map<uint64_t, uint16_t> kmer_counts;
    if(threshold > 0)
        kmer_counts = count_kmers(input, freq_skmers);

    if(use_ht) {
        if(loc) {
            if(number_shapes > 0) {
                for(int i = 0; i < number_shapes; i++) {
                    FlatMapBuilder builder;
                    process_freq_kmers([&](uint64_t key, uint32_t pos) {builder.add(key, pos);},
                                       input, freq_skmers, kmer_counts, shapes.shapes[i]);
                    hashmaps.push_back(std::move(builder).build());
                }
            }
            else {
                FlatMapBuilder builder, builder_rc;
                process_freq_kmers([&](uint64_t key, uint32_t pos) {builder.add(key, pos);},
                                   input, freq_skmers, kmer_counts);
                hashmaps.push_back(std::move(builder).build());
            }
        }
        else {
            if(number_shapes > 0) {
                for(int i = 0; i < number_shapes; i++) {
                    gtl::flat_hash_set<uint64_t> hashset, hashset_rc;
                    process_freq_kmers([&](uint64_t key, uint32_t) {hashset.insert(key);},
                                       input, freq_skmers, kmer_counts, shapes.shapes[i]);
                    hashsets.push_back(hashset);
                }
            }
            else {
                gtl::flat_hash_set<uint64_t> hashset, hashset_rc;
                process_freq_kmers([&](uint64_t key, uint32_t) {hashset.insert(key);},
                                   input, freq_skmers, kmer_counts);
                hashsets.push_back(hashset);
            }
        }
    }
    else {
        if(loc) {
            gtl::flat_hash_map<uint64_t, std::vector<uint64_t>> freq_kmers_map; // todo: 32 bit pos if possible

            process_freq_kmers([&](uint64_t key, uint32_t pos) {freq_kmers_map[key].push_back(pos);},
                               input, freq_skmers, kmer_counts);

            build_level<4>(freq_kmers_map, r4, offsets4);
        }
        else {
            gtl::flat_hash_set<uint64_t> freq_kmers_map;

            process_freq_kmers([&](uint64_t key, uint32_t) {freq_kmers_map.insert(key);},
                               input, freq_skmers, kmer_counts);

            auto build_ef = [&](auto& set, auto& ef) {
                std::vector<uint64_t> keys(set.begin(), set.end());
                set.clear();
                std::ranges::sort(keys);
                ef = sux::bits::EliasFano(keys, 1ULL << (2 * window_size));
            };

            build_ef(freq_kmers_map, r4);
        }
    }

}


void print_freq_skmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> & input,
    const std::vector<SkmerInfo> & freq_skmers)
{
    std::cout << "frequent skmers: " << freq_skmers.size() << "\n\n";
    seqan3::sequence_file_output fout{"../frequent_skmers.fasta"};
    for (size_t i = 0; i < freq_skmers.size(); ++i) {
        auto const & skmer_info = freq_skmers[i];
        // std::cout << ">skmer_" << i << " seq=" << skmer_info.seq_id << " start=" << skmer_info.start << " end=" << skmer_info.end << '\n';
        auto skmer = input[skmer_info.seq_id] | std::views::drop(skmer_info.start) | std::views::take(skmer_info.end - skmer_info.start);
        // for (char c : skmer | seqan3::views::to_char)
        //     std::cout << c;
        // std::cout << '\n';
        std::string id = "seq=" + std::to_string(skmer_info.seq_id) + " start=" + std::to_string(skmer_info.start) + " end=" + std::to_string(skmer_info.end);
        fout.emplace_back(skmer, id);
    }
}


void RSHash::build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input)
{
    no_text_kmers = mark_sequences(input, k, endpoints);
    const size_t text_length = endpoints.size();
    const size_t log_text_length = std::bit_width(text_length);

    std::vector<uint64_t> minimizers1;
    std::vector<uint16_t> minimizers1_occurences;
    std::vector<SkmerInfo> freq_skmers;
    uint64_t no_skmers1;
    if(log_text_length <= 32)
        no_skmers1 = get_minimizers<1, MinimizerInfo32>(input, freq_skmers, minimizers1, minimizers1_occurences);
    else
        no_skmers1 = get_minimizers<1, MinimizerInfo64>(input, freq_skmers, minimizers1, minimizers1_occurences);

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimizers1, M1);
    minimizers1.clear();

    std::cout << "mark skmers1...\n";
    mark_occurences<1>(no_skmers1, minimizers1_occurences);
    minimizers1_occurences.clear();

    std::cout << "get frequent skmers...\n";
    size_t freq_kmers;
    freq_kmers = get_frequent_skmers<1>(input, freq_skmers, freq_skmers);

    if(level > 1)
    {
        std::cout << "count minimizers2...\n";
        std::vector<uint64_t> minimizers2;
        std::vector<uint16_t> minimizers2_occurences;
        uint64_t no_skmers2;
        if(log_text_length <= 32)
            no_skmers2 = get_minimizers<2, MinimizerInfo32>(input, freq_skmers, minimizers2, minimizers2_occurences);
        else
            no_skmers2 = get_minimizers<2, MinimizerInfo64>(input, freq_skmers, minimizers2, minimizers2_occurences);

        std::cout << "build R_2...\n";
        const uint64_t M2 = 1ULL << (m2+m2);
        r2 = sux::bits::EliasFano(minimizers2, M2);
        minimizers2.clear();

        std::cout << "mark skmers2...\n";
        mark_occurences<2>(no_skmers2, minimizers2_occurences);
        minimizers2_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers2;
        freq_kmers = get_frequent_skmers<2>(input, freq_skmers, freq_skmers2);
        freq_skmers.clear();
        freq_skmers = freq_skmers2;
    }

    if(level > 2)
    {
        std::cout << "count minimizers3...\n";
        std::vector<uint64_t> minimizers3;
        std::vector<uint16_t> minimizers3_occurences;
        uint64_t no_skmers3;
        if(log_text_length <= 32)
            no_skmers3 = get_minimizers<3, MinimizerInfo32>(input, freq_skmers, minimizers3, minimizers3_occurences);
        else
            no_skmers3 = get_minimizers<3, MinimizerInfo64>(input, freq_skmers, minimizers3, minimizers3_occurences);

        std::cout << "build R_3...\n";
        const uint64_t M3 = 1ULL << (m3+m3);
        r3 = sux::bits::EliasFano(minimizers3, M3);
        minimizers3.clear();

        std::cout << "mark skmers3...\n";
        mark_occurences<3>(no_skmers3, minimizers3_occurences);
        minimizers3_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers3;
        freq_kmers = get_frequent_skmers<3>(input, freq_skmers, freq_skmers3);
        freq_skmers.clear();
        freq_skmers = freq_skmers3;
    }

    // print_freq_skmers(input, freq_skmers);

    std::cout << "build last level...\n" << freq_kmers << " frequent kmers\n";
    last_level(input, freq_skmers);

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    print_info();
}
