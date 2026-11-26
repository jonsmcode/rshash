#include "rshash.hpp"
#include "minimiser_views.hpp"
#include "kxsort.h"
#include <optional>
#include <thread>


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
    uint64_t max_occurences;
    if constexpr (level == 1)
        max_occurences = m_thres1;
    if constexpr (level == 2)
        max_occurences = m_thres2;
    if constexpr (level == 3)
        max_occurences = m_thres3;

    // k-mer/shape frequency threshold: a k-mer/shape stored on this level must occur less than threshold times, so the
    // minimisers occurring threshold times or more are passed on to the next level; the k-mers/shapes that reach the
    // last level are counted exactly (see last_level). the build view emits all positions of a minimiser in a window
    // (ties), so every occurrence of a k-mer/shape value is found from its own position of the minimiser: the leftmost
    // one on its strand (forward), the mirrored one (reverse complement). a forward and a reverse complement occurrence
    // could only share this position if it held an m-mer equal to its reverse complement, which is impossible for odd
    // m (required with a threshold, see build): frequency <= minimiser occurrences
    size_t no_threshold_minimizers = 0;
    auto exceeds_threshold = [&](const size_t occurences) {
        return threshold > 0 && occurences >= threshold;
    };

    size_t write_pos_idx = 0, write_min_idx = 0;

    // keeps the distinct positions of a minimiser if there are at most threshold of them; the radix sort does not
    // order positions within small buckets, and at levels > 1 two super-k-mers can report the same position
    auto keep_minimizer = [&](const size_t start, const size_t end) {
        const uint64_t minimizer = minimizers[start].first;
        std::sort(minimizers.begin() + start, minimizers.begin() + end, [](auto const &a, auto const &b) { return a.second < b.second; });
        size_t occurences = 0;
        for (size_t j = start; j < end; j++)
            occurences += (j == start || minimizers[j].second != minimizers[j-1].second);
        if(occurences > max_occurences)
            return;
        if(exceeds_threshold(occurences)) {
            no_threshold_minimizers++;
            return;
        }

        minimizers[write_min_idx++].first = minimizer;
        uint64_t previous = 0;
        for (size_t j = start; j < end; j++) {
            const uint64_t position = minimizers[j].second;
            if (j == start || position != previous)
                minimizers[write_pos_idx++].second = position;
            previous = position;
        }
        counts.push_back(static_cast<uint16_t>(occurences));
    };

    size_t start = 0;
    for(size_t i = 1; i < minimizers.size(); i++) {
        if(minimizers[i].first != minimizers[start].first) {
            keep_minimizer(start, i);
            start = i;
        }
    }
    keep_minimizer(start, minimizers.size());

    if(no_threshold_minimizers > 0)
        std::cout << no_threshold_minimizers << " minimizers passed on because of the k-mer threshold " << threshold << "\n";

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

template <typename Count>
static void mark_occurences(const size_t total_occs, const std::vector<Count> &occurences,
    bit_vector &s, sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s_select)
{
    s = bit_vector(total_occs+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < occurences.size(); i++) {
        j += occurences[i];
        s[j] = 1;
    }

    s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
}

template<int level>
void RSHash::mark_occurences(const size_t total_occs, const std::vector<uint16_t> &occurences)
{
    if constexpr (level == 1)
        ::mark_occurences(total_occs, occurences, s1, s1_select);
    if constexpr (level == 2)
        ::mark_occurences(total_occs, occurences, s2, s2_select);
    if constexpr (level == 3)
        ::mark_occurences(total_occs, occurences, s3, s3_select);
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


// calls f(window, text position) for every shape window (of shape.length bases) of the frequent super-k-mers
template <typename F>
void RSHash::for_each_shape_window(F&& f, const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers, const Shape64 &shape)
{
    // a super-k-mer (of kernels) is extended by the overlap of the shape window on both sides: the kernel lies at the
    // start of the run of a forward shape and at the end of the run of a reverse shape, which is shapes.overlap
    // (>= shape.overlap + run - kernel) bases away from the shape ends. so the windows of neighbouring super-k-mers
    // can overlap: every window (text position) is only processed once, the super-k-mers are sorted by text position
    const size_t extension = shapes.overlap;
    uint64_t next_pos = 0, last_start = 0;
    for (const auto& skmer_info : freq_skmers) {
        const uint64_t sequence_begin = endpoints.select(skmer_info.seq_id + 1);
        size_t s = (skmer_info.start > extension) ? skmer_info.start - extension : 0;
        const size_t e = std::min(skmer_info.end + extension, input[skmer_info.seq_id].size());
        if(sequence_begin + s < last_start)
            throw std::runtime_error("super-k-mers of the last level are not sorted");
        last_start = sequence_begin + s;
        s = std::max<uint64_t>(sequence_begin + s, next_pos) - sequence_begin;
        if(s + shape.length > e)
            continue;
        // text position of the current shape window
        uint64_t pos = sequence_begin + s;
        next_pos = sequence_begin + e - shape.length + 1;

        auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);

        for (auto&& kmer : skmer | rshash::views::kmer_view({.window_size = shape.length}))
            f(kmer, pos++);
    }
}


template <typename AddForward, typename AddReverse, typename Keep>
void RSHash::process_freq_kmers(AddForward&& addfwd, AddReverse&& addrev, Keep&& keep,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers, const Shape64 &shape)
{
    constexpr bool filter = !std::is_same_v<std::decay_t<Keep>, KeepAll>;
    for_each_shape_window([&](const auto &kmer, const uint64_t pos) {
        const uint64_t shape_fwd = _pext_u64(kmer.value, shape.mask.lo) | (_pext_u64(kmer.value_hi, shape.mask.hi) << shape.lo_weight);
        const uint64_t shape_rev = _pext_u64(kmer.value, shape.mask_rev.lo) | (_pext_u64(kmer.value_hi, shape.mask_rev.hi) << shape.rev_lo_weight);
        if constexpr (filter) {
            // the k-mer/shape threshold is applied to the shape values of both strands: the forward structure
            // stores the shape of the window, the reverse structure (reversed shape) the shape of its reverse
            // complement, which is found by the queries with this shape value
            const uint64_t shape_of_rc = _pext_u64(kmer.value_rev, shape.mask.lo) | (_pext_u64(kmer.value_rev_hi, shape.mask.hi) << shape.lo_weight);
            if(keep(shape_fwd))
                addfwd(shape_fwd, pos);
            if(keep(shape_of_rc))
                addrev(shape_rev, pos);
        }
        else {
            addfwd(shape_fwd, pos);
            addrev(shape_rev, pos);
        }
    }, input, freq_skmers, shape);
}


// symmetric shapes (see symmetric_placement): canonical values (minimum of the shape on the window and on its reverse
// complement) with their text positions. the threshold of the canonical value is the one of both values (the counts
// of a value and of its reverse complement are the same: forward plus reverse occurrences)
template <typename Add, typename Keep>
void RSHash::process_canonical_shapes(Add&& add, Keep&& keep,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers, const Shape64 &shape)
{
    constexpr bool filter = !std::is_same_v<std::decay_t<Keep>, KeepAll>;
    for_each_shape_window([&](const auto &kmer, const uint64_t pos) {
        const uint64_t shape_fwd = _pext_u64(kmer.value, shape.mask.lo) | (_pext_u64(kmer.value_hi, shape.mask.hi) << shape.lo_weight);
        const uint64_t shape_of_rc = _pext_u64(kmer.value_rev, shape.mask.lo) | (_pext_u64(kmer.value_rev_hi, shape.mask.hi) << shape.lo_weight);
        const uint64_t canonical = std::min(shape_fwd, shape_of_rc);
        if constexpr (filter) {
            if(!keep(canonical))
                return;
        }
        // add can also take the orientation (whether the canonical value is the forward value)
        if constexpr (std::is_invocable_v<Add&, uint64_t, uint32_t, bool>)
            add(canonical, pos, shape_fwd <= shape_of_rc);
        else
            add(canonical, pos);
    }, input, freq_skmers, shape);
}


// canonical k-mers (minimum of the k-mer and its reverse complement) with their text positions, k <= 32
template <typename Add, typename Keep>
void RSHash::process_canonical_kmers(Add&& add, Keep&& keep,
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    constexpr bool filter = !std::is_same_v<std::decay_t<Keep>, KeepAll>;
    for (const auto& skmer_info : freq_skmers) {
        uint32_t pos = endpoints.select(skmer_info.seq_id + 1) + skmer_info.start;
        auto skmer = input[skmer_info.seq_id] | std::views::drop(skmer_info.start) | std::views::take(skmer_info.end - skmer_info.start);

        for (auto&& kmer : skmer | rshash::views::kmer_view({.window_size = window_size})) {
            const uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
            if constexpr (filter) {
                if(!keep(canonical_kmer)) {
                    ++pos;
                    continue;
                }
            }
            // add can also take the orientation (whether the canonical k-mer is the forward k-mer)
            if constexpr (std::is_invocable_v<Add&, uint64_t, uint32_t, bool>)
                add(canonical_kmer, pos, kmer.value <= kmer.value_rev);
            else
                add(canonical_kmer, pos);
            ++pos;
        }
    }
}


// keys have to be sorted
static inline sux::bits::EliasFano<sux::util::AllocType::MALLOC> build_ef(const std::vector<uint64_t> &keys)
{
    // shape values can use all 64 bits, so the universe is bounded by the largest key
    // todo: a key of UINT64_MAX (e.g. poly-T for weight 32) cannot be represented
    const uint64_t universe = keys.empty() ? 1 : std::max(keys.back(), keys.back() + 1);
    return sux::bits::EliasFano(keys, universe);
}


// fills the select bitvector and offsets of structure idx from pairs sorted by sort_unique and returns its rank structure
RSHash::EF RSHash::build_level(const size_t idx, const std::vector<KeyValue>& pairs)
{
    std::vector<uint64_t> keys;
    std::vector<uint64_t> occurrences; // k-mers in repeats can occur more than 2^16 times
    bits::compact_vector::builder builder;
    builder.resize(pairs.size(), std::bit_width(endpoints.size()));

    for (size_t i = 0; i < pairs.size(); i++) {
        if (i == 0 || pairs[i].key != pairs[i-1].key) {
            keys.push_back(pairs[i].key);
            occurrences.push_back(0);
        }
        occurrences.back()++;
        builder.push_back(pairs[i].value);
    }

    ::mark_occurences(pairs.size(), occurrences, si[idx], si_select[idx]);
    builder.build(offsetsi[idx]);

    return build_ef(keys);
}


void RSHash::last_level(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    if(!store_last_level)
        return;

    // shapes: a forward and a reverse structure per shape (with hash tables one canonical structure per symmetric
    // shape, see number_symmetric), k-mers: one structure of the canonical k-mers
    const int no_shapes = std::max(number_shapes, 1);

    // calls add/add_rc for the (forward/reverse) values of shape i (the canonical k-mers) that keep accepts
    auto process_filtered = [&](const int i, auto &&add, auto &&add_rc, auto &&keep) {
        if(number_shapes > 0)
            process_freq_kmers(add, add_rc, keep, input, freq_skmers, shapes.shapes[i]);
        else
            process_canonical_kmers(add, keep, input, freq_skmers);
    };

    // k-mer/shape threshold: the frequencies (number of text positions of the canonical k-mer, of the shape value on
    // both strands) of the k-mers/shapes on the last level, counted up to threshold. all occurrences of a k-mer/shape
    // with frequent minimisers on all levels are on the last level, so these are their exact frequencies
    std::vector<gtl::flat_hash_map<uint64_t, uint16_t>> counts(threshold > 0 ? no_shapes : 0);
    if(threshold > 0) {
        std::cout << "count last level k-mers/shapes...\n";
        auto skip = [](uint64_t, uint32_t) {};
        auto count_shape = [&](const int i) {
            auto &count = counts[i];
            process_filtered(i, skip, skip, [&](const uint64_t value) {
                uint16_t &c = count[value];
                c += (c < threshold);
                return false;
            });
        };
        if(parallel_last_level) {
            std::vector<std::jthread> workers;
            for(int i = 1; i < no_shapes; i++)
                workers.emplace_back([&count_shape, i]() { count_shape(i); });
            count_shape(0);
        }
        else
            for(int i = 0; i < no_shapes; i++)
                count_shape(i);

        for(int i = 0; i < no_shapes; i++) {
            size_t frequent = 0;
            for(const auto &[value, c] : counts[i])
                frequent += (c >= threshold);
            std::cout << (number_shapes > 0 ? "shape " + std::to_string(shape_order[i]) + ": " : "") << frequent << " of " << counts[i].size()
                      << " distinct last level k-mers/shapes with frequency >= " << threshold << " are not stored\n";
        }
    }

    // calls f(keep) with the keep predicate of shape i (the k-mers)
    auto with_keep = [&](const int i, auto &&f) {
        if(threshold > 0) {
            const auto &count = counts[i];
            f([&count, t = threshold](const uint64_t value) {
                const auto it = count.find(value);
                return it == count.end() || it->second < t;
            });
        }
        else
            f(KeepAll{});
    };
    auto process = [&](const int i, auto &&add, auto &&add_rc) {
        with_keep(i, [&](auto &&keep) { process_filtered(i, add, add_rc, keep); });
    };

    // the shape (0 for k-mers) of structure j and which of its values the structure holds: first the canonical
    // structures (the k-mers or the symmetric shapes), then the forward and reverse structure of each other shape
    enum class Kind { forward, reverse, canonical };
    std::vector<std::pair<int, Kind>> structures;
    if(number_shapes == 0)
        structures.emplace_back(0, Kind::canonical);
    for(int i = 0; i < number_symmetric; i++)
        structures.emplace_back(i, Kind::canonical);
    for(int i = number_symmetric; i < number_shapes; i++) {
        structures.emplace_back(i, Kind::forward);
        structures.emplace_back(i, Kind::reverse);
    }
    const size_t no_structures = structures.size();
    const size_t no_canonical = (number_shapes == 0) ? 1 : number_symmetric;

    // calls add(key, position[, forward]) for the values of structure j: the forward or reverse values of a shape,
    // the canonical values of a symmetric shape or the canonical k-mers
    auto process_structure = [&](const size_t j, auto &&add) {
        const int i = structures[j].first;
        const Kind kind = structures[j].second;
        auto skip = [](uint64_t, uint32_t, bool = true) {};
        if(kind == Kind::forward)
            process(i, add, skip);
        else if(kind == Kind::reverse)
            process(i, skip, add);
        else if(number_shapes == 0)
            process(i, add, skip);
        else
            with_keep(i, [&](auto &&keep) { process_canonical_shapes(add, keep, input, freq_skmers, shapes.shapes[i]); });
    };

    // with locate: the (key, position) pairs of structure j are collected and sorted instead of hashed one by one
    auto collect_sorted = [&](const size_t j) {
        std::vector<KeyValue> pairs;
        process_structure(j, [&](uint64_t key, uint32_t pos) { pairs.push_back({key, pos}); });
        sort_unique(pairs);
        return pairs;
    };

    // runs task(j) for all structures, one thread each if parallel_last_level
    auto for_each_structure = [&](auto &&task) {
        if(!parallel_last_level) {
            for(size_t j = 0; j < no_structures; j++)
                task(j);
            return;
        }
        std::vector<std::jthread> workers;
        for(size_t j = 1; j < no_structures; j++)
            workers.emplace_back([&task, j]() { task(j); });
        task(0);
    };

    if(use_ht) {
        if(loc) {
            std::vector<FlatMap> maps(no_structures);
#ifdef FLATMAP_ENCODED
            // positional information as in the u2o table of piscem-style rshash (encode_entry of ContigTable):
            // (sequence id << (ref_len_bits + 1)) | (position in sequence << 1) | forward, with
            // ref_len_bits = ceil(log2(max sequence length + 1)) and ceil(log2(no sequences + 1)) bits for the id
            uint64_t max_len = 0;
            for(const auto &record : input)
                max_len = std::max<uint64_t>(max_len, record.size());
            const uint64_t ref_len_bits = std::bit_width(max_len);
            const uint64_t ref_id_bits = std::bit_width(input.size());
            const uint64_t width = ref_len_bits + ref_id_bits + 1;
            if(width > 57) // compact_vector::access() reads at most 57 bits
                throw std::runtime_error("encoded positions need more than 57 bits");
            // canonical k-mers/shapes: the orientation is carried in the lowest bit of the (shifted) key until encoding
            for(const auto &[i, kind] : structures)
                if(kind == Kind::canonical && 2 * (number_shapes == 0 ? window_size : shapes.shapes[i].weight) > 63)
                    throw std::runtime_error("encoded positions need k <= 31 (shapes: weight <= 31 for symmetric shapes)");
            std::cout << "encoded positions: " << ref_id_bits << " + " << ref_len_bits << " + 1 = " << width << " bits\n";

            for_each_structure([&](const size_t j) {
                const int key_shift = (structures[j].second == Kind::canonical) ? 1 : 0;
                std::vector<KeyValue> pairs;
                process_structure(j, [&](uint64_t key, uint32_t pos, bool fw = true) {
                    pairs.push_back({key_shift ? (key << 1 | fw) : key, pos});
                });
                sort_unique(pairs);
                // shapes: the forward and reverse structures
                const bool structure_fw = (structures[j].second == Kind::forward);
                maps[j] = FlatMap::from_sorted_encoded(pairs, width, key_shift, ref_len_bits, [&](const KeyValue &p) {
                    const bool fw = key_shift ? (p.key & 1) : structure_fw;
                    // sequence i starts at text position endpoints.select(i + 1) (after the padding at 0)
                    const uint64_t seq_id = endpoints.rank(uint64_t(p.value) + 1) - 2;
                    const uint64_t seq_pos = p.value - endpoints.select(seq_id + 1);
                    return (seq_id << (ref_len_bits + 1)) | (seq_pos << 1) | uint64_t(fw);
                });
            });
#else
            for_each_structure([&](const size_t j) { maps[j] = FlatMap::from_sorted(collect_sorted(j)); });
#endif
            for(size_t j = 0; j < no_structures; j++)
                (j < no_canonical ? canonical_maps : hashmaps).push_back(std::move(maps[j]));
        }
        else {
            std::vector<gtl::flat_hash_set<uint64_t>> sets(no_structures);
            for_each_structure([&](const size_t j) {
                process_structure(j, [&set = sets[j]](uint64_t key, uint32_t) { set.insert(key); });
            });
            for(size_t j = 0; j < no_structures; j++)
                (j < no_canonical ? canonical_sets : hashsets).push_back(std::move(sets[j]));
        }
    }
    else {
        ri.reserve(no_structures);
        if(loc) {
            // sized up front, the select structures point into si; each structure fills its own index
            si.resize(no_structures);
            si_select.resize(no_structures);
            offsetsi.resize(no_structures);

            std::vector<std::optional<EF>> efs(no_structures);
            for_each_structure([&](const size_t j) { efs[j] = build_level(j, collect_sorted(j)); });
            for(auto &ef : efs)
                ri.push_back(std::move(*ef));
        }
        else {
            for(int i = 0; i < no_shapes; i++) {
                gtl::flat_hash_set<uint64_t> set, set_rc;
                process(i, [&](uint64_t key, uint32_t) {set.insert(key);},
                           [&](uint64_t key, uint32_t) {set_rc.insert(key);});

                for(auto *s : {&set, &set_rc}) {
                    if(s == &set_rc && number_shapes == 0)
                        break;
                    std::vector<uint64_t> keys(s->begin(), s->end());
                    s->clear();
                    std::ranges::sort(keys);
                    ri.push_back(build_ef(keys));
                }
            }
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
    // with a threshold, the minimiser occurrences bound the k-mer/shape frequencies only for odd minimiser lengths
    // (see filter_freq_minimizers)
    if(threshold > 0) {
        const uint64_t ms[] = {m1, m2, m3};
        for(uint64_t i = 0; i < level; i++)
            if(ms[i] % 2 == 0)
                throw std::invalid_argument("the k-mer/shape threshold needs minimisers of odd length (m" + std::to_string(i + 1) + " = " + std::to_string(ms[i]) + ")");
    }
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
    // no minimiser on this level (e.g. a threshold <= the minimiser occurrences): an empty structure as for an unused
    // level, the upper bits of an empty Elias-Fano would take M1 bits
    r1 = sux::bits::EliasFano(minimizers1, minimizers1.empty() ? 1 : M1);
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
        r2 = sux::bits::EliasFano(minimizers2, minimizers2.empty() ? 1 : M2);
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
        r3 = sux::bits::EliasFano(minimizers3, minimizers3.empty() ? 1 : M3);
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

    std::cout << (store_last_level ? "build last level...\n" : "no last level, dropped...\n") << freq_kmers << " frequent kmers\n";
    last_level(input, freq_skmers);

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    print_info();
}
