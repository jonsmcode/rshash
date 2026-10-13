#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>
#include "rshash.hpp"



template <int no_shapes, bool use_ht>
void RSHash::initialise_locatefn_impl()
{
    switch (level) {
    case 1:
        streaming_locate_fn = &RSHash::streaming_locate1<no_shapes, use_ht>;
        locate_fn           = &RSHash::locate1<no_shapes, use_ht>;
        break;

    case 2:
        streaming_locate_fn = &RSHash::streaming_locate2<no_shapes, use_ht>;
        locate_fn           = &RSHash::locate2<no_shapes, use_ht>;
        break;

    case 3:
        streaming_locate_fn = &RSHash::streaming_locate3<no_shapes, use_ht>;
        locate_fn           = &RSHash::locate3<no_shapes, use_ht>;
        break;

    default:
        streaming_locate_fn = nullptr;
        locate_fn = nullptr;
    }
}

template<int no_shapes>
void RSHash::initialise_locatefn_dispatch()
{
    if (use_ht)
        initialise_locatefn_impl<no_shapes, true>();
    else
        initialise_locatefn_impl<no_shapes, false>();
}

void RSHash::initialise_locatefn()
{
    switch (number_shapes) {
    case 0:
        initialise_locatefn_dispatch<0>();
        break;
    case 1:
        initialise_locatefn_dispatch<1>();
        break;
    case 2:
        initialise_locatefn_dispatch<2>();
        break;
    case 3:
        initialise_locatefn_dispatch<3>();
        break;
    default:
        throw std::runtime_error("Unsupported number of shapes");
    }
}

uint64_t RSHash::streaming_locate(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions) {
    return (this->*streaming_locate_fn)(query, positions, no_positions);
}

uint64_t RSHash::locate(const std::vector<mask128_t> &kmers, std::vector<uint64_t> &positions, uint64_t &no_positions) {
    return (this->*locate_fn)(kmers, positions, no_positions);
}


template<int no_shapes, bool use_ht>
inline bool RSHash::locate_last_level(const uint64_t kmer, const uint64_t kmer_hi, const uint64_t kmer_rc, const uint64_t kmer_rc_hi,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    // structure j holds the forward (even j) or reverse (odd j) values of shape j/2;
    // all occurrences of all structures are counted (windows found by several shapes are counted several times)
    bool found = false;

    // a structure stores the starts of the shape windows, the (aligned) window starts dist bases earlier
    auto count = [&](const size_t j, const uint64_t key, [[maybe_unused]] const uint64_t dist) {
        if constexpr (use_ht) {
            const auto values = hashmaps[j].find(key);
            if(values.empty())
                return;
            // for(uint32_t v : values) positions.push_back(v - dist + shapes.overlap); // position of the kernel
            no_positions += values.size();
        }
        else {
            uint64_t rank;
            if(!ri[j].contains(key, rank))
                return;
            const size_t p = si_select[j].select(rank);
            // for(size_t o = p; o < si_select[j].select(rank + 1); o++) positions.push_back(offsetsi[j].access(o) - dist + shapes.overlap);
            no_positions += si_select[j].select(rank + 1) - p;
        }
        found = true;
    };

    if constexpr (no_shapes > 0) {
        for(int i = 0; i < no_shapes; ++i) {
            const Shape64 &shape = shapes.shapes[i];
            // the shape values are only extracted when the last level is reached
            const uint64_t shape_fwd = _pext_u64(kmer, shape.w_mask.lo) | (_pext_u64(kmer_hi, shape.w_mask.hi) << shape.w_lo_weight);
            const uint64_t shape_rev = _pext_u64(kmer_rc, shape.w_mask_rev.lo) | (_pext_u64(kmer_rc_hi, shape.w_mask_rev.hi) << shape.w_rev_lo_weight);
            count(2*i, shape_fwd, shape.w_dist_right);
            count(2*i + 1, shape_rev, shape.w_rev_dist_right);
        }
    }
    else
        count(0, std::min<uint64_t>(kmer, kmer_rc), 0);

    return found;
}


// k-mers given one by one (k <= 32 or shape windows), cascading through the levels like the streaming locate
template<int max_level, int no_shapes, bool use_ht>
uint64_t RSHash::locate_kmers(const std::vector<mask128_t> &kmers, std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    struct Window { uint64_t value, value_hi, value_rev, value_rev_hi; };
    uint64_t* offsets = new uint64_t[std::max({m_thres1, m_thres2, m_thres3})];
    uint64_t minimiser, minimiser_rank;
    uint64_t kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;
    uint64_t shapes_fwd[2*std::max(no_shapes, 1)], shapes_rev[2*std::max(no_shapes, 1)];
    uint64_t found_kmers = 0;

    // looks the k-mer up in the given level; false if its minimiser is not stored there
    auto locate_level = [&]<int level>(auto &r, auto &s_select) -> bool {
        minimiser = find_minimiser<level>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(!r.contains(minimiser, minimiser_rank))
            return false;
        const size_t p = s_select.select(minimiser_rank);
        const size_t no_minimiser = s_select.select(minimiser_rank+1) - p;
        found_kmers += locate_offsets<level, no_shapes>(offsets, p, no_minimiser, shapes_fwd, shapes_rev,
            left_minimiser_position, right_minimiser_position, positions, no_positions);
        return true;
    };

    for(const mask128_t &kmer : kmers) {
        const mask128_t kmer_rc = crc128(kmer, window_size);
        locate_query_words<no_shapes>(Window{kmer.lo, kmer.hi, kmer_rc.lo, kmer_rc.hi}, shapes_fwd, shapes_rev, kernel, kernel_rev);

        if(locate_level.template operator()<1>(r1, s1_select))
            continue;
        if constexpr (max_level >= 2)
            if(locate_level.template operator()<2>(r2, s2_select))
                continue;
        if constexpr (max_level >= 3)
            if(locate_level.template operator()<3>(r3, s3_select))
                continue;
        found_kmers += locate_last_level<no_shapes, use_ht>(kmer.lo, kmer.hi, kmer_rc.lo, kmer_rc.hi, positions, no_positions);
    }

    delete[] offsets;

    return found_kmers;
}

template<int no_shapes, bool use_ht>
uint64_t RSHash::locate1(const std::vector<mask128_t> &kmers, std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    return locate_kmers<1, no_shapes, use_ht>(kmers, positions, no_positions);
}

template<int no_shapes, bool use_ht>
uint64_t RSHash::locate2(const std::vector<mask128_t> &kmers, std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    return locate_kmers<2, no_shapes, use_ht>(kmers, positions, no_positions);
}

template<int no_shapes, bool use_ht>
uint64_t RSHash::locate3(const std::vector<mask128_t> &kmers, std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    return locate_kmers<3, no_shapes, use_ht>(kmers, positions, no_positions);
}



// checks one text window (lo, hi): forward against the query shapes, reverse complement against the reverse query
// shapes; in_sequence(pos, length) is only asked after the masks matched.
// returns the number of shape occurrences in the window: every matching (shape, orientation) counts, like one
// structure each on the last level; a k-mer window counts at most once, like the canonical k-mers on the last level
template<int no_shapes, typename InSequence>
inline uint64_t RSHash::match_window(const uint64_t window_lo, const uint64_t window_hi, const uint64_t window_pos,
    const bool fwd, const bool rc, const uint64_t *shapes_fwd, const uint64_t *shapes_rev, InSequence &&in_sequence)
{
    // the shape (not the whole window) has to lie within one sequence
    uint64_t matches = 0;
    if constexpr (no_shapes > 0) {
        for(int j = 0; j < no_shapes; ++j) {
            const Shape64 &shape = shapes.shapes[j];
            matches += fwd && shapes_fwd[2*j] == (window_lo & shape.w_mask.lo) && shapes_fwd[2*j + 1] == (window_hi & shape.w_mask.hi) &&
                       in_sequence(window_pos + shape.w_dist_right, shape.length);
            matches += rc && shapes_rev[2*j] == (window_lo & shape.w_mask_rev.lo) && shapes_rev[2*j + 1] == (window_hi & shape.w_mask_rev.hi) &&
                       in_sequence(window_pos + shape.w_rev_dist_right, shape.length);
        }
    }
    else {
        matches = ((fwd && shapes_fwd[0] == window_lo && shapes_fwd[1] == window_hi) ||
                   (rc && shapes_rev[0] == window_lo && shapes_rev[1] == window_hi)) &&
                  in_sequence(window_pos, window_size);
    }

    return matches;
}


// counts the shape occurrences (k-mers) matching the query among the candidate windows of all minimiser offsets;
// window_at(i, pos) gives the (lo, hi) words of the window at position pos of offset i,
// in_sequence(i, pos, length) whether pos .. pos+length-1 lies within one sequence
template<int level, int no_shapes, typename WindowAt, typename InSequence>
inline bool RSHash::locate_windows(const uint64_t *offsets, const size_t no_minimiser,
    const uint64_t *shapes_fwd, const uint64_t *shapes_rev,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    std::vector<uint64_t> &positions, uint64_t &no_positions, WindowAt &&window_at, InSequence &&in_sequence)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;

    // windows of a super-k-mer that can hold the query: forward at span-1-left / right,
    // reverse complement at left / span-1-right; the same window is only checked once
    struct Candidate { size_t pos; bool fwd, rc, fwd_right, rc_right; };
    Candidate candidates[4];
    size_t no_candidates = 0;
    auto add_candidate = [&](const size_t pos, const bool fwd, const bool right) {
        for(size_t c = 0; c < no_candidates; c++) {
            if(candidates[c].pos == pos) {
                (fwd ? candidates[c].fwd_right : candidates[c].rc_right) = right && !(fwd ? candidates[c].fwd : candidates[c].rc);
                (fwd ? candidates[c].fwd : candidates[c].rc) = true;
                return;
            }
        }
        candidates[no_candidates++] = Candidate{pos, fwd, !fwd, fwd && right, !fwd && right};
    };
    add_candidate(span-1-left_minimiser_pos, true, false);
    add_candidate(left_minimiser_pos, false, false);

    // if the minimiser occurs more than once in the query, a window has to be searched from its leftmost and its
    // rightmost minimiser occurrence; a window found from the rightmost one is only counted if its leftmost
    // occurrence (d bases to the left, forward, or to the right, reverse complement) is not an offset itself
    const size_t d = (span-1-right_minimiser_pos) - left_minimiser_pos;
    if(d > 0) {
        add_candidate(right_minimiser_pos, true, true);
        add_candidate(span-1-right_minimiser_pos, false, true);
    }
    auto is_offset = [&](const uint64_t offset) {
        return std::binary_search(offsets, offsets + no_minimiser, offset); // offsets are sorted
    };

    // todo: SIMD buffer check
    bool found = false;
    for(size_t i = 0; i < no_minimiser; i++) {
        uint64_t window_start = offsets[i];
        if constexpr (no_shapes > 0)
            window_start -= shapes.overlap;

        for(size_t c = 0; c < no_candidates; c++) {
            const Candidate &candidate = candidates[c];
            const mask128_t window = window_at(i, candidate.pos);
            const uint64_t window_pos = window_start + candidate.pos;
            auto matches = [&](const bool fwd, const bool rc) {
                return match_window<no_shapes>(window.lo, window.hi, window_pos, fwd, rc, shapes_fwd, shapes_rev,
                    [&](const uint64_t pos, const uint64_t length) { return in_sequence(i, pos, length); });
            };

            // a match from a left minimiser occurrence always counts, a match from a right one only if
            // the window's left occurrence is not an offset (only looked up after the window matched);
            // shapes count per orientation, a k-mer window only once
            const bool fwd_left = candidate.fwd && !candidate.fwd_right;
            const bool rc_left = candidate.rc && !candidate.rc_right;
            uint64_t match = (fwd_left || rc_left) ? matches(fwd_left, rc_left) : 0;
            if((no_shapes > 0 || !match) && candidate.fwd_right) {
                const uint64_t m = matches(true, false);
                if(m && !is_offset(offsets[i] - d))
                    match += m;
            }
            if((no_shapes > 0 || !match) && candidate.rc_right) {
                const uint64_t m = matches(false, true);
                if(m && !is_offset(offsets[i] + d))
                    match += m;
            }

            if(match) {
                // positions.push_back(window_pos + shapes.overlap); // position of the kernel
                no_positions += match;
                found = true;
            }
        }
    }

    return found;
}


// streaming: the windows of all offsets are in the buffer, the sequence of each minimiser is known
template<int level, int no_shapes>
inline void RSHash::locate_buffer(uint64_t *buffer, uint64_t *offsets, uint64_t *sequence_ends, const size_t no_minimiser,
    const uint64_t *shapes_fwd, const uint64_t *shapes_rev,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    std::vector<uint64_t> &positions, uint64_t &no_positions, uint64_t &found_kmers)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;

    found_kmers += locate_windows<level, no_shapes>(offsets, no_minimiser, shapes_fwd, shapes_rev,
        left_minimiser_pos, right_minimiser_pos, positions, no_positions,
        [&](const size_t i, const size_t pos) { return mask128_t{buffer[2*(i*span + pos)], buffer[2*(i*span + pos) + 1]}; },
        [&](const size_t i, const uint64_t pos, const uint64_t length) {
            return pos >= sequence_ends[2*i] && pos + length - 1 < sequence_ends[2*i+1];
        });
}


// single k-mers: only the candidate windows are read from the text, the sequence is only looked up for matches
template<int level, int no_shapes>
inline bool RSHash::locate_offsets(uint64_t *offsets, const size_t p, const size_t no_minimiser,
    const uint64_t *shapes_fwd, const uint64_t *shapes_rev,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;

    for(size_t i = 0; i < no_minimiser; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i) + 1 - span;
        if constexpr (level == 2)
            offsets[i] = offsets2.access(p+i) + 1 - span;
        if constexpr (level == 3)
            offsets[i] = offsets3.access(p+i) + 1 - span;
    }

    return locate_windows<level, no_shapes>(offsets, no_minimiser, shapes_fwd, shapes_rev,
        left_minimiser_pos, right_minimiser_pos, positions, no_positions,
        [&](const size_t i, const size_t pos) {
            uint64_t window_start = offsets[i] + pos;
            if constexpr (no_shapes > 0)
                window_start -= shapes.overlap;
            return access(window_start);
        },
        [&](const size_t, const uint64_t pos, const uint64_t length) {
            uint64_t start_pos, end_pos;
            return check_overlap(pos, start_pos, end_pos, length);
        });
}


// query words of a window: the masked (lo, hi) words per shape for the buffer checks and the kernel for the minimisers
template<int no_shapes, typename Window>
inline void RSHash::locate_query_words(const Window &window, uint64_t *shapes_fwd, uint64_t *shapes_rev,
    uint64_t &kernel, uint64_t &kernel_rev)
{
    if constexpr (no_shapes > 0) {
        for(int i = 0; i < no_shapes; ++i) {
            const Shape64 &shape = shapes.shapes[i];
            shapes_fwd[2*i] =     window.value & shape.w_mask.lo;
            shapes_fwd[2*i + 1] = window.value_hi & shape.w_mask.hi;
            shapes_rev[2*i] =     window.value_rev & shape.w_mask_rev.lo;
            shapes_rev[2*i + 1] = window.value_rev_hi & shape.w_mask_rev.hi;
        }
        kernel = (window.value & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
        kernel_rev = (window.value_rev & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_rev_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
    }
    else {
        shapes_fwd[0] = window.value;
        shapes_fwd[1] = window.value_hi;
        shapes_rev[0] = window.value_rev;
        shapes_rev[1] = window.value_rev_hi;
        kernel = window.value;
        kernel_rev = window.value_rev;
    }
}


template<int no_shapes, bool use_ht>
uint64_t RSHash::streaming_locate1(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t found_kmers = 0;
    uint64_t current_pos_minimiser1 = INF, current_neg_minimiser1 = INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* kmer_buffer1 = new uint64_t[2*(m_thres1 * span1)];
    uint64_t* sequence_ends1 = new uint64_t[2*m_thres1];
    size_t no_minimiser1;
    uint64_t minimiser1, minimiser_rank1;
    size_t left_minimiser_position1, right_minimiser_position1;
    bool rolling1 = false;
    uint64_t kernel, kernel_rev;
    uint64_t shapes_fwd[2*std::max(no_shapes, 1)], shapes_rev[2*std::max(no_shapes, 1)];

    for (auto&& window : query | rshash::views::kmer_view({.window_size = window_size}))
    {
        locate_query_words<no_shapes>(window, shapes_fwd, shapes_rev, kernel, kernel_rev);

        if(!rolling1) {
            minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position1, right_minimiser_position1);
            rolling1 = true;
        }
        else
            update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser_position1, right_minimiser_position1);

        if(minimiser1 == current_pos_minimiser1) {
            locate_buffer<1, no_shapes>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, shapes_fwd, shapes_rev, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
        }
        else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser_rank1)) {
            const size_t minimiser_position1 = s1_select.select(minimiser_rank1);
            no_minimiser1 = s1_select.select(minimiser_rank1+1) - minimiser_position1;

            fill_buffer2<1>(offsets1, kmer_buffer1, sequence_ends1, minimiser_position1, no_minimiser1);
            locate_buffer<1, no_shapes>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, shapes_fwd, shapes_rev, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
            current_pos_minimiser1 = minimiser1;
        }
        else {
            found_kmers += locate_last_level<no_shapes, use_ht>(window.value, window.value_hi, window.value_rev, window.value_rev_hi, positions, no_positions);
            current_neg_minimiser1 = minimiser1;
        }
    }

    delete[] kmer_buffer1;
    delete[] offsets1;
    delete[] sequence_ends1;

    return found_kmers;
}


template<int no_shapes, bool use_ht>
uint64_t RSHash::streaming_locate2(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t found_kmers = 0;
    uint64_t current_pos_minimiser1 = INF, current_neg_minimiser1 = INF;
    uint64_t current_pos_minimiser2 = INF, current_neg_minimiser2 = INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* kmer_buffer1 = new uint64_t[2*(m_thres1 * span1)];
    uint64_t* kmer_buffer2 = new uint64_t[2*(m_thres2 * span2)];
    uint64_t* sequence_ends1 = new uint64_t[2*m_thres1];
    uint64_t* sequence_ends2 = new uint64_t[2*m_thres2];
    size_t no_minimiser1, no_minimiser2;
    uint64_t minimiser1, minimiser_rank1;
    uint64_t minimiser2, minimiser_rank2;
    size_t left_minimiser_position1, right_minimiser_position1;
    size_t left_minimiser_position2, right_minimiser_position2;
    bool rolling1 = false;
    bool rolling2 = false;
    uint64_t kernel, kernel_rev;
    uint64_t shapes_fwd[2*std::max(no_shapes, 1)], shapes_rev[2*std::max(no_shapes, 1)];

    for (auto&& window : query | rshash::views::kmer_view({.window_size = window_size}))
    {
        locate_query_words<no_shapes>(window, shapes_fwd, shapes_rev, kernel, kernel_rev);

        if(!rolling1) {
            minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position1, right_minimiser_position1);
            rolling1 = true;
        }
        else
            update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser_position1, right_minimiser_position1);

        if(minimiser1 == current_pos_minimiser1) {
            locate_buffer<1, no_shapes>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, shapes_fwd, shapes_rev, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
            rolling2 = false;
        }
        else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser_rank1)) {
            const size_t minimiser_position1 = s1_select.select(minimiser_rank1);
            no_minimiser1 = s1_select.select(minimiser_rank1+1) - minimiser_position1;

            fill_buffer2<1>(offsets1, kmer_buffer1, sequence_ends1, minimiser_position1, no_minimiser1);
            locate_buffer<1, no_shapes>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, shapes_fwd, shapes_rev, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
            current_pos_minimiser1 = minimiser1;
            rolling2 = false;
        }
        else {
            if(!rolling2) {
                minimiser2 = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position2, right_minimiser_position2);
                rolling2 = true;
            }
            else
                update_minimiser<2>(kernel, kernel_rev, minimiser2, left_minimiser_position2, right_minimiser_position2);

            if(minimiser2 == current_pos_minimiser2) {
                locate_buffer<2, no_shapes>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, shapes_fwd, shapes_rev, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
            }
            else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser_rank2)) {
                const size_t minimiser_position2 = s2_select.select(minimiser_rank2);
                no_minimiser2 = s2_select.select(minimiser_rank2+1) - minimiser_position2;

                fill_buffer2<2>(offsets2, kmer_buffer2, sequence_ends2, minimiser_position2, no_minimiser2);
                locate_buffer<2, no_shapes>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, shapes_fwd, shapes_rev, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
                current_pos_minimiser2 = minimiser2;
                current_neg_minimiser1 = minimiser1;
            }
            else {
                found_kmers += locate_last_level<no_shapes, use_ht>(window.value, window.value_hi, window.value_rev, window.value_rev_hi, positions, no_positions);
                current_neg_minimiser1 = minimiser1;
                current_neg_minimiser2 = minimiser2;
            }
        }
    }

    delete[] kmer_buffer1;
    delete[] kmer_buffer2;
    delete[] offsets1;
    delete[] offsets2;
    delete[] sequence_ends1;
    delete[] sequence_ends2;

    return found_kmers;
}


template<int no_shapes, bool use_ht>
uint64_t RSHash::streaming_locate3(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t found_kmers = 0;
    uint64_t current_pos_minimiser1 = INF, current_neg_minimiser1 = INF;
    uint64_t current_pos_minimiser2 = INF, current_neg_minimiser2 = INF;
    uint64_t current_pos_minimiser3 = INF, current_neg_minimiser3 = INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* offsets3 = new uint64_t[m_thres3];
    uint64_t* kmer_buffer1 = new uint64_t[2*(m_thres1 * span1)];
    uint64_t* kmer_buffer2 = new uint64_t[2*(m_thres2 * span2)];
    uint64_t* kmer_buffer3 = new uint64_t[2*(m_thres3 * span3)];
    uint64_t* sequence_ends1 = new uint64_t[2*m_thres1];
    uint64_t* sequence_ends2 = new uint64_t[2*m_thres2];
    uint64_t* sequence_ends3 = new uint64_t[2*m_thres3];
    size_t no_minimiser1, no_minimiser2, no_minimiser3;
    uint64_t minimiser1, minimiser_rank1;
    uint64_t minimiser2, minimiser_rank2;
    uint64_t minimiser3, minimiser_rank3;
    size_t left_minimiser_position1, right_minimiser_position1;
    size_t left_minimiser_position2, right_minimiser_position2;
    size_t left_minimiser_position3, right_minimiser_position3;
    bool rolling1 = false;
    bool rolling2 = false;
    bool rolling3 = false;
    uint64_t kernel, kernel_rev;
    uint64_t shapes_fwd[2*std::max(no_shapes, 1)], shapes_rev[2*std::max(no_shapes, 1)];

    for (auto&& window : query | rshash::views::kmer_view({.window_size = window_size}))
    {
        locate_query_words<no_shapes>(window, shapes_fwd, shapes_rev, kernel, kernel_rev);

        if(!rolling1) {
            minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position1, right_minimiser_position1);
            rolling1 = true;
        }
        else
            update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser_position1, right_minimiser_position1);

        if(minimiser1 == current_pos_minimiser1) {
            locate_buffer<1, no_shapes>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, shapes_fwd, shapes_rev, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
            rolling2 = false;
            rolling3 = false;
        }
        else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser_rank1)) {
            const size_t minimiser_position1 = s1_select.select(minimiser_rank1);
            no_minimiser1 = s1_select.select(minimiser_rank1+1) - minimiser_position1;

            fill_buffer2<1>(offsets1, kmer_buffer1, sequence_ends1, minimiser_position1, no_minimiser1);
            locate_buffer<1, no_shapes>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, shapes_fwd, shapes_rev, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
            current_pos_minimiser1 = minimiser1;
            rolling2 = false;
            rolling3 = false;
        }
        else {
            if(!rolling2) {
                minimiser2 = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position2, right_minimiser_position2);
                rolling2 = true;
            }
            else
                update_minimiser<2>(kernel, kernel_rev, minimiser2, left_minimiser_position2, right_minimiser_position2);

            if(minimiser2 == current_pos_minimiser2) {
                locate_buffer<2, no_shapes>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, shapes_fwd, shapes_rev, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
                rolling3 = false;
            }
            else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser_rank2)) {
                const size_t minimiser_position2 = s2_select.select(minimiser_rank2);
                no_minimiser2 = s2_select.select(minimiser_rank2+1) - minimiser_position2;

                fill_buffer2<2>(offsets2, kmer_buffer2, sequence_ends2, minimiser_position2, no_minimiser2);
                locate_buffer<2, no_shapes>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, shapes_fwd, shapes_rev, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
                current_pos_minimiser2 = minimiser2;
                current_neg_minimiser1 = minimiser1;
                rolling3 = false;
            }
            else {
                if(!rolling3) {
                    minimiser3 = find_minimiser<3>(kernel, kernel_rev, left_minimiser_position3, right_minimiser_position3);
                    rolling3 = true;
                }
                else
                    update_minimiser<3>(kernel, kernel_rev, minimiser3, left_minimiser_position3, right_minimiser_position3);

                if(minimiser3 == current_pos_minimiser3) {
                    locate_buffer<3, no_shapes>(kmer_buffer3, offsets3, sequence_ends3, no_minimiser3, shapes_fwd, shapes_rev, left_minimiser_position3, right_minimiser_position3, positions, no_positions, found_kmers);
                }
                else if(minimiser3 != current_neg_minimiser3 && r3.contains(minimiser3, minimiser_rank3)) {
                    const size_t minimiser_position3 = s3_select.select(minimiser_rank3);
                    no_minimiser3 = s3_select.select(minimiser_rank3+1) - minimiser_position3;

                    fill_buffer2<3>(offsets3, kmer_buffer3, sequence_ends3, minimiser_position3, no_minimiser3);
                    locate_buffer<3, no_shapes>(kmer_buffer3, offsets3, sequence_ends3, no_minimiser3, shapes_fwd, shapes_rev, left_minimiser_position3, right_minimiser_position3, positions, no_positions, found_kmers);
                    current_pos_minimiser3 = minimiser3;
                    current_neg_minimiser1 = minimiser1;
                    current_neg_minimiser2 = minimiser2;
                }
                else {
                    found_kmers += locate_last_level<no_shapes, use_ht>(window.value, window.value_hi, window.value_rev, window.value_rev_hi, positions, no_positions);
                    current_neg_minimiser1 = minimiser1;
                    current_neg_minimiser2 = minimiser2;
                    current_neg_minimiser3 = minimiser3;
                }
            }
        }
    }

    delete[] kmer_buffer1;
    delete[] kmer_buffer2;
    delete[] kmer_buffer3;
    delete[] offsets1;
    delete[] offsets2;
    delete[] offsets3;
    delete[] sequence_ends1;
    delete[] sequence_ends2;
    delete[] sequence_ends3;

    return found_kmers;
}
