#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>
#include "rshash.hpp"


template <int no_shapes, bool use_ht, bool locate>
void RSHash::initialise_lookupfn_impl()
{
    switch (level) {
    case 1:
        streaming_lookup_fn = &RSHash::streaming_lookup1<no_shapes, use_ht, locate>;
        lookup_fn           = &RSHash::lookup1<no_shapes, use_ht, locate>;
        break;

    case 2:
        streaming_lookup_fn = &RSHash::streaming_lookup2<no_shapes, use_ht, locate>;
        lookup_fn           = &RSHash::lookup2<no_shapes, use_ht, locate>;
        break;

    case 3:
        streaming_lookup_fn = &RSHash::streaming_lookup3<no_shapes, use_ht, locate>;
        lookup_fn           = &RSHash::lookup3<no_shapes, use_ht, locate>;
        break;

    default:
        streaming_lookup_fn = nullptr;
        lookup_fn = nullptr;
    }
}


template<int no_shapes>
void RSHash::initialise_lookupfn_dispatch()
{
    if (use_ht) {
        if (loc)
            initialise_lookupfn_impl<no_shapes, true, true>();
        else
            initialise_lookupfn_impl<no_shapes, true, false>();
    } else {
        if (loc)
            initialise_lookupfn_impl<no_shapes, false, true>();
        else
            initialise_lookupfn_impl<no_shapes, false, false>();
    }
}

void RSHash::initialise_lookupfn()
{
    switch (number_shapes) {
    case 0:
        initialise_lookupfn_dispatch<0>();
        break;
    case 1:
        initialise_lookupfn_dispatch<1>();
        break;
    case 2:
        initialise_lookupfn_dispatch<2>();
        break;
    case 3:
        initialise_lookupfn_dispatch<3>();
        break;
    default:
        throw std::runtime_error("Unsupported number of shapes");
    }
}


uint64_t RSHash::streaming_lookup(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions) {
    return (this->*streaming_lookup_fn)(query, extensions);
}

uint64_t RSHash::lookup(const std::vector<mask128_t> &kmers) {
    return (this->*lookup_fn)(kmers);
}


template <int no_shapes, typename Forward>
inline bool contains_impl(Forward &forward, const std::vector<Shape64>& shapes,
    const uint64_t kmer, const uint64_t kmer_hi, const uint64_t kmer_rc, const uint64_t kmer_rc_hi)
{
    if constexpr (no_shapes > 0) {
        for(int i = 0; i < no_shapes; ++i) {
            Shape64 shape = shapes[i];
            const uint64_t shape_fwd = _pext_u64(kmer, shape.w_mask.lo) | (_pext_u64(kmer_hi, shape.w_mask.hi) << shape.w_lo_weight);
            const uint64_t shape_rev  = _pext_u64(kmer_rc, shape.w_mask_rev.lo) | (_pext_u64(kmer_rc_hi, shape.w_mask_rev.hi) << shape.w_rev_lo_weight);

            if(forward[2*i].contains(shape_fwd) || forward[2*i+1].contains(shape_rev))
                return true;
        }

        return false;
    }
    else
        return forward[0].contains(std::min<uint64_t>(kmer, kmer_rc));
}

template<int no_shapes, bool use_ht, bool locate>
inline bool RSHash::lookup_last_level(const uint64_t kmer, const uint64_t kmer_hi, const uint64_t kmer_rc, const uint64_t kmer_rc_hi)
{
    if constexpr (use_ht) {
        if constexpr (locate)
            return contains_impl<no_shapes>(hashmaps, shapes.shapes, kmer, kmer_hi, kmer_rc, kmer_rc_hi);
        else
            return contains_impl<no_shapes>(hashsets, shapes.shapes, kmer, kmer_hi, kmer_rc, kmer_rc_hi);
    }
    return false;
    // todo: multiple Rs for multiple shapes!
    // else {
    //     return contains_impl<no_shapes>(r4, r5, kmer, kmer_rc, shapes_fwd, shapes_rev);
    // }
}


template<int no_shapes, bool use_ht, bool locate>
uint64_t RSHash::lookup1(const std::vector<mask128_t> &kmers)
{
    uint64_t* offsets = new uint64_t[m_thres1];
    uint64_t minimiser, minimiser_rank, kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;
    constexpr size_t n = std::max(2, 2*no_shapes);
    uint64_t shapes_fwd[n], shapes_rev[n];
    uint64_t occurences = 0;
    
    for(mask128_t kmer : kmers) {
        mask128_t kmer_rc = crc128(kmer, window_size);

        if constexpr (no_shapes > 0) {
            kernel = (kmer.lo & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((kmer.hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
            kernel_rev = (kmer_rc.lo & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((kmer_rc.hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);

            for(int i = 0; i < no_shapes; ++i) {
                Shape64 &shape = shapes.shapes[i];
                shapes_fwd[2*i] = kmer.lo & shape.w_mask.lo;
                shapes_fwd[2*i+1] = kmer.hi & shape.w_mask.hi;
                shapes_rev[2*i] = kmer_rc.lo & shape.w_mask_rev.lo;
                shapes_rev[2*i+1] = kmer_rc.hi & shape.w_mask_rev.hi;
            }
        }
        else {
            kernel = kmer.lo;
            kernel_rev = kmer_rc.lo;
            shapes_fwd[0] = kmer.lo;
            shapes_fwd[1] = kmer.hi;
            shapes_rev[0] = kmer_rc.lo;
            shapes_rev[1] = kmer_rc.hi;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1, no_shapes>(shapes_fwd, shapes_rev, offsets, p, left_minimiser_position, right_minimiser_position, no_minimiser);
        }
        else {
            occurences += lookup_last_level<no_shapes, use_ht, locate>(kmer.lo, kmer.hi, kmer_rc.lo, kmer_rc.hi);
        }
    }

    delete[] offsets;

    return occurences;
}


template<int no_shapes, bool use_ht, bool locate>
uint64_t RSHash::lookup2(const std::vector<mask128_t> &kmers)
{
    uint64_t* offsets = new uint64_t[std::max(m_thres1, m_thres2)];
    uint64_t minimiser, minimiser_rank, kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;
    constexpr size_t n = std::max(2, 2*no_shapes);
    uint64_t shapes_fwd[n], shapes_rev[n];
    uint64_t occurences = 0;

    for(mask128_t kmer : kmers) {
        mask128_t kmer_rc = crc128(kmer, window_size);

        if constexpr (no_shapes > 0) {
            kernel = (kmer.lo & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((kmer.hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
            kernel_rev = (kmer_rc.lo & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((kmer_rc.hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);

            for(int i = 0; i < no_shapes; ++i) {
                Shape64 &shape = shapes.shapes[i];
                shapes_fwd[2*i] = kmer.lo & shape.w_mask.lo;
                shapes_fwd[2*i+1] = kmer.hi & shape.w_mask.hi;
                shapes_rev[2*i] = kmer_rc.lo & shape.w_mask_rev.lo;
                shapes_rev[2*i+1] = kmer_rc.hi & shape.w_mask_rev.hi;
            }
        }
        else {
            kernel = kmer.lo;
            kernel_rev = kmer_rc.lo;
            shapes_fwd[0] = kmer.lo;
            shapes_fwd[1] = kmer.hi;
            shapes_rev[0] = kmer_rc.lo;
            shapes_rev[1] = kmer_rc.hi;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1, no_shapes>(shapes_fwd, shapes_rev, offsets, p, left_minimiser_position, right_minimiser_position, no_minimiser);
            continue;
        }

        minimiser = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r2.contains(minimiser, minimiser_rank)) {
            size_t p = s2_select.select(minimiser_rank);
            size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
            occurences += check<2, no_shapes>(shapes_fwd, shapes_rev, offsets, p, left_minimiser_position, right_minimiser_position, no_minimiser);
        }
        else {
            occurences += lookup_last_level<no_shapes, use_ht, locate>(kmer.lo, kmer.hi, kmer_rc.lo, kmer_rc.hi);
        }
    }

    delete[] offsets;

    return occurences;
}


template<int no_shapes, bool use_ht, bool locate>
uint64_t RSHash::lookup3(const std::vector<mask128_t> &kmers)
{
    uint64_t* offsets = new uint64_t[std::max({m_thres1, m_thres2, m_thres3})];
    uint64_t minimiser, minimiser_rank, kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;
    constexpr size_t n = std::max(2, 2*no_shapes);
    uint64_t shapes_fwd[n], shapes_rev[n];
    uint64_t occurences = 0;

    for(mask128_t kmer : kmers) {
        mask128_t kmer_rc = crc128(kmer, window_size);

        if constexpr (no_shapes > 0) {
            kernel = (kmer.lo & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((kmer.hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
            kernel_rev = (kmer_rc.lo & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((kmer_rc.hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);

            for(int i = 0; i < no_shapes; ++i) {
                Shape64 &shape = shapes.shapes[i];
                shapes_fwd[2*i] = kmer.lo & shape.w_mask.lo;
                shapes_fwd[2*i+1] = kmer.hi & shape.w_mask.hi;
                shapes_rev[2*i] = kmer_rc.lo & shape.w_mask_rev.lo;
                shapes_rev[2*i+1] = kmer_rc.hi & shape.w_mask_rev.hi;
            }
        }
        else {
            kernel = kmer.lo;
            kernel_rev = kmer_rc.lo;
            shapes_fwd[0] = kmer.lo;
            shapes_fwd[1] = kmer.hi;
            shapes_rev[0] = kmer_rc.lo;
            shapes_rev[1] = kmer_rc.hi;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1, no_shapes>(shapes_fwd, shapes_rev, offsets, p, left_minimiser_position, right_minimiser_position, no_minimiser);
            continue;
        }

        minimiser = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r2.contains(minimiser, minimiser_rank)) {
            size_t p = s2_select.select(minimiser_rank);
            size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
            occurences += check<2, no_shapes>(shapes_fwd, shapes_rev, offsets, p, left_minimiser_position, right_minimiser_position, no_minimiser);
            continue;
        }

        minimiser = find_minimiser<3>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r3.contains(minimiser, minimiser_rank)) {
            size_t p = s3_select.select(minimiser_rank);
            size_t no_minimiser = s3_select.select(minimiser_rank+1) - p;
            occurences += check<3, no_shapes>(shapes_fwd, shapes_rev, offsets, p, left_minimiser_position, right_minimiser_position, no_minimiser);
        }
        else {
            occurences += lookup_last_level<no_shapes, use_ht, locate>(kmer.lo, kmer.hi, kmer_rc.lo, kmer_rc.hi);
        }
    }

    delete[] offsets;

    return occurences;
}


template<int level, int no_shapes>
inline bool RSHash::check(const uint64_t* shapes_fwd, const uint64_t* shapes_rev, uint64_t* offsets,
    const size_t p, const size_t left_minimiser_position, const size_t right_minimiser_position, const size_t no_minimiser)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;

    auto matches = [&](uint64_t pos, uint64_t pos_rc) -> bool {
        if constexpr (no_shapes > 0) {
            pos -= shapes.overlap;
            pos_rc -= shapes.overlap;
        }
        const uint64_t word_fwd_lo = get_word64(pos);
        const uint64_t word_fwd_hi = get_word64(pos + 32);
        const uint64_t word_rc_lo  = get_word64(pos_rc);
        const uint64_t word_rc_hi  = get_word64(pos_rc + 32);
        // one cache line?!

        if constexpr (no_shapes > 0) {
            for(int j = 0; j < no_shapes; ++j) {
                const Shape64& shape = shapes.shapes[j];

                if((shapes_fwd[2 * j] == (word_fwd_lo & shape.w_mask.lo) && shapes_fwd[2 * j + 1] == (word_fwd_hi & shape.w_mask.hi)) ||
                   (shapes_rev[2 * j] == (word_rc_lo & shape.w_mask_rev.lo) && shapes_rev[2 * j + 1] == (word_rc_hi & shape.w_mask_rev.hi)))
                    return true;
            }
        }
        else {
            if((shapes_fwd[0] == word_fwd_lo && shapes_fwd[1] == (word_fwd_hi & windowmask.hi)) || 
               (shapes_rev[0] == word_rc_lo && shapes_rev[1] == (word_rc_hi & windowmask.hi)))
                return true;
        }

        return false;
    };

    for(size_t i = 0; i < no_minimiser; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i)-span+1;
        if constexpr (level == 2)
            offsets[i] = offsets2.access(p+i)-span+1;
        if constexpr (level == 3)
            offsets[i] = offsets3.access(p+i)-span+1;
    }

    if (left_minimiser_position == span-1-right_minimiser_position) {
        for (size_t i = 0; i < no_minimiser; ++i) {
            const uint64_t offset = offsets[i];

            const uint64_t pos = offset + span-1-left_minimiser_position;
            const uint64_t pos_rc = offset + left_minimiser_position;

            if (matches(pos, pos_rc))
                return true;
        }
    }
    else {
        for (size_t i = 0; i < no_minimiser; ++i) {
            const uint64_t offset = offsets[i];
            
            const uint64_t pos1 = offset + span-1-left_minimiser_position;
            const uint64_t pos_rc1 = offset + left_minimiser_position;
            const uint64_t pos2 = offset + right_minimiser_position;
            const uint64_t pos_rc2 = offset + span-1-right_minimiser_position;

            if (matches(pos1, pos_rc1) || matches(pos2, pos_rc2))
                return true;
        }
    }

    return false;
}



template<int no_shapes>
inline bool RSHash::extend_in_text(uint64_t &text_pos, uint64_t start, uint64_t end,
    bool forward, const uint64_t query, const uint64_t query_rc, uint64_t* shapes_fwd, uint64_t* shapes_rev,
    uint64_t &window_lo, uint64_t &window_hi)
{
    if constexpr (no_shapes == 0) {
        if(forward) {
            if(++text_pos < end) {
                const uint64_t new_rank = get_base(text_pos);
                return new_rank == (query >> windowshift);
            }
        }
        else {
            if(--text_pos >= start) {
                const uint64_t new_rank = get_base(text_pos);
                return new_rank == (query_rc & 0b11);
            }
        }
        return false;
    }
    else {
        if(forward) {
            if(++text_pos < end - shapes.w_dist) {
                const uint64_t new_rank = get_base(text_pos);
                window_lo = (window_lo >> 2) | ((window_hi & 3ULL) << 62);
                window_hi = (window_hi >> 2) | (new_rank << windowshift);

                for(int i = 0; i < no_shapes; i++) {
                    Shape64 &shape = shapes.shapes[i];
                    if(shapes_fwd[2*i] == (window_lo & shape.w_mask.lo) && shapes_fwd[2*i + 1] == (window_hi & shape.w_mask.hi))
                        return true;
                }
            }
        }
        else {
            if(--text_pos >= start + shapes.w_dist) {
                const uint64_t new_rank = get_base(text_pos);
                window_hi = ((window_hi << 2) | (window_lo >> 62)) & windowmask.hi;
                window_lo = (window_lo << 2) | new_rank;
            
                for(int i = 0; i < no_shapes; i++) {
                    Shape64 &shape = shapes.shapes[i];
                    if(shapes_rev[2*i] == (window_lo & shape.w_mask_rev.lo) && shapes_rev[2*i + 1] == (window_hi & shape.w_mask_rev.hi))
                        return true;
                }
            }
        }
        return false;
    }
}


template<int level, int no_shapes>
inline bool RSHash::check_shape(uint64_t *buffer, uint64_t offset,
    const uint64_t kmer, const uint64_t kmer_rc, uint64_t* shapes_fwd, uint64_t* shapes_rev,
    const size_t s, const size_t minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos,
    uint64_t &text_kmer, uint64_t &text_kmer_hi)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;

    auto check_candidate_fwd = [&](uint64_t candidate_lo, uint64_t candidate_hi, uint64_t query_lo, uint64_t query_hi,
        uint64_t window, uint64_t window_hi,
        uint64_t window_pos, uint64_t shape_pos, uint64_t shape_length) -> bool
    {
        if (candidate_lo == query_lo && candidate_hi == query_hi && check_overlap(shape_pos, start_pos, end_pos, shape_length)) {
            forward = true;
            text_pos = window_pos + window_size - 1;
            text_kmer = window;
            text_kmer_hi = window_hi;
            return true;
        }
        return false;
    };
    auto check_candidate_rev = [&](uint64_t candidate_lo, uint64_t candidate_hi, uint64_t query_lo, uint64_t query_hi,
        uint64_t window, uint64_t window_hi,
        uint64_t window_pos, uint64_t shape_pos, uint64_t shape_length) -> bool
    {
        if (candidate_lo == query_lo && candidate_hi == query_hi && check_overlap(shape_pos, start_pos, end_pos, shape_length)) {
            forward = false;
            text_pos = window_pos;
            text_kmer = window;
            text_kmer_hi = window_hi;
            return true;
        }
        return false;
    };

    const uint64_t pos = span-1-minimiser_pos;
    const uint64_t kmer_pos = offset + pos;
    const uint64_t pos_rc = minimiser_pos;
    const uint64_t kmer_pos_rc = offset + pos_rc;
    const uint64_t window_lo = buffer[2*(s + pos)];
    const uint64_t window_hi = buffer[2*(s + pos) + 1];
    const uint64_t window_rc_lo = buffer[2*(s + pos_rc)];
    const uint64_t window_rc_hi = buffer[2*(s + pos_rc) + 1];

    if constexpr (no_shapes > 0) {
        const uint64_t window_pos = kmer_pos - shapes.overlap;
        const uint64_t window_pos_rc = kmer_pos_rc - shapes.overlap;

        for(int i = 0; i < no_shapes; ++i) {
            const Shape64 &shape = shapes.shapes[i];
            uint64_t candidate_lo = window_lo & shape.w_mask.lo;
            uint64_t candidate_hi = window_hi & shape.w_mask.hi;
            if(check_candidate_fwd(candidate_lo, candidate_hi, shapes_fwd[2*i], shapes_fwd[2*i + 1], window_lo, window_hi, window_pos, window_pos + shape.w_dist_right, shape.length))
                return true;

            uint64_t candidate_rc_lo = window_rc_lo & shape.w_mask_rev.lo;
            uint64_t candidate_rc_hi = window_rc_hi & shape.w_mask_rev.hi;
            if(check_candidate_rev(candidate_rc_lo, candidate_rc_hi, shapes_rev[2*i], shapes_rev[2*i + 1], window_rc_lo, window_rc_hi, window_pos_rc, window_pos_rc + shape.w_rev_dist_right, shape.length))
                return true;
        }
    }
    else {
        return check_candidate_fwd(window_lo, window_hi, kmer, kmer, kmer, kmer, kmer_pos, kmer_pos, window_size)
            || check_candidate_rev(window_rc_lo, window_rc_hi, kmer_rc, kmer_rc, kmer_rc, kmer_rc, kmer_pos_rc, kmer_pos_rc, window_size);
    }

    return false;
}


template<int level, int no_shapes>
inline bool RSHash::check_shape2(uint64_t *buffer, uint64_t offset,
    const uint64_t kmer, const uint64_t kmer_rc, uint64_t* shapes_fwd, uint64_t* shapes_rev,
    const size_t s, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos,
    uint64_t &text_kmer_lo, uint64_t &text_kmer_hi)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;

    auto check_candidate_fwd = [&](uint64_t candidate_lo, uint64_t candidate_hi, uint64_t query_lo, uint64_t query_hi,
        uint64_t window_lo, uint64_t window_hi,
        uint64_t window_pos, uint64_t shape_pos, uint64_t shape_length) -> bool
    {
        if (candidate_lo == query_lo && candidate_hi == query_hi && check_overlap(shape_pos, start_pos, end_pos, shape_length)) {
            forward = true;
            text_pos = window_pos + window_size - 1;
            text_kmer_lo = window_lo;
            text_kmer_hi = window_hi;
            return true;
        }
        return false;
    };
    auto check_candidate_rev = [&](uint64_t candidate_lo, uint64_t candidate_hi, uint64_t query_lo, uint64_t query_hi,
        uint64_t window_lo, uint64_t window_hi,
        uint64_t window_pos, uint64_t shape_pos, uint64_t shape_length) -> bool
    {
        if (candidate_lo == query_lo && candidate_hi == query_hi && check_overlap(shape_pos, start_pos, end_pos, shape_length)) {
            forward = false;
            text_pos = window_pos;
            text_kmer_lo = window_lo;
            text_kmer_hi = window_hi;
            return true;
        }
        return false;
    };

    uint64_t left_pos = span-1-left_minimiser_pos;
    uint64_t left_kmer_pos = offset + left_pos;
    uint64_t left_pos_rc = left_minimiser_pos;
    uint64_t left_kmer_pos_rc = offset + left_pos_rc;
    uint64_t right_pos = right_minimiser_pos;
    uint64_t right_kmer_pos = offset + right_pos;
    uint64_t right_pos_rc = span-1-right_minimiser_pos;
    uint64_t right_kmer_pos_rc = offset + right_pos_rc;
    uint64_t left_window_lo = buffer[2*(s + left_pos)];
    uint64_t left_window_hi = buffer[2*(s + left_pos) + 1];
    uint64_t left_window_rev_lo = buffer[2*(s + left_pos_rc)];
    uint64_t left_window_rev_hi = buffer[2*(s + left_pos_rc) + 1];
    uint64_t right_window_lo = buffer[2*(s + right_pos)];
    uint64_t right_window_hi = buffer[2*(s + right_pos) + 1];
    uint64_t right_window_rev_lo = buffer[2*(s + right_pos_rc)];
    uint64_t right_window_rev_hi = buffer[2*(s + right_pos_rc) + 1];

    if constexpr (no_shapes > 0) {
        const uint64_t left_window_pos = left_kmer_pos - shapes.overlap;
        const uint64_t left_window_pos_rc = left_kmer_pos_rc - shapes.overlap;
        const uint64_t right_window_pos = right_kmer_pos - shapes.overlap;
        const uint64_t right_window_pos_rc = right_kmer_pos_rc - shapes.overlap;
        for(int i = 0; i < no_shapes; ++i) {
            const Shape64 &shape = shapes.shapes[i];

            uint64_t left_candidate_lo = left_window_lo & shape.w_mask.lo;
            uint64_t left_candidate_hi = left_window_hi & shape.w_mask.hi;
            if(check_candidate_fwd(left_candidate_lo, left_candidate_hi, shapes_fwd[2*i], shapes_fwd[2*i + 1], left_window_lo, left_window_hi, left_window_pos, left_window_pos + shape.w_dist_right, shape.length))
                return true;

            uint64_t left_candidate_rc_lo = left_window_rev_lo & shape.w_mask_rev.lo;
            uint64_t left_candidate_rc_hi = left_window_rev_hi & shape.w_mask_rev.hi;
            if(check_candidate_rev(left_candidate_rc_lo, left_candidate_rc_hi, shapes_rev[2*i], shapes_rev[2*i + 1], left_window_rev_lo, left_window_rev_hi, left_window_pos_rc, left_window_pos_rc + shape.w_rev_dist_right, shape.length))
                return true;

            uint64_t right_candidate_lo = right_window_lo & shape.w_mask.lo;
            uint64_t right_candidate_hi = right_window_hi & shape.w_mask.hi;
            if(check_candidate_fwd(right_candidate_lo, right_candidate_hi, shapes_fwd[2*i], shapes_fwd[2*i + 1], right_window_lo, right_window_hi, right_window_pos, right_window_pos + shape.w_dist_right, shape.length))
                return true;

            uint64_t right_candidate_rc_lo = right_window_rev_lo & shape.w_mask_rev.lo;
            uint64_t right_candidate_rc_hi = right_window_rev_hi & shape.w_mask_rev.hi;
            if(check_candidate_rev(right_candidate_rc_lo, right_candidate_rc_hi, shapes_rev[2*i], shapes_rev[2*i + 1], right_window_rev_lo, right_window_rev_hi, right_window_pos_rc, right_window_pos_rc + shape.w_rev_dist_right, shape.length))
                return true;
        }
    }
    else {
        return check_candidate_fwd(left_window_lo, left_window_hi, kmer, kmer, kmer, kmer, left_kmer_pos, left_kmer_pos, window_size) ||
                check_candidate_rev(left_window_rev_lo, left_window_rev_hi, kmer_rc, kmer_rc, kmer_rc, kmer_rc, left_kmer_pos_rc, left_kmer_pos_rc, window_size) ||
                check_candidate_fwd(right_window_lo, right_window_hi, kmer, kmer, kmer, kmer, right_kmer_pos, right_kmer_pos, window_size) ||
                check_candidate_rev(right_window_rev_lo, right_window_rev_hi, kmer_rc, kmer_rc,kmer_rc, kmer_rc, right_kmer_pos_rc, right_kmer_pos_rc, window_size);
    }

    return false;
}



template<int level, int no_shapes>
inline bool RSHash::lookup_buffer128(uint64_t* buffer, uint64_t *offsets, const size_t no_skmers,
    const uint64_t kmer, const uint64_t kmer_rc, uint64_t* shapes_fwd, uint64_t* shapes_rev,
    uint64_t &text_pos, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &start_pos, uint64_t &end_pos, uint64_t &text_kmer_lo, uint64_t &text_kmer_hi)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;

    size_t s = 0;
    if(left_minimiser_pos != span-1-right_minimiser_pos) {
        for(size_t i = 0; i < no_skmers; i++) {
            if(check_shape2<level, no_shapes>(buffer, offsets[i], kmer, kmer_rc, shapes_fwd, shapes_rev, s, left_minimiser_pos, right_minimiser_pos, forward, text_pos, start_pos, end_pos, text_kmer_lo, text_kmer_hi))
                return true;
            s += span;
        }
    }
    else {
        for(size_t i = 0; i < no_skmers; i++) {
            if(check_shape<level, no_shapes>(buffer, offsets[i], kmer, kmer_rc, shapes_fwd, shapes_rev, s, left_minimiser_pos, forward, text_pos, start_pos, end_pos, text_kmer_lo, text_kmer_hi))
                return true;
            s += span;
        }
    }
    
    return false;
}


template<int no_shapes, bool use_ht, bool locate>
uint64_t RSHash::streaming_lookup1(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF;
    uint64_t current_neg_minimiser1=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* buffer1 = new uint64_t[2*(m_thres1 * span1)];
    size_t no_skmers1;
    uint64_t sequence_begin, sequence_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    uint64_t kernel, kernel_rev;
    uint64_t shapes_fwd[2*no_shapes], shapes_rev[2*no_shapes];
    uint64_t text_kmer_lo, text_kmer_hi;

    uint64_t occurences = 0;
    for(auto && window : query | rshash::views::kmer_view({.window_size = window_size}))
    {
        if constexpr (no_shapes > 0) {
            for(int i = 0; i < no_shapes; ++i) {
                const Shape64 &shape = shapes.shapes[i];
                shapes_fwd[2*i] =     window.value & shape.w_mask.lo;
                shapes_fwd[2*i + 1] = window.value_hi & shape.w_mask.hi;
                shapes_rev[2*i] =     window.value_rev & shape.w_mask_rev.lo;
                shapes_rev[2*i + 1] = window.value_rev_hi & shape.w_mask_rev.hi;
            }
        }

        if(found && extend_in_text<no_shapes>(text_pos, sequence_begin, sequence_end, forward, window.value, window.value_rev, shapes_fwd, shapes_rev, text_kmer_lo, text_kmer_hi)) {
            occurences++;
            extensions++;
            rolling = false;
        }
        else {
            if constexpr (no_shapes > 0) {
                kernel = (window.value & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
                kernel_rev = (window.value_rev & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_rev_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
            }
            else {
                kernel = window.value;
                kernel_rev = window.value_rev;
            }

            if(rolling)
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser1_position, right_minimiser1_position);
                rolling = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer128<1, no_shapes>(buffer1, offsets1, no_skmers1, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                occurences += found;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer128<1>(offsets1, buffer1, p, no_skmers1);
                found = lookup_buffer128<1, no_shapes>(buffer1, offsets1, no_skmers1, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                occurences += found;
                current_minimiser1 = minimiser1;
            }    
            else {
                occurences += lookup_last_level<no_shapes, use_ht, locate>(window.value, window.value_hi, window.value_rev, window.value_rev_hi);
                found = false;
                current_neg_minimiser1 = minimiser1;
            }
        }
    }

    delete[] offsets1;
    delete[] buffer1;
    
    return occurences;
}


template<int no_shapes, bool use_ht, bool locate>
uint64_t RSHash::streaming_lookup2(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* buffer1 = new uint64_t[2*(m_thres1 * span1)];
    uint64_t* buffer2 = new uint64_t[2*(m_thres2 * span2)];
    size_t no_skmers1, no_skmers2, no_skmers3;
    uint64_t sequence_begin, sequence_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling1 = false;
    bool rolling2 = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    size_t left_minimiser2_position, right_minimiser2_position;
    uint64_t minimiser2, minimiser2_rank;
    uint64_t kernel, kernel_rev;
    uint64_t shapes_fwd[2*no_shapes], shapes_rev[2*no_shapes];
    uint64_t text_kmer_lo, text_kmer_hi;

    uint64_t occurences = 0;
    for(auto && window : query | rshash::views::kmer_view({.window_size = window_size}))
    {
        if constexpr (no_shapes > 0) {
            for(int i = 0; i < no_shapes; ++i) {
                const Shape64 &shape = shapes.shapes[i];
                shapes_fwd[2*i] =     window.value & shape.w_mask.lo;
                shapes_fwd[2*i + 1] = window.value_hi & shape.w_mask.hi;
                shapes_rev[2*i] =     window.value_rev & shape.w_mask_rev.lo;
                shapes_rev[2*i + 1] = window.value_rev_hi & shape.w_mask_rev.hi;
            }
        }

        if(found && extend_in_text<no_shapes>(text_pos, sequence_begin, sequence_end, forward, window.value, window.value_rev, shapes_fwd, shapes_rev, text_kmer_lo, text_kmer_hi)) {
            occurences++;
            extensions++;
            rolling1 = false;
            rolling2 = false;
        }
        else {
            if constexpr (no_shapes > 0) {
                kernel = (window.value & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
                kernel_rev = (window.value_rev & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_rev_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
            }
            else {
                kernel = window.value;
                kernel_rev = window.value_rev;
            }

            if(rolling1)
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser1_position, right_minimiser1_position);
                rolling1 = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer128<1, no_shapes>(buffer1, offsets1, no_skmers1, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                occurences += found;
                rolling2 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer128<1>(offsets1, buffer1, p, no_skmers1);
                found = lookup_buffer128<1, no_shapes>(buffer1, offsets1, no_skmers1, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                occurences += found;
                current_minimiser1 = minimiser1;
                rolling2 = false;
            }
            else {
                if(rolling2)
                    update_minimiser<2>(kernel, kernel_rev, minimiser2, left_minimiser2_position, right_minimiser2_position);
                else {
                    minimiser2 = find_minimiser<2>(kernel, kernel_rev, left_minimiser2_position, right_minimiser2_position);
                    rolling2 = true;
                }

                if(minimiser2 == current_minimiser2) {
                    found = lookup_buffer128<2, no_shapes>(buffer2, offsets2, no_skmers2, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                    occurences += found;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser2_rank)) {
                    const size_t p = s2_select.select(minimiser2_rank);
                    no_skmers2 = s2_select.select(minimiser2_rank+1) - p;

                    fill_buffer128<2>(offsets2, buffer2, p, no_skmers2);
                    found = lookup_buffer128<2, no_shapes>(buffer2, offsets2, no_skmers2, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                    occurences += found;
                    current_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                }   
                else {
                    occurences += lookup_last_level<no_shapes, use_ht, locate>(window.value, window.value_hi, window.value_rev, window.value_rev_hi);
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


template<int no_shapes, bool use_ht, bool locate>
uint64_t RSHash::streaming_lookup3(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF, current_minimiser3=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF, current_neg_minimiser3=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* offsets3 = new uint64_t[m_thres3];
    uint64_t* buffer1 = new uint64_t[2*(m_thres1 * span1)];
    uint64_t* buffer2 = new uint64_t[2*(m_thres2 * span2)];
    uint64_t* buffer3 = new uint64_t[2*(m_thres3 * span3)];
    size_t no_skmers1, no_skmers2, no_skmers3;
    uint64_t sequence_begin, sequence_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling1 = false;
    bool rolling2 = false;
    bool rolling3 = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    size_t left_minimiser2_position, right_minimiser2_position;
    uint64_t minimiser2, minimiser2_rank;
    size_t left_minimiser3_position, right_minimiser3_position;
    uint64_t minimiser3, minimiser3_rank;
    uint64_t text_kmer_lo, text_kmer_hi;
    uint64_t kernel, kernel_rev;
    uint64_t shapes_fwd[2*no_shapes], shapes_rev[2*no_shapes];

    uint64_t occurences = 0;
    for(auto && window : query | rshash::views::kmer_view({.window_size = window_size}))
    {
        if constexpr (no_shapes > 0) {
            for(int i = 0; i < no_shapes; ++i) {
                const Shape64 &shape = shapes.shapes[i];
                shapes_fwd[2*i] =     window.value & shape.w_mask.lo;
                shapes_fwd[2*i + 1] = window.value_hi & shape.w_mask.hi;
                shapes_rev[2*i] =     window.value_rev & shape.w_mask_rev.lo;
                shapes_rev[2*i + 1] = window.value_rev_hi & shape.w_mask_rev.hi;
            }
        }

        if(found && extend_in_text<no_shapes>(text_pos, sequence_begin, sequence_end, forward, window.value, window.value_rev, shapes_fwd, shapes_rev, text_kmer_lo, text_kmer_hi)) {
            occurences++;
            extensions++;
            rolling1 = false;
            rolling2 = false;
            rolling3 = false;
        }
        else {
            if constexpr (no_shapes > 0) {
                kernel = (window.value & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
                kernel_rev = (window.value_rev & shapes.kernel_mask.lo) >> 2*shapes.overlap | ((window.value_rev_hi & shapes.kernel_mask.hi) << shapes.kernel_length_lo);
            }
            else {
                kernel = window.value;
                kernel_rev = window.value_rev;
            }

            if(rolling1)
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser1_position, right_minimiser1_position);
                rolling1 = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer128<1, no_shapes>(buffer1, offsets1, no_skmers1, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                occurences += found;
                rolling2 = false;
                rolling3 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer128<1>(offsets1, buffer1, p, no_skmers1);
                found = lookup_buffer128<1, no_shapes>(buffer1, offsets1, no_skmers1, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                occurences += found;
                current_minimiser1 = minimiser1;
                rolling2 = false;
                rolling3 = false;
            }
            else {
                if(rolling2)
                    update_minimiser<2>(kernel, kernel_rev, minimiser2, left_minimiser2_position, right_minimiser2_position);
                else {
                    minimiser2 = find_minimiser<2>(kernel, kernel_rev, left_minimiser2_position, right_minimiser2_position);
                    rolling2 = true;
                }

                if(minimiser2 == current_minimiser2) {
                    found = lookup_buffer128<2, no_shapes>(buffer2, offsets2, no_skmers2, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                    occurences += found;
                    rolling3 = false;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser2_rank)) {
                    const size_t p = s2_select.select(minimiser2_rank);
                    no_skmers2 = s2_select.select(minimiser2_rank+1) - p;

                    fill_buffer128<2>(offsets2, buffer2, p, no_skmers2);
                    found = lookup_buffer128<2, no_shapes>(buffer2, offsets2, no_skmers2, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                    occurences += found;
                    current_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                    rolling3 = false;
                }
                else {
                    if(rolling3)
                        update_minimiser<3>(kernel, kernel_rev, minimiser3, left_minimiser3_position, right_minimiser3_position);
                    else {
                        minimiser3 = find_minimiser<3>(kernel, kernel_rev, left_minimiser3_position, right_minimiser3_position);
                        rolling3 = true;
                    }

                    if(minimiser3 == current_minimiser3) {
                        found = lookup_buffer128<3, no_shapes>(buffer3, offsets3, no_skmers3, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser3_position, right_minimiser3_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                        occurences += found;
                    }
                    else if(minimiser3 != current_neg_minimiser3 && r3.contains(minimiser3, minimiser3_rank)) {
                        const size_t p = s3_select.select(minimiser3_rank);
                        no_skmers3 = s3_select.select(minimiser3_rank+1) - p;

                        fill_buffer128<3>(offsets3, buffer3, p, no_skmers3);
                        found = lookup_buffer128<3, no_shapes>(buffer3, offsets3, no_skmers3, window.value, window.value_rev, shapes_fwd, shapes_rev, text_pos, left_minimiser3_position, right_minimiser3_position, forward, sequence_begin, sequence_end, text_kmer_lo, text_kmer_hi);
                        occurences += found;
                        current_minimiser3 = minimiser3;
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                    }
                    else {
                        occurences += lookup_last_level<no_shapes, use_ht, locate>(window.value, window.value_hi, window.value_rev, window.value_rev_hi);
                        found = false;
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                        current_neg_minimiser3 = minimiser3;
                    }
                }
            }
        }
    }

    delete[] offsets1;
    delete[] offsets2;
    delete[] offsets3;
    delete[] buffer1;
    delete[] buffer2;
    delete[] buffer3;
    
    return occurences;
}