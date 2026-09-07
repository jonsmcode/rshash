#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>
#include "rshash.hpp"



template <bool use_shape, bool use_ht>
void RSHash::initialise_locatefn_impl()
{
    switch (level) {
    case 1:
        streaming_locate_fn = &RSHash::streaming_locate1<use_shape, use_ht>;
        locate_fn           = &RSHash::locate1<use_shape, use_ht>;
        break;

    case 2:
        streaming_locate_fn = &RSHash::streaming_locate2<use_shape, use_ht>;
        locate_fn           = &RSHash::locate2<use_shape, use_ht>;
        break;

    case 3:
        streaming_locate_fn = &RSHash::streaming_locate3<use_shape, use_ht>;
        locate_fn           = &RSHash::locate3<use_shape, use_ht>;
        break;

    default:
        streaming_locate_fn = nullptr;
        locate_fn = nullptr;
    }
}

void RSHash::initialise_locatefn()
{
    const bool use_shape = shapes.shapes[0].value != std::numeric_limits<uint32_t>::max();
    if (use_shape) {
        if (use_ht)
            initialise_locatefn_impl<true, true>();
        else
            initialise_locatefn_impl<true, false>();
    }
    else {
        if (use_ht)
            initialise_locatefn_impl<false, true>();
        else
            initialise_locatefn_impl<false, false>();
    }
}

uint64_t RSHash::streaming_locate(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions) {
    return (this->*streaming_locate_fn)(query, positions, no_positions);
}

uint64_t RSHash::locate(const std::vector<uint64_t> &kmers, std::vector<uint64_t> &positions) {
    return (this->*locate_fn)(kmers, positions);
}


template<bool use_shape, bool use_ht>
inline bool RSHash::locate_last_level(const uint64_t kmer, const uint64_t kmer_rc,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    bool found = false;

    if constexpr (use_ht) {
        auto locate = [&](auto &map, uint64_t key) {
            if(auto it = map.find(key); !it.empty()) {
                for(uint32_t v : it)
                    no_positions++;
                    // positions.push_back(v);
                // std::memcpy(positions + pos, span.data(), span.size_bytes());
                found = true;
            }
        };
        // if constexpr (use_shape) {
        //     locate(hashmap, kmer);
        //     if(kmer != kmer_rc)
        //         locate(hashmap_rc, kmer_rc);
        // }
        // else
        //     locate(hashmap, std::min<uint64_t>(kmer, kmer_rc));
    }
    else {
        auto locate = [&](auto &ef, auto &sel, uint64_t key) {
            uint64_t rank;
            if(ef.contains(key, rank)) {
                const size_t pos = sel.select(rank);
                const size_t occs = sel.select(rank + 1) - pos;
                no_positions += occs;
                found = true;
            }
        };
        if constexpr (use_shape) {
            locate(r4, s4_select, kmer);
            if(kmer != kmer_rc)
                locate(r5, s5_select, kmer_rc);
        }
        else
            locate(r4, s4_select, std::min<uint64_t>(kmer, kmer_rc));
    }

    return found;
}


template<bool use_shape, bool use_ht>
uint64_t RSHash::locate1(const std::vector<uint64_t> &kmers, std::vector<uint64_t> &positions)
{
    uint64_t no_positions;
    uint64_t* offsets = new uint64_t[m_thres1];
    uint64_t minimiser, minimiser_rank, kmer_rank;;
    uint64_t kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes

    for(uint64_t kmer : kmers) {
        uint64_t kmer_rc = crc(kmer, window_size);

        if constexpr (use_shape) {
            kernel = (kmer & shapes.kernel_mask) >> 2*shapes.overlap;
            kernel_rev = (kmer_rc & shapes.kernel_mask) >> 2*shapes.overlap;
            kmer = _pext_u64(kmer, shape.mask);
            kmer_rc = _pext_u64(kmer_rc, shape.mask);
        }
        else {
            kernel = kmer;
            kernel_rev = kmer_rc;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            no_positions += check_pos<1, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            locate_last_level<use_shape, use_ht>(kmer, kmer_rc, positions, no_positions);
        }
    }

    delete[] offsets;

    return no_positions;
}

template<bool use_shape, bool use_ht>
uint64_t RSHash::locate2(const std::vector<uint64_t> &kmers, std::vector<uint64_t> &positions)
{
    uint64_t no_positions;
    uint64_t* offsets = new uint64_t[m_thres2];
    uint64_t minimiser, minimiser_rank, kmer_rank;
    uint64_t kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes

    for(uint64_t kmer : kmers) {
        uint64_t kmer_rc = crc(kmer, window_size);

        if constexpr (use_shape) {
            kernel = (kmer & shapes.kernel_mask) >> 2*shapes.overlap;
            kernel_rev = (kmer_rc & shapes.kernel_mask) >> 2*shapes.overlap;
            kmer = _pext_u64(kmer, shape.mask);
            kmer_rc = _pext_u64(kmer_rc, shape.mask);
        }
        else {
            kernel = kmer;
            kernel_rev = kmer_rc;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            no_positions += check_pos<1, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            minimiser = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
            if(r2.contains(minimiser, minimiser_rank)) {
                size_t p = s2_select.select(minimiser_rank);
                size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
                no_positions += check_pos<2, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
            }
            else {
                locate_last_level<use_shape, use_ht>(kmer, kmer_rc, positions, no_positions);
            }
        }
    }

    delete[] offsets;

    return no_positions;
}

template<bool use_shape, bool use_ht>
uint64_t RSHash::locate3(const std::vector<uint64_t> &kmers, std::vector<uint64_t> &positions)
{
    uint64_t no_positions;
    uint64_t* offsets = new uint64_t[m_thres3];
    uint64_t minimiser, minimiser_rank, kmer_rank;
    uint64_t kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes

    for(uint64_t kmer : kmers) {
        uint64_t kmer_rc = crc(kmer, window_size);

        if constexpr (use_shape) {
            kernel = (kmer & shapes.kernel_mask) >> 2*shapes.overlap;
            kernel_rev = (kmer_rc & shapes.kernel_mask) >> 2*shapes.overlap;
            kmer = _pext_u64(kmer, shape.mask);
            kmer_rc = _pext_u64(kmer_rc, shape.mask);
        }
        else {
            kernel = kmer;
            kernel_rev = kmer_rc;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);

        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            no_positions += check_pos<1, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            minimiser = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
            if(r2.contains(minimiser, minimiser_rank)) {
                size_t p = s2_select.select(minimiser_rank);
                size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
                no_positions += check_pos<2, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
            }
            else {
                minimiser = find_minimiser<3>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
                if(r3.contains(minimiser, minimiser_rank)) {
                    size_t p = s3_select.select(minimiser_rank);
                    size_t no_minimiser = s3_select.select(minimiser_rank+1) - p;
                    no_positions += check_pos<3, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
                }
                else {
                    locate_last_level<use_shape, use_ht>(kmer, kmer_rc, positions, no_positions);
                }
            }
        }
    }

    delete[] offsets;

    return no_positions;
}


template<int level, bool use_shape>
inline uint64_t RSHash::check_pos(const uint64_t kmer, const uint64_t kmer_rc,
    uint64_t* offsets, const size_t p, const size_t no_skmers,
    const size_t left_minimiser_position, const size_t right_minimiser_position)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;
    
    for(size_t i = 0; i < no_skmers; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i)-span+1;
        if constexpr (level == 2)
            offsets[i] = offsets2.access(p+i)-span+1;
        if constexpr (level == 3)
            offsets[i] = offsets3.access(p+i)-span+1;
    }

    uint64_t positions = 0;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes
    for(size_t i = 0; i < no_skmers; i++) {
        uint64_t offset = offsets[i];
        uint64_t pos = span-1-left_minimiser_position;
        uint64_t pos_rc = left_minimiser_position;

        uint64_t hash_fwd = get_word64(offset + pos - shape.overlap) & windowmask;
        uint64_t hash_rc = get_word64(offset + pos_rc - shape.overlap) & windowmask;

        if constexpr (use_shape) {
            hash_fwd = _pext_u64(hash_fwd, shape.mask);
            hash_rc = _pext_u64(hash_rc, shape.mask);
        }

        positions += kmer == hash_fwd;
        positions += kmer_rc == hash_rc;

        if(left_minimiser_position != k-m1-right_minimiser_position) {
            pos = right_minimiser_position;
            pos_rc = span-1-right_minimiser_position;

            hash_fwd = get_word64(offset + pos - shape.overlap) & windowmask;
            hash_rc = get_word64(offset + pos_rc - shape.overlap) & windowmask;

            if constexpr (use_shape) {
                hash_fwd = _pext_u64(hash_fwd, shape.mask);
                hash_rc = _pext_u64(hash_rc, shape.mask);
            }

            positions += kmer == hash_fwd;
            positions += kmer_rc == hash_rc;
        }

    }

    return positions;
}


template<int level, bool use_shape>
inline bool RSHash::report_minimiser_pos(uint64_t *buffer, uint64_t offset,
    const uint64_t query, const uint64_t queryrc, uint64_t &no_positions, const size_t s,
    const size_t minimiser_pos, const uint64_t start_pos, const uint64_t end_pos,
    std::vector<uint64_t> &positions)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;

    uint64_t candidate, candidate_rc, pos, pos_rc;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes
    if constexpr (use_shape) {
        pos = overlap - shape.overlap + span-1-minimiser_pos;
        pos_rc = overlap - shape.overlap + minimiser_pos;
        offset -= overlap;
        candidate = _pext_u64(buffer[s + pos], shape.mask);
        candidate_rc = _pext_u64(buffer[s + pos_rc], shape.mask);
    }
    else {
        pos = span-1-minimiser_pos;
        pos_rc = minimiser_pos;
        candidate = buffer[s + pos];
        candidate_rc = buffer[s + pos_rc];
    }

    if(candidate == query && offset + pos >= start_pos && offset + pos + window_size - 1 < end_pos) {
        // forward = true;
        // positions.push_back(offset + pos + window_size - 1);
        no_positions++;
        return true;
    }
    if(candidate_rc == queryrc && offset + pos_rc >= start_pos && offset + pos_rc + window_size - 1 < end_pos) {
        // forward = false;
        // positions.push_back(offset + pos_rc);
        no_positions++;
        return true;
    }

    return false;
}


template<int level, bool use_shape>
inline bool RSHash::report_minimiser_pos2(uint64_t *buffer, uint64_t offset,
    const uint64_t query, const uint64_t queryrc, uint64_t &no_positions, const size_t s,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    const uint64_t start_pos, const uint64_t end_pos, std::vector<uint64_t> &positions)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;

    uint64_t left_candidate, left_candidate_rc, right_candidate, right_candidate_rc;
    uint64_t left_pos, left_pos_rc, right_pos, right_pos_rc;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes
    if constexpr (use_shape) {
        left_pos = overlap - shape.overlap + span-1-left_minimiser_pos;
        left_pos_rc = overlap - shape.overlap + left_minimiser_pos;
        right_pos = overlap - shape.overlap + right_minimiser_pos;
        right_pos_rc = overlap - shape.overlap + span-1-right_minimiser_pos;
        left_candidate = _pext_u64(buffer[s + left_pos], shape.mask);
        left_candidate_rc = _pext_u64(buffer[s + left_pos_rc], shape.mask);
        right_candidate = _pext_u64(buffer[s + right_pos], shape.mask);
        right_candidate_rc = _pext_u64(buffer[s + right_pos_rc], shape.mask);
        offset -= overlap;
    }
    else {
        left_pos = span-1-left_minimiser_pos;
        left_pos_rc = left_minimiser_pos;
        right_pos = right_minimiser_pos;
        right_pos_rc = span-1-right_minimiser_pos;
        left_candidate = buffer[s + left_pos];
        left_candidate_rc = buffer[s + left_pos_rc];
        right_candidate = buffer[s + right_pos];
        right_candidate_rc = buffer[s + right_pos_rc];
    }

    if(left_candidate == query && offset + left_pos >= start_pos && offset + left_pos + window_size - 1 < end_pos) {
        // forward = true;
        // positions.push_back(offset + left_pos + window_size - 1);
        no_positions++;
        return true;
    }
    if(left_candidate_rc == queryrc && offset + left_pos_rc >= start_pos && offset + left_pos_rc + window_size - 1 < end_pos) {
        // forward = false;
        // positions.push_back(offset + left_pos_rc);
        no_positions++;
        return true;
    }
    if(right_candidate == query && offset + right_pos >= start_pos && offset + right_pos + window_size - 1 < end_pos) {
        // forward = true;
        // positions.push_back(offset + right_pos + window_size - 1);
        no_positions++;
        return true;
    }
    if(right_candidate_rc == queryrc && offset + right_pos_rc >= start_pos && offset + right_pos_rc + window_size - 1 < end_pos) {
        // forward = false;
        // positions.push_back(offset + right_pos_rc);
        no_positions++;
        return true;
    }

    return false;
}


template<int level, bool use_shape>
inline void RSHash::locate_buffer(uint64_t *buffer, uint64_t *offsets, uint64_t *sequence_ends, const size_t no_minimiser,
    const uint64_t query, const uint64_t queryrc,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    std::vector<uint64_t> &positions, uint64_t &no_positions, uint64_t &found_kmers)
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
    if constexpr (level == 3) {
        span = span3;
        m = m3;
    }

    // todo: SIMD buffer check
    bool found = false;
    size_t s = 0;
    if(left_minimiser_pos != k-m-right_minimiser_pos) {
        for(size_t i = 0; i < no_minimiser; i++) {
            found |= report_minimiser_pos2<level, use_shape>(buffer, offsets[i], query, queryrc, no_positions, s, left_minimiser_pos, right_minimiser_pos, sequence_ends[2*i], sequence_ends[2*i+1], positions);
            s += span;
        }
    }
    else {
        for(size_t i = 0; i < no_minimiser; i++) {
            found |= report_minimiser_pos<level, use_shape>(buffer, offsets[i], query, queryrc, no_positions, s, left_minimiser_pos, sequence_ends[2*i], sequence_ends[2*i+1], positions);
            s += span;
        }
    }
    
    found_kmers += found;
}



template<bool use_shape, bool use_ht>
uint64_t RSHash::streaming_locate1(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    uint64_t found_kmers = 0;
    uint64_t current_pos_minimiser=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser=std::numeric_limits<uint64_t>::max();
    uint64_t* offsets = new uint64_t[m_thres1];
    uint64_t* sequence_ends = new uint64_t[2*m_thres1];
    uint64_t* kmer_buffer = new uint64_t[m_thres1 * (span1 + overlap)];
    size_t no_minimiser;
    uint64_t minimiser, minimiser_rank;
    uint64_t kmer_rank;
    size_t left_minimiser_position, right_minimiser_position;
    bool begin = true;
    uint64_t kernel, kernel_rev, kmer, kmer_rc;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes

    for (auto&& window : query | rshash::views::kmerview({.window_size = window_size}))
    {
        if constexpr (use_shape) {
            kernel = (kmer & shapes.kernel_mask) >> 2*shapes.overlap;
            kernel_rev = (kmer_rc & shapes.kernel_mask) >> 2*shapes.overlap;
            kmer = _pext_u64(window.value, shape.mask);
            kmer_rc = _pext_u64(window.value_rev, shape.mask);
        }
        else {
            kernel = window.value;
            kmer = window.value;
            kernel_rev = window.value_rev;
            kmer_rc = window.value_rev;
        }

        if(begin) {
            minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
            begin = false;
        }
        else
            update_minimiser<1>(kernel, kernel_rev, minimiser, left_minimiser_position, right_minimiser_position);

        if(minimiser == current_pos_minimiser)
            locate_buffer<1, use_shape>(kmer_buffer, offsets, sequence_ends, no_minimiser, kmer, kmer_rc, left_minimiser_position, right_minimiser_position, positions, no_positions, found_kmers);
        else if(minimiser != current_neg_minimiser && r1.contains(minimiser, minimiser_rank)) {
            const size_t minimiser_position = s1_select.select(minimiser_rank);
            no_minimiser = s1_select.select(minimiser_rank+1) - minimiser_position;

            fill_buffer2<1>(offsets, kmer_buffer, sequence_ends, minimiser_position, no_minimiser);
            locate_buffer<1, use_shape>(kmer_buffer, offsets, sequence_ends, no_minimiser, kmer, kmer_rc, left_minimiser_position, right_minimiser_position, positions, no_positions, found_kmers);
            current_pos_minimiser = minimiser;
        }
        else {
            found_kmers += locate_last_level<use_shape, use_ht>(kmer, kmer_rc, positions, no_positions);
            current_neg_minimiser = minimiser;
        }

    }

    delete[] kmer_buffer;
    delete[] offsets;
    delete[] sequence_ends;

    return found_kmers;
}

template<bool use_shape, bool use_ht>
uint64_t RSHash::streaming_locate2(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    uint64_t found_kmers = 0;
    uint64_t current_pos_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_pos_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* kmer_buffer1 = new uint64_t[(m_thres1) * (span1 + overlap)];
    uint64_t* kmer_buffer2 = new uint64_t[(m_thres2) * (span2 + overlap)];
    uint64_t* sequence_ends1 = new uint64_t[2*m_thres1];
    uint64_t* sequence_ends2 = new uint64_t[2*m_thres2];
    size_t no_minimiser1, no_minimiser2;
    uint64_t minimiser1, minimiser_rank1;
    uint64_t minimiser2, minimiser_rank2;
    uint64_t kmer_rank;
    size_t left_minimiser_position1, right_minimiser_position1;
    size_t left_minimiser_position2, right_minimiser_position2;
    bool rolling1 = false;
    bool rolling2 = false;
    uint64_t kernel, kernel_rev, kmer, kmer_rc;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes

    for (auto&& window : query | rshash::views::kmerview({.window_size = window_size}))
    {
        if constexpr (use_shape) {
            kernel = (kmer & shapes.kernel_mask) >> 2*shapes.overlap;
            kernel_rev = (kmer_rc & shapes.kernel_mask) >> 2*shapes.overlap;
            kmer = _pext_u64(window.value, shape.mask);
            kmer_rc = _pext_u64(window.value_rev, shape.mask);
        }
        else {
            kernel = window.value;
            kmer = window.value;
            kernel_rev = window.value_rev;
            kmer_rc = window.value_rev;
        }

            if(!rolling1) {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position1, right_minimiser_position1);
                rolling1 = true;
            }
            else
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser_position1, right_minimiser_position1);

            if(minimiser1 == current_pos_minimiser1) {
                locate_buffer<1, use_shape>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, kmer, kmer_rc, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
                rolling2 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser_rank1)) {
                const size_t minimiser_position1 = s1_select.select(minimiser_rank1);
                no_minimiser1 = s1_select.select(minimiser_rank1+1) - minimiser_position1;

                fill_buffer2<1>(offsets1, kmer_buffer1, sequence_ends1, minimiser_position1, no_minimiser1);
                locate_buffer<1, use_shape>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, kmer, kmer_rc, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
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

                if (minimiser2 == current_pos_minimiser2) {
                    locate_buffer<2, use_shape>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, kmer, kmer_rc, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser_rank2)) {
                    const size_t minimiser_position2 = s2_select.select(minimiser_rank2);
                    no_minimiser2 = s2_select.select(minimiser_rank2+1) - minimiser_position2;

                    fill_buffer2<2>(offsets2, kmer_buffer2, sequence_ends2, minimiser_position2, no_minimiser2);
                    locate_buffer<2, use_shape>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, kmer, kmer_rc, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
                    current_pos_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                }
                else {
                    found_kmers += locate_last_level<use_shape, use_ht>(kmer, kmer_rc, positions, no_positions);
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

template<bool use_shape, bool use_ht>
uint64_t RSHash::streaming_locate3(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<uint64_t> &positions, uint64_t &no_positions)
{
    uint64_t found_kmers = 0;
    uint64_t current_pos_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_pos_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_pos_minimiser3=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser3=std::numeric_limits<uint64_t>::max();
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* offsets3 = new uint64_t[m_thres3];
    uint64_t* kmer_buffer1 = new uint64_t[(m_thres1) * (span1 + overlap)];
    uint64_t* kmer_buffer2 = new uint64_t[(m_thres2) * (span2 + overlap)];
    uint64_t* kmer_buffer3 = new uint64_t[(m_thres3) * (span3 + overlap)];
    uint64_t* sequence_ends1 = new uint64_t[2*m_thres1];
    uint64_t* sequence_ends2 = new uint64_t[2*m_thres2];
    uint64_t* sequence_ends3 = new uint64_t[2*m_thres3];
    size_t no_minimiser1, no_minimiser2, no_minimiser3;
    uint64_t minimiser1, minimiser_rank1;
    uint64_t minimiser2, minimiser_rank2;
    uint64_t minimiser3, minimiser_rank3;
    uint64_t kmer_rank;
    size_t left_minimiser_position1, right_minimiser_position1;
    size_t left_minimiser_position2, right_minimiser_position2;
    size_t left_minimiser_position3, right_minimiser_position3;
    bool rolling1 = false;
    bool rolling2 = false;
    bool rolling3 = false;
    uint64_t kernel, kernel_rev, kmer, kmer_rc;
    Shape32 shape = shapes.shapes[0]; // todo: multiple shapes

    for (auto&& window : query | rshash::views::kmerview({.window_size = window_size}))
    {
        if constexpr (use_shape) {
            kernel = (window.value & shapes.kernel_mask) >> 2*shapes.overlap;
            kernel_rev = (window.value_rev & shapes.kernel_mask) >> 2*shapes.overlap;
            kmer = _pext_u64(window.value, shape.mask);
            kmer_rc = _pext_u64(window.value_rev, shape.mask);
        }
        else {
            kernel = window.value;
            kmer = window.value;
            kernel_rev = window.value_rev;
            kmer_rc = window.value_rev;
        }
        
            if(!rolling1) {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position1, right_minimiser_position1);
                rolling1 = true;
            }
            else
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser_position1, right_minimiser_position1);

            if(minimiser1 == current_pos_minimiser1) {
                locate_buffer<1, use_shape>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, kmer, kmer_rc, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
                rolling2 = false;
                rolling3 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser_rank1)) {
                const size_t minimiser_position1 = s1_select.select(minimiser_rank1);
                no_minimiser1 = s1_select.select(minimiser_rank1+1) - minimiser_position1;

                fill_buffer2<1>(offsets1, kmer_buffer1, sequence_ends1, minimiser_position1, no_minimiser1);
                locate_buffer<1, use_shape>(kmer_buffer1, offsets1, sequence_ends1, no_minimiser1, kmer, kmer_rc, left_minimiser_position1, right_minimiser_position1, positions, no_positions, found_kmers);
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

                if (minimiser2 == current_pos_minimiser2) {
                    locate_buffer<2, use_shape>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, kmer, kmer_rc, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
                    rolling3 = false;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser_rank2)) {
                    const size_t minimiser_position2 = s2_select.select(minimiser_rank2);
                    no_minimiser2 = s2_select.select(minimiser_rank2+1) - minimiser_position2;

                    fill_buffer2<2>(offsets2, kmer_buffer2, sequence_ends2, minimiser_position2, no_minimiser2);
                    locate_buffer<2, use_shape>(kmer_buffer2, offsets2, sequence_ends2, no_minimiser2, kmer, kmer_rc, left_minimiser_position2, right_minimiser_position2, positions, no_positions, found_kmers);
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
                        locate_buffer<3, use_shape>(kmer_buffer3, offsets3, sequence_ends3, no_minimiser3, kmer, kmer_rc, left_minimiser_position3, right_minimiser_position3, positions, no_positions, found_kmers);
                    }
                    else if(minimiser3 != current_neg_minimiser3 && r3.contains(minimiser3, minimiser_rank3)) {
                        const size_t minimiser_position3 = s3_select.select(minimiser_rank3);
                        no_minimiser3 = s3_select.select(minimiser_rank3+1) - minimiser_position3;

                        fill_buffer2<3>(offsets3, kmer_buffer3, sequence_ends3, minimiser_position3, no_minimiser3);
                        locate_buffer<3, use_shape>(kmer_buffer3, offsets3, sequence_ends3, no_minimiser3, kmer, kmer_rc, left_minimiser_position3, right_minimiser_position3, positions, no_positions, found_kmers);
                        current_pos_minimiser3 = minimiser3;
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                    }
                    else {
                        found_kmers += locate_last_level<use_shape, use_ht>(kmer, kmer_rc, positions, no_positions);
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
