#include <cereal/cereal.hpp>
// #include <gtl/phmap.h>
// #include <sux/bits/EliasFano.hpp>
// #include <cereal/types/vector.hpp>

namespace cereal {

template <class Archive>
void save(Archive& ar, const gtl::flat_hash_set<uint64_t>& set) {
    ar(static_cast<std::size_t>(set.size()));
    for (const auto& elem : set) {
        ar(elem);
    }
}

template <class Archive>
void load(Archive& ar, gtl::flat_hash_set<uint64_t>& set) {
    std::size_t size;
    ar(size);
    set.clear();
    set.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        uint64_t elem;
        ar(elem);
        set.insert(std::move(elem));
    }
}

} // namespace cereal


// template <sux::util::AllocType AT = sux::util::AllocType::MALLOC>
// struct EliasFanoCereal {
//     std::vector<uint64_t> lower_bits;
//     std::vector<uint64_t> upper_bits;
//     uint64_t num_bits;
//     uint64_t num_ones;
//     int l;
//     int block_size;
//     int block_length;
//     uint64_t block_size_mask;
//     uint64_t lower_l_bits_mask;

//     EliasFanoCereal() = default;

//     EliasFanoCereal(const sux::bits::EliasFano<AT>& ef)
//         : lower_bits(ef.lower_bits.begin(), ef.lower_bits.end()),
//           upper_bits(ef.upper_bits.begin(), ef.upper_bits.end()),
//           num_bits(ef.num_bits),
//           num_ones(ef.num_ones),
//           l(ef.l),
//           block_size(ef.block_size),
//           block_length(ef.block_length),
//           block_size_mask(ef.block_size_mask),
//           lower_l_bits_mask(ef.lower_l_bits_mask)
//     {}

//     template<class Archive>
//     void serialize(Archive& ar) {
//         ar(lower_bits, upper_bits, num_bits, num_ones, l, block_size, block_length, block_size_mask, lower_l_bits_mask);
//     }
// };


// template <sux::util::AllocType AT = sux::util::AllocType::MALLOC>
// class EliasFanoRestorer : public sux::bits::EliasFano<AT> {
// public:
//     using Base = sux::bits::EliasFano<AT>;
//     using Base::lower_bits;
//     using Base::upper_bits;
//     using Base::num_bits;
//     using Base::num_ones;
//     using Base::l;
//     using Base::block_size;
//     using Base::block_length;
//     using Base::block_size_mask;
//     using Base::lower_l_bits_mask;
//     using Base::select_upper;
//     using Base::selectz_upper;

//     EliasFanoRestorer(const EliasFanoCereal<AT>& efc)
//         : Base() // default construct
//     {
//         lower_bits.resize(efc.lower_bits.size());
//         for (size_t i = 0; i < efc.lower_bits.size(); ++i)
//             lower_bits[i] = efc.lower_bits[i];

//         upper_bits.resize(efc.upper_bits.size());
//         for (size_t i = 0; i < efc.upper_bits.size(); ++i)
//             upper_bits[i] = efc.upper_bits[i];

//         num_bits = efc.num_bits;
//         num_ones = efc.num_ones;
//         l = efc.l;
//         block_size = efc.block_size;
//         block_length = efc.block_length;
//         block_size_mask = efc.block_size_mask;
//         lower_l_bits_mask = efc.lower_l_bits_mask;

//         // Rebuild select structures
//         select_upper = sux::bits::SimpleSelectHalf<AT>(&upper_bits, num_ones + (num_bits >> l) + 1);
//         selectz_upper = sux::bits::SimpleSelectZeroHalf<AT>(&upper_bits, num_ones + (num_bits >> l) + 1);
//     }
// };


