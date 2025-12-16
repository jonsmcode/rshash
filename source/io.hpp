#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>

// #include <gtl/phmap.h>
// #include <sux/bits/EliasFano.hpp>


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



