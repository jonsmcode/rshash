#include <bitset>
#include <stdint.h>
#include <immintrin.h>
#include "util.hpp"



static inline uint64_t compute_shape_mask(uint32_t const shape) {
    uint64_t x = _pdep_u64(shape, 0x5555555555555555ULL);
    return x | (x << 1);
}

static inline mask128_t compute_shape_mask(uint64_t shape)
{
    uint32_t lo_shape = static_cast<uint32_t>(shape);
    uint32_t hi_shape = static_cast<uint32_t>(shape >> 32);

    uint64_t lo = _pdep_u64(lo_shape, 0x5555555555555555ULL);
    uint64_t hi = _pdep_u64(hi_shape, 0x5555555555555555ULL);

    return {
        lo | (lo << 1),
        hi | (hi << 1)
    };
}

static inline constexpr uint32_t bit_length(uint32_t x) {
    return x == 0 ? 0 : 32 - __builtin_clz(x);
}

static inline constexpr uint32_t bit_length(uint64_t x) {
    return std::bit_width(x);
}

static inline constexpr uint32_t reverse32(uint32_t x) {
    x = ((x >> 1)  & 0x55555555) | ((x & 0x55555555) << 1);
    x = ((x >> 2)  & 0x33333333) | ((x & 0x33333333) << 2);
    x = ((x >> 4)  & 0x0F0F0F0F) | ((x & 0x0F0F0F0F) << 4);
    x = ((x >> 8)  & 0x00FF00FF) | ((x & 0x00FF00FF) << 8);
    x = (x >> 16) | (x << 16);
    return x;
}

static inline constexpr uint64_t reverse64(uint64_t x)
{
    x = ((x >> 1)  & 0x5555555555555555ULL) | ((x & 0x5555555555555555ULL) << 1);
    x = ((x >> 2)  & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
    x = ((x >> 4)  & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
    x = ((x >> 8)  & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
    x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);

    x = (x >> 32) | (x << 32);

    return x;
}


static inline constexpr uint32_t reverse_shape(uint32_t x) {
    assert(x != 0);
    uint32_t len = 32 - __builtin_clz(x);
    return reverse32(x) >> (32 - len);
}

static inline constexpr uint64_t reverse_shape(uint64_t x) {
    assert(x != 0);
    uint64_t len = std::bit_width(x);
    return reverse64(x) >> (64 - len);
}


typedef struct {
    unsigned start;
    unsigned end;
    unsigned len;
} run_t;

static inline constexpr run_t find_long_run(uint32_t shape)
{
    run_t best{32, 32, 0};

    unsigned pos = 0;
    uint32_t x = shape;

    while (x) {
        unsigned zeros = std::countr_zero(x);
        pos += zeros;
        x >>= zeros;

        unsigned len = std::countr_zero(~x);

        if (len > best.len) {
            best.start = pos;
            best.len = len;
            best.end = pos + len;
        }

        pos += len;
        x >>= len;
    }

    return best;
}

static inline constexpr run_t find_long_run(uint64_t shape)
{
    run_t best{64, 64, 0};

    unsigned pos = 0;
    uint64_t x = shape;

    while (x) {
        unsigned zeros = std::countr_zero(x);
        pos += zeros;
        x >>= zeros;

        unsigned len = std::countr_zero(~x);

        if (len > best.len) {
            best.start = pos;
            best.len = len;
            best.end = pos + len;
        }

        pos += len;
        x >>= len;
    }

    return best;
}


static inline constexpr bool canonical_shape(uint32_t x)
{
    if (x == 0)
        return true;

    unsigned n = std::bit_width(x);

    for (unsigned i = 0; i < n / 2; ++i) {
        if (((x >> i) & 1) != ((x >> (n - 1 - i)) & 1))
            return false;
    }

    return true;
}

static inline constexpr bool canonical_shape(uint64_t x)
{
    if (x == 0)
        return true;

    unsigned n = std::bit_width(x);

    return x == reverse64(x) >> (64 - n);
}



typedef struct {
    uint32_t value;
    uint64_t mask;
    uint64_t w_mask;
    unsigned weight;
    unsigned length;
    unsigned kernel_length;
    unsigned overlap;
    unsigned overlap_left;
    unsigned overlap_right;
    bool is_canonical;
} Shape32;

typedef struct {
    std::vector<Shape32> shapes;
    unsigned length;
    unsigned overlap;
    unsigned kernel_length;
    uint64_t kernel_mask;
} Shapes32;


typedef struct {
    uint64_t value;
    mask128_t mask;
    mask128_t w_mask;
    mask128_t mask_rev;
    mask128_t w_mask_rev;
    unsigned weight;
    unsigned lo_weight;
    unsigned rev_lo_weight;
    unsigned w_lo_weight;
    unsigned w_rev_lo_weight;
    unsigned w_dist_right;
    unsigned w_dist_left;
    unsigned w_rev_dist_right;
    unsigned w_rev_dist_left;
    unsigned length;
    unsigned kernel_length;
    unsigned overlap;
    unsigned overlap_left;
    unsigned overlap_right;
    bool is_canonical;
} Shape64;

typedef struct {
    std::vector<Shape64> shapes;
    unsigned length;
    unsigned overlap;
    unsigned kernel_length;
    mask128_t kernel_mask;
    unsigned kernel_length_lo;
    unsigned w_dist;
} Shapes64;



static inline void print_shape(const Shape32 &shape) {
    std::cout << "Shape value: " << std::bitset<32>(shape.value) << "\n";
    std::cout << "Mask: " << std::bitset<64>(shape.mask) << "\n";
    std::cout << "aligned ask: " << std::bitset<64>(shape.w_mask) << "\n";
    std::cout << "Weight: " << shape.weight << "\n";
    std::cout << "Length: " << shape.length << "\n";
    std::cout << "Overlap: " << shape.overlap<< "\n";
    std::cout << "Overlap Left: " << shape.overlap_left << "\n";
    std::cout << "Overlap Right: " << shape.overlap_right << "\n";
}

static inline void print_shapes(const Shapes32 &shapes) {
    for (const Shape32 &shape : shapes.shapes)
        print_shape(shape);
    std::cout << "Shapes length: " << shapes.length << "\n";
    std::cout << "Shapes overlap: " << shapes.overlap << "\n";
    std::cout << "Shapes kernel length: " << shapes.kernel_length << "\n";
    std::cout << "Shapes kernel mask: " << std::bitset<64>(shapes.kernel_mask) << "\n";
}

static inline void print_mask128(const mask128_t &mask) {
    std::cout << std::bitset<64>(mask.hi) << " " << std::bitset<64>(mask.lo);
}

static inline void print_shape(const Shape64 &shape)
{
    std::cout << "Shape value: "  << std::bitset<64>(shape.value) << "\n";
    std::cout << "Mask: ";
    print_mask128(shape.mask);
    std::cout << "\n";
    std::cout << "aligned mask: ";
    print_mask128(shape.w_mask);
    std::cout << "\n";
    std::cout << "Mask rev: ";
    print_mask128(shape.mask_rev);
    std::cout << "\n";
    std::cout << "aligned rev mask: ";
    print_mask128(shape.w_mask_rev);
    std::cout << "\n";
    std::cout << "Weight: " << shape.weight << "\n";
    std::cout << "Weight mask low: " << shape.w_lo_weight << "\n";
    std::cout << "Weight rev mask low: " << shape.w_rev_lo_weight << "\n";
    std::cout << "left dist aligned mask: " << shape.w_dist_left << "\n";
    std::cout << "right dist aligned mask: " << shape.w_dist_right << "\n";
    std::cout << "left dist aligned rev mask: " << shape.w_rev_dist_left << "\n";
    std::cout << "right dist aligned rev mask: " << shape.w_rev_dist_right << "\n";
    std::cout << "Length: " << shape.length << "\n";
    std::cout << "Kernel length: " << shape.kernel_length << "\n";
    std::cout << "Overlap: " << shape.overlap << "\n";
    std::cout << "Overlap Left: " << shape.overlap_left << "\n";
    std::cout << "Overlap Right: " << shape.overlap_right << "\n";
    run_t run = find_long_run(shape.value);
    std::cout << "Shapes run start: " << run.start << "\n";
    std::cout << "Shapes run end: " << run.end << "\n";
}

static inline void print_shapes(const Shapes64 &shapes)
{
    for (const Shape64 &shape : shapes.shapes)
        print_shape(shape);

    std::cout << "Shapes length: " << shapes.length << "\n";
    std::cout << "Shapes overlap: " << shapes.overlap << "\n";
    std::cout << "Shapes kernel length: " << shapes.kernel_length << "\n";
    std::cout << "Shapes kernel length low: " << shapes.kernel_length_lo << "\n";
    std::cout << "Shapes kernel mask: ";
    print_mask128(shapes.kernel_mask);
    std::cout << "\n";
}


static inline Shape32 shape32_create(uint32_t value) {
    Shape32 shape;
    shape.value = value;

    if(value == std::numeric_limits<uint32_t>::max()) {
        shape.mask = std::numeric_limits<uint64_t>::max();
        shape.weight = 32;
        shape.length = 32;
        shape.overlap = 0;
        return shape;
    }
    
    shape.mask = compute_shape_mask(value);
    shape.weight = __builtin_popcount(value);
    shape.length = bit_length(value);
    run_t run = find_long_run(value);
    shape.kernel_length = run.len;
    shape.overlap_right = run.start;
    shape.overlap_left = shape.length - run.end;
    shape.overlap = std::max(shape.overlap_left, shape.overlap_right);
    shape.is_canonical = canonical_shape(value);

    return shape;
}

static inline Shape64 shape64_create(uint64_t value)
{
    Shape64 shape;
    shape.value = value;

    if (value == UINT64_MAX) {
        shape.mask.lo = UINT64_MAX;
        shape.mask.hi = UINT64_MAX;
        shape.weight = 64;
        shape.length = 64;
        shape.overlap = 0;
        shape.w_mask = shape.mask;
        shape.w_lo_weight = 0;
        return shape;
    }

    shape.weight = std::popcount(value);
    if(shape.weight > 32) {
        std::cerr << "shape weight > 32 not supported\n";
        exit(1);
    }

    shape.mask = compute_shape_mask(value);
    shape.mask_rev = compute_shape_mask(reverse_shape(value));
    shape.lo_weight = std::popcount(shape.mask.lo);
    shape.rev_lo_weight = std::popcount(shape.mask_rev.lo);
    shape.length = std::bit_width(value);
    run_t run = find_long_run(value);
    shape.kernel_length = run.len;
    shape.overlap_right = run.start;
    shape.overlap_left = shape.length - run.end;
    shape.overlap = std::max(shape.overlap_left, shape.overlap_right);
    shape.is_canonical = canonical_shape(value);
    
    return shape;
}

static inline constexpr mask128_t mask128_shl(mask128_t x, unsigned shift) {
    if (shift == 0)
        return x;
    if (shift < 64)
        return {x.lo << shift, (x.hi << shift) | (x.lo >> (64 - shift))};
    if (shift < 128)
        return {0, x.lo << (shift - 64)};
    return {0, 0};
}

static inline void align_shapes(Shapes32 &shapes) {
    for (Shape32 &shape : shapes.shapes) {
        unsigned shift = 2 * (shapes.overlap - shape.overlap_right);
        shape.w_mask = shape.mask << shift;
    }
}

static inline void align_shapes(Shapes64 &shapes)
{
    for (Shape64 &shape : shapes.shapes) {
        shape.w_dist_right = shapes.overlap - shape.overlap_right;
        shape.w_mask = mask128_shl(shape.mask, 2 * shape.w_dist_right);
        shape.w_dist_left = (64 -std::bit_width(shape.w_mask.hi))/2;
        shape.w_lo_weight = std::popcount(shape.w_mask.lo);

        shape.w_rev_dist_right = shapes.overlap - shape.overlap_left;
        shape.w_mask_rev = mask128_shl(shape.mask_rev, 2 * shape.w_rev_dist_right);
        shape.w_rev_lo_weight = std::popcount(shape.w_mask_rev.lo);
        shape.w_rev_dist_left = (64 -std::bit_width(shape.w_mask_rev.hi))/2;

        shapes.w_dist = std::max(shape.w_dist_left, shape.w_dist_right);
        shapes.w_dist = std::max(shapes.w_dist, shape.w_rev_dist_right);
        shapes.w_dist = std::max(shapes.w_dist, shape.w_rev_dist_left);
    }
}

static inline Shapes32 shape32_create(const std::vector<uint32_t> &values)
{
    Shapes32 window;
    window.overlap = 0;
    window.kernel_length = 32;
    for(uint32_t shape_val : values) {
        // if(!canonical_shape(shape_val)) {
        //     std::cerr << "shape " << std::bitset<32>(shape_val) << " is not canonical\n";
        //     exit(1);
        // }
        Shape32 shape = shape32_create(shape_val);
        window.shapes.emplace_back(shape);
        window.overlap = std::max(window.overlap, shape.overlap);
        window.kernel_length = std::min(window.kernel_length, shape.kernel_length);
    }
    
    window.length = window.kernel_length + 2*window.overlap;
    if(window.length > 32) {
        std::cerr << "shapes length > 32 not supported\n";
        exit(1);
    }
    window.kernel_mask = compute_mask(2u * window.kernel_length) << (2 * window.overlap);

    align_shapes(window);
    print_shapes(window);

    return window;
}

static inline Shapes64 shape64_create(const std::vector<uint64_t> &values)
{
    Shapes64 window;
    window.overlap = 0;
    window.kernel_length = 64;
    for (uint64_t shape_val : values) {
        // if (!canonical_shape(shape_val)) {
        //     std::cerr << "shape " << std::bitset<64>(shape_val) << " is not canonical\n";
        //     exit(1);
        // }

        Shape64 shape = shape64_create(shape_val);
        window.shapes.emplace_back(shape);
        window.overlap = std::max(window.overlap, shape.overlap);

        window.kernel_length = std::min(window.kernel_length, shape.kernel_length);
    }

    window.length = window.kernel_length + 2 * window.overlap;

    if (window.length > 64) {
        std::cerr << "aligned shapes length > 64 not supported\n";
        exit(1);
    }

    window.kernel_mask = mask128_shl(compute_mask128(2 * window.kernel_length), 2 * window.overlap);
    window.kernel_length_lo = std::popcount(window.kernel_mask.lo);
    align_shapes(window);
    print_shapes(window);

    return window;
}


namespace cereal {

template <class Archive>
void serialize(Archive& ar, Shape32& shape) {
    ar(shape.value,
       shape.mask,
       shape.w_mask,
       shape.weight,
       shape.length,
       shape.overlap,
       shape.is_canonical);
}

template <class Archive>
void serialize(Archive& ar, Shapes32& window) {
    ar(window.shapes,
       window.overlap,
       window.length,
       window.kernel_length,
       window.kernel_mask);
}

template <class Archive>
void serialize(Archive& ar, Shape64& shape) {
    ar(shape.value,
       shape.mask.lo,
       shape.mask.hi,
       shape.w_mask.lo,
       shape.w_mask.hi,
       shape.mask_rev.lo,
       shape.mask_rev.hi,
       shape.w_mask_rev.lo,
       shape.w_mask_rev.hi,
       shape.w_lo_weight,
       shape.w_rev_lo_weight,
       shape.lo_weight,
       shape.rev_lo_weight,
       shape.w_dist_right,
       shape.w_dist_left,
       shape.w_rev_dist_right,
       shape.w_rev_dist_left,
       shape.weight,
       shape.length,
       shape.overlap,
       shape.overlap_left,
       shape.overlap_right,
       shape.is_canonical);
}

template <class Archive>
void serialize(Archive& ar, Shapes64& window) {
    ar(window.shapes,
       window.overlap,
       window.length,
       window.kernel_length,
       window.kernel_mask.lo,
       window.kernel_mask.hi,
       window.kernel_length_lo,
       window.w_dist);
}

}