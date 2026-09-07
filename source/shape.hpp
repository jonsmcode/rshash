#include <bitset>
#include <stdint.h>
#include <immintrin.h>
#include "util.hpp"


static inline uint64_t compute_shape_mask(uint32_t const shape) {
    uint64_t x = _pdep_u64(shape, 0x5555555555555555ULL);
    return x | (x << 1);
}

static inline constexpr uint32_t bit_length(uint32_t x) {
    return x == 0 ? 0 : 32 - __builtin_clz(x);
}

static inline constexpr uint32_t reverse32(uint32_t x) {
    x = ((x >> 1)  & 0x55555555) | ((x & 0x55555555) << 1);
    x = ((x >> 2)  & 0x33333333) | ((x & 0x33333333) << 2);
    x = ((x >> 4)  & 0x0F0F0F0F) | ((x & 0x0F0F0F0F) << 4);
    x = ((x >> 8)  & 0x00FF00FF) | ((x & 0x00FF00FF) << 8);
    x = (x >> 16) | (x << 16);
    return x;
}

static inline constexpr uint32_t reverse_shape(uint32_t x) {
    assert(x != 0);

    uint32_t len = 32 - __builtin_clz(x);
    return reverse32(x) >> (32 - len);
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


static inline constexpr unsigned sum_runs(uint32_t shape, unsigned m) {
    unsigned sum = 0;
    while(shape) {
        unsigned start = __builtin_ctz(shape);
        shape >>= start;
        unsigned len = __builtin_ctz(~shape);
        if (len > m)
            sum += len;

        shape >>= len;
    }
    return sum;
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



typedef struct {
    uint32_t value;
    uint64_t mask;
    uint64_t w_mask;
    size_t weight;
    size_t length;
    size_t overlap;
} Shape32;


typedef struct {
    std::vector<Shape32> shapes;
    size_t length;
    size_t overlap;
    size_t kernel_length;
    uint64_t kernel_mask;
} Shapes32;



static inline void print_shape(const Shape32 &shape) {
    std::cout << "Shape value: " << std::bitset<32>(shape.value) << "\n";
    std::cout << "Mask: " << std::bitset<64>(shape.mask) << "\n";
    std::cout << "aligned ask: " << std::bitset<64>(shape.w_mask) << "\n";
    std::cout << "Weight: " << shape.weight << "\n";
    std::cout << "Length: " << shape.length << "\n";
    std::cout << "Overlap: " << shape.overlap<< "\n";
}

static inline void print_shapes(const Shapes32 &shapes) {
    for (const Shape32 &shape : shapes.shapes)
        print_shape(shape);
    std::cout << "Shapes length: " << shapes.length << "\n";
    std::cout << "Shapes overlap: " << shapes.overlap << "\n";
    std::cout << "Shapes kernel length: " << shapes.kernel_length << "\n";
    std::cout << "Shapes kernel mask: " << std::bitset<64>(shapes.kernel_mask) << "\n";
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
    shape.overlap = shape.length - run.end;

    return shape;
}

static inline void align_shapes(Shapes32 &shapes) {
    for (Shape32 &shape : shapes.shapes) {
        shape.w_mask = shape.mask << (2 * (shapes.overlap - shape.overlap));
    }
}

static inline Shapes32 shape32_create(const std::vector<uint32_t> &values)
{
    Shapes32 window;
    window.overlap = 0;
    window.kernel_length = 32;
    for(uint32_t shape_val : values) {
        if(!canonical_shape(shape_val)) {
            std::cerr << "shape " << std::bitset<32>(shape_val) << " is not canonical\n";
            exit(1);
        }
        Shape32 shape = shape32_create(shape_val);
        window.shapes.emplace_back(shape);
        window.overlap = std::max(window.overlap, shape.overlap);
        window.kernel_length = std::min(window.kernel_length, shape.length - 2*shape.overlap);
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


namespace cereal {

template <class Archive>
void serialize(Archive& ar, Shape32& shape) {
    ar(shape.value,
       shape.mask,
       shape.w_mask,
       shape.weight,
       shape.length,
       shape.overlap);
}

template <class Archive>
void serialize(Archive& ar, Shapes32& window) {
    ar(window.shapes,
       window.overlap,
       window.length,
       window.kernel_length,
       window.kernel_mask);
}

}