# RSHash

**RSHash** is a compressed, exact data structure for k-mers
(strings of length k over the DNA alphabet {A,C,G,T}) based on bitvectors with **R**ank and **S**elect support to **Hash** k-mers.

Given a set of DNA-strings it efficiently allows to:
- **Lookup** a k-mer x, i.e., answer if x is present in the input strings
- **Streaming Lookup** all k-mers in a query DNA-string
- **Access** a k-mer given a position in the input strings


### Setup
Checkout
```
git clone --recursive https://github.com/jonsmcode/rshash
```
Compile
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -D CMAKE_CXX_COMPILER=g++-14
make
```


### Run

See usage, command line options, and parameters with
```
.build/source/rshash --help
```
e.g.
```
.build/source/rshash build -i ../datasets/genome.fa.ust.fa.gz -d ../datasets/genome.rshash -k 31 --m1 18 --t1 64 -l 1
.build/source/rshash query -d ../datasets/genome.rshash -l 1 -q ../datasets/query.fastq.gz
```

### Benchmarks

Get datasets used by Pibiri (2022) with
```
./benchmarks/download_and_preprocess_datasets.sh
```

Run our benchmarks with script
```
./benchmarks/benchmark_rshash.sh
```

