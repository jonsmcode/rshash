# RSHash

**RSHash** is a compressed, exact data structure for k-mers
(strings of length k over the DNA alphabet {A,C,G,T}) based on bitvectors with **R**ank and **S**elect support to **Hash** k-mers.
Given a set of DNA-strings it efficiently allows to:
- **Lookup**(x) a k-mer x, i.e., answer if x is present in a string
- **Streaming Lookup**(Q) all k-mers in a query DNA-string Q longer than k
- **Access**(p) a k-mer at position p in the input strings


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
```


### Run
```
.build/source/rshash build -i ../datasets/cod.k31.unitigs.fa.ust.fa.gz -d ../datasets/cod.k31.unitigs.fa.ust.dict -k 31 --m1 18 --t1 64 -l 1
.build/source/rshash query -d ../datasets/cod.k31.unitigs.fa.ust.dict -l 1 -q ../datasets/SRR16288382_1.fastq.gz
```

### Benchmarks

Datasets

```
mkdir datasets
cd datasets
wget https://zenodo.org/records/17582116/data.fa.gz
```

