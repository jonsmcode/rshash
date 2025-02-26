# kmerdict


### Setup
To checkout run:
```
git clone --recurse-submodules git@github.com:jonasschultemattler/kmerdict.git
```

### Compile
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -D CMAKE_CXX_COMPILER=g++-14
```

### Run
E.g.
```
./source/main bq -i ../data/salmonella_100_k31_ust.fa.gz -q ../data/SRR5833294.10K.fastq.gz -k 31 -m 15
```
