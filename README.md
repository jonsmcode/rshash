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

cmake -DCMAKE_BUILD_TYPE=Release .. -D DCMAKE_CXX_COMPILER=$(brew --prefix llvm@19)/bin/clang++
```

### Datasets

```
mkdir datasets
cd datasets
wget https://zenodo.org/records/7239205/files/celegans.k31.unitigs.fa.ust.fa.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR162/082/SRR16288382/SRR16288382_1.fastq.gz
```


### Run
```
.build/source/lookup build -i ../datasets/celegans.k31.unitigs.fa.ust.fa.gz -d ../datasets/celegans.k31.unitigs.fa.ust.dict -k 31 -m 15
.build/source/lookup query -i ../datasets/celegans.k31.unitigs.fa.ust.fa.gz -d ../datasets/celegans.k31.unitigs.fa.ust.dict -q ../datasets/SRR16288382_1.fastq.gz -k 31 -m 15
```
