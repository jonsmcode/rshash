# RSHash


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

### Datasets

```
mkdir datasets
cd datasets
wget https://zenodo.org/records/17582116/data.fa.gz
```


### Run
```
.build/source/lookup build -i ../datasets/cod.k31.unitigs.fa.ust.fa.gz -d ../datasets/cod.k31.unitigs.fa.ust.dict -k 31 --m1 18 --t1 64 -l 1
.build/source/lookup query -d ../datasets/cod.k31.unitigs.fa.ust.dict -l 1 -q ../datasets/SRR16288382_1.fastq.gz
```
