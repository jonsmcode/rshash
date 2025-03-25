#!/bin/bash

PROGRAM="../build/source/locate"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/unitigs" >/dev/null 2>&1 && pwd )"
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/datasets" >/dev/null 2>&1 && pwd )"

m=15

run()
{
  FILES=$(ls $1*.fa.gz)
  # FILES=$(ls $1*_text.fasta)

  for f in $FILES
  do
    # BASENAME=$(echo "$f" | sed 's/_text\.fasta$//')
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')

    if [[ "$BASENAME" == *"human"* ]]; then
        if [ "$k" -le 32 ]; then
          echo $BASENAME

          # $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -m $m

          parent_dir=$(dirname "$f")
          file_name=$(basename "$f")
          
          query="${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz"

          # $PROGRAM query -i "$f" -d "${BASENAME}.dict" -q $query -k $k -m $m > "${BASENAME}.positions" 2>&1
          ../build/test/verify -i $f -q $query -p "${BASENAME}.positions" -k $k > correct.txt

          cat correct.txt

        fi
    fi
  done
}


# for data in $(find $DIR -mindepth 1 -maxdepth 1 -type d); do
for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  run $data/
done
