#!/bin/bash

PROGRAM="../build/source/lookup"

today=$(date +%Y-%m-%d-%H-%M-%S)

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../test/datasets" >/dev/null 2>&1 && pwd )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/unitigs" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="ourspace-results-$today.csv"

m=15

run()
{
  FILES=$(ls $1*.fa.gz)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')

    if [[ "$BASENAME" == *"bacterial"* || "$BASENAME" == *"human"* || "$BASENAME" == *"cod"* || "$BASENAME" == *"kestrel"* ]]; then
        if [ "$k" -le 32 ]; then
          echo $BASENAME

          echo $f >> $LOG

          /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -m $m > prog_out.txt 2>&1

          cat prog_out.txt >> $LOG

          textlength=$(grep "text length: " prog_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
          # no_distinct_kmers=$(awk '/no distinct kmers/ {print $4}' prog_out.txt | tr -d ',')
          no_minimiser=$(awk '/no minimiser/ {print $3}' prog_out.txt)
          # no_kmers=$(awk '/no kmers/ {print $3}' prog_out.txt)

          buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
          buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
          file_size=$(stat -f%z "${BASENAME}.dict")
          
          # density_r=$(awk '/density r/ {print $3}' prog_out.txt | tr -d '%')
          # freq_kmers=$(awk '/freq kmers/ {print $3}' prog_out.txt)
          # density_s=$(awk '/density s/ {print $3}' prog_out.txt | tr -d '%')
          # offset_width=$(awk '/offset width/ {print $3}' prog_out.txt)
          # span_width=$(awk '/span width/ {print $3}' prog_out.txt)
          # bytes_o=$(awk '/allocated/ {print $2}' prog_out.txt)

          parent_dir=$(dirname "$f")
          file_name=$(basename "$f")

          echo "$f,$textlength,$no_minimiser,$buildtime,$buildmem",$file_size >> "$CSV"
        fi
    fi
  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,textlength,minimiser,buildtime,buildmem,filesize" > "$CSV"
  run $data/
done
