#!/bin/bash

PROGRAM="../build/source/lookup"

today=$(date +%Y-%m-%d-%H-%M-%S)

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../test/datasets" >/dev/null 2>&1 && pwd )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/unitigs" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="ourlookup-results-$today.csv"

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

          file_size=$(stat -f%z "${BASENAME}.dict")
          buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
          buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
          textlength=$(grep "text length: " prog_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
          # no_distinct_kmers=$(awk '/no distinct kmers/ {print $4}' prog_out.txt | tr -d ',')
          # density_r=$(awk '/density r/ {print $3}' prog_out.txt | tr -d '%')
          # no_kmers=$(awk '/no kmers/ {print $3}' prog_out.txt)
          # no_minimiser=$(awk '/no minimiser/ {print $3}' prog_out.txt)
          # freq_kmers=$(awk '/freq kmers/ {print $3}' prog_out.txt)
          # density_s=$(awk '/density s/ {print $3}' prog_out.txt | tr -d '%')
          # bytes_o=$(awk '/allocated/ {print $2}' prog_out.txt)

          parent_dir=$(dirname "$f")
          file_name=$(basename "$f")

          # if [[ "$BASENAME" == *"bacterial"* ]]; then
          #   queries=("${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz" "${parent_dir/unitigs/queries}/SRR5833294.fastq.gz")
          # elif [[ "$BASENAME" == *"cod"* ]]; then
          #   queries=("${parent_dir/unitigs/queries}/SRR12858649.fastq.gz" "${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz")
          # elif [[ "$BASENAME" == *"kestrel"* ]]; then
          #   queries=("${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz" "${parent_dir/unitigs/queries}/SRR12858649.fastq.gz")
          # elif [[ "$BASENAME" == *"human"* ]]; then
          #   queries=("${parent_dir/unitigs/queries}/SRR5833294.fastq.gz" "${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz")
          # else
          #   queries=()
          # fi
          # low hitrates
          if [[ "$BASENAME" == *"bacterial"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5833294.fastq.gz")
          elif [[ "$BASENAME" == *"cod"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz")
          elif [[ "$BASENAME" == *"kestrel"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR12858649.fastq.gz")
          elif [[ "$BASENAME" == *"human"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz")
          else
            queries=()
          fi

          for query in "${queries[@]}"; do
              echo $f >> $LOG
              echo $query >> $LOG
              echo $query
              /usr/bin/time -l -o time.txt $PROGRAM query -i "$f" -d "${BASENAME}.dict" -q $query -k $k -m $m > prog_out.txt 2>&1

              cat prog_out.txt >> $LOG
              
              querytime=$(cat time.txt | grep "real" | awk '{print $1}')
              querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
              k_mers=$(grep "k-mers" prog_out.txt | sed -E 's/.*k-mers: ([0-9]+).*/\1/')
              found=$(grep "found" prog_out.txt | sed -E 's/.*found: ([0-9]+).*/\1/')
              
              # echo "$f,$query,$buildtime,$buildmem",$file_size,$density_r,$density_s,$freq_kmers,$bytes_o,$querytime,$querymem,$k_mers",$found" >> "$CSV"
              totalspace=
              ((totalspace = file_size + textlength / 4))
              echo "$f,$query,$buildtime,$buildmem",$file_size,$totalspace,$querytime,$querymem,$k_mers",$found" >> "$CSV"

              # rm -f time.txt
              # rm -f prog_out.txt
          done
        fi
    fi
  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  # echo "textfile,queryfile,buildtime,buildmem,space,density_r,density_s,freq_kmers,bytes_o,querytime,querymem,kmers,found" > "$CSV"
  echo "textfile,queryfile,buildtime,buildmem,filesize,totalspace,querytime,querymem,kmers,found" > "$CSV"
  run $data/
done
