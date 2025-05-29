#!/bin/bash

PROGRAM="../../sshash/build/sshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/synthetic" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="sshash-results-synthetic$today.csv"


run()
{
  FILES=$(ls ${1}*text.fasta | awk '{ print length, $0 }' | sort -n | cut -d" " -f2-)

  for f in $FILES
  do
    BASENAME=$(basename "$f" .fasta)
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')
    length=$(awk -F'_n|_o' '{print $2}' <<< "$BASENAME")

    ms=()
    for i in {1..2}; do
        m=$(echo "l($length)/l(4)+$i" | bc -l)
        m=$(printf "%.0f" "$m")
        ms+=("$m")
    done

    for m in "${ms[@]}"; do
      echo $m
      echo $f >> $LOG

      # /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -o "${DIR}/${BASENAME}.index" -k $k -m $m > prog_out.txt 2>&1
      /usr/bin/time -l -o time.txt $PROGRAM build -i $f -o "${DIR}/${BASENAME}.index" -k $k -m $m > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG

      file_size=$(stat -f%z "${DIR}/${BASENAME}.index")
      buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
      buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')

      parent_dir=$(dirname "$f")
      file_name=$(basename "$f")

      basename2=$(echo "$f" | sed 's/\_text.fasta$//')
      query=(${basename2}_query.fasta)
      # query=(${basename2}_query.fastq.gz)

      echo $f >> $LOG
      echo $query >> $LOG
      echo $query

      /usr/bin/time -l -o time.txt $PROGRAM query -i "${DIR}/${BASENAME}.index" -q $query > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG
      
      querytime=$(cat time.txt | grep "real" | awk '{print $1}')
      querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
      space=$(grep "total:" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
      k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
      found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
      
      echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$querytime,$querymem,$k_mers",$found" >> "$CSV"
    done

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,buildtime,buildmem,indexsize,querytime,querymem,kmers,found" > "$CSV"
  run $data/
done
