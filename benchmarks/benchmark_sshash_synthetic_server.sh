#!/bin/bash

PROGRAM="../../sshash/build/sshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/synthetic" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="sshash-results-server-synthetic$today.csv"


run()
{
  FILES=$(ls ${1}*text.fasta | awk '{ print length, $0 }' | sort -n | cut -d" " -f2-)

  for f in $FILES
  do
    BASENAME=$(basename "$f" .fasta)
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')
    length=$(awk -F'_n|_o' '{print $2}' <<< "$BASENAME")

    ms=()
    for ((i=-2; i<=2; i++)); do
        m=$(echo "l($length)/l(4)+$i" | bc -l)
        m=$(printf "%.0f" "$m")
        ms+=("$m")
    done

    for m in "${ms[@]}"; do
      echo $m
      echo $f >> $LOG

      /usr/bin/time -v -o time.txt $PROGRAM build -i $f -o "${DIR}/${BASENAME}.index" -k $k -m $m > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG

      file_size=$(stat -c%s "${DIR}/${BASENAME}.index")
      buildtime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
      buildmem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
      space=$(grep "total:" prog_out.txt | sed -E 's/^ *total: *([0-9.eE+-]+).*/\1/')

      parent_dir=$(dirname "$f")
      file_name=$(basename "$f")

      basename2=$(echo "$f" | sed 's/\_text.fasta$//')
      query=(${basename2}_query.fasta)

      echo $f >> $LOG
      echo $query >> $LOG
      echo $query

      /usr/bin/time -l -o time.txt $PROGRAM query -i "${DIR}/${BASENAME}.index" -q $query > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG
      
      querytime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
      querymem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
      k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
      found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
      
      echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$space,$querytime,$querymem,$k_mers",$found" >> "$CSV"
    done

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,buildtime [s],buildmem [MB],indexsize [MB], space [bits\kmer],querytime [s],querymem [MB],kmers,found" > "$CSV"
  run $data/
done
