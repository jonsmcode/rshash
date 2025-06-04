#!/bin/bash

PROGRAM="../build/source/lookup_unitigs"
# PROGRAM="../build/source/comp_lookup"

today=$(date +%Y-%m-%d-%H-%M-%S)

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/synthetic" >/dev/null 2>&1 && pwd )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../test/datasets/synthetic" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
# CSV="comp_lookup-results-synthetic$today.csv"
CSV="lookup-results-synthetic$today.csv"
# CSV="sshash-results-synthetic$today.csv"
compression=9

run()
{
  FILES=$(ls ${1}*text.fasta | awk '{ print length, $0 }' | sort -n | cut -d" " -f2-)

  for f in $FILES
  do
    BASENAME=$(basename "$f" .fasta)
    length=$(awk -F'_n|_o' '{print $2}' <<< "$BASENAME")

    ms=()
    for ((i=-2; i<=2; i++)); do
        m=$(echo "l($length)/l(4)+$i" | bc -l)
        m=$(printf "%.0f" "$m")
        ms+=("$m")
    done

    for m in "${ms[@]}"; do
      echo $BASENAME
      echo $m
      echo $f >> $LOG

      # k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')
      k=$((2 * compression + m - 2))
      k=$(( k > 31 ? 31 : k ))

      # /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -d "${DIR}/${BASENAME}.dict" -k $k -m $m > prog_out.txt 2>&1
      /usr/bin/time -l -o time.txt $PROGRAM build -i $f -d "${DIR}/${BASENAME}.dict" -k $k -m $m > prog_out.txt 2>&1
      # /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -o "${DIR}/${BASENAME}.index" -k $k -m $m > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG

      file_size=$(stat -f%z "${DIR}/${BASENAME}.dict")
      buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
      buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
      textlength=$(grep "text length: " prog_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
      no_distinct_kmers=$(awk '/no distinct kmers/ {print $4}' prog_out.txt | tr -d ',')
      density_r=$(awk '/density r/ {print $3}' prog_out.txt | tr -d '%')
      density_s=$(awk '/density s/ {print $3}' prog_out.txt | tr -d '%')
      no_kmers=$(awk '/no kmers/ {print $3}' prog_out.txt)
      no_minimiser=$(awk '/no minimiser/ {print $3}' prog_out.txt)
      freq_kmers=$(awk '/freq kmers/ {print $3}' prog_out.txt)
      spaceoffsets=$(grep '^offsets:' prog_out.txt | cut -d':' -f2 | xargs)
      spacer=$(grep '^R:' prog_out.txt | cut -d':' -f2 | xargs)
      spaces=$(grep '^S:' prog_out.txt | cut -d':' -f2 | xargs)
      spacetotal=$(grep '^total:' prog_out.txt | cut -d':' -f2 | xargs)

      parent_dir=$(dirname "$f")
      file_name=$(basename "$f")

      basename2=$(echo "$f" | sed 's/\_text.fasta$//')
      query=(${basename2}_query.fasta)

      echo $f >> $LOG
      echo $query >> $LOG
      echo $query
      /usr/bin/time -l -o time.txt $PROGRAM query -d "${DIR}/${BASENAME}.dict" -q $query > prog_out.txt 2>&1
      # /usr/bin/time -l -o time.txt $PROGRAM query -i "${DIR}/${BASENAME}.index" -q $query > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG
      
      querytime=$(cat time.txt | grep "real" | awk '{print $1}')
      querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
      k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
      found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')

      querytimekmer=$(echo "scale=10; $querytime * 1e9 / $k_mers" | bc)
      
      # echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$querytime,$querymem,$k_mers",$found" >> "$CSV"
      echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$spaceoffsets,$spacer,$spaces,$spacetotal,$density_r,$density_s,$no_minimiser,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"
    done

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,buildtime [s],buildmem [B],indexsize [B],spaceoffsets [bits/kmer],spaceR [bits/kmer],spaceS [bits/kmer],spacetotal [bits/kmer],density_r [%],density_s [%], no minimizer, querytime [s],querymem [B],kmers,found" > "$CSV"
  # echo "textfile,queryfile,k,m,buildtime,buildmem,filesize,querytime,querymem,kmers,found" > "$CSV"
  run $data/
done
