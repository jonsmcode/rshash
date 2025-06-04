#!/bin/bash

PROGRAM="../../sshash/build/sshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../test/datasets/synthetic" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="sshash-results-server-synthetic$today.csv"
compression=9

run()
{
  FILES=$(ls ${1}*text.fasta | awk '{ print length, $0 }' | sort -n | cut -d" " -f2-)

  for f in $FILES
  do
    BASENAME=$(basename "$f" .fasta)
    # k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')
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

      k=$((2 * compression + m - 2))
      k=$(( k > 31 ? 31 : k ))

      /usr/bin/time -v -o time.txt $PROGRAM build -i $f -o "${DIR}/${BASENAME}.index" -k $k -m $m > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG

      file_size=$(stat -c%s "${DIR}/${BASENAME}.index")
      buildtime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
      buildmem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
      num_super_kmers=$(grep "^num_super_kmers" prog_out.txt | awk '{print $2}')
      space=$(grep "total:" prog_out.txt | sed -E 's/^ *total: *([0-9.eE+-]+).*/\1/')
      space_o=$(grep "offsets:" prog_out.txt | sed -E 's/^ *offsets: *([0-9.eE+-]+).*/\1/')
      space_m=$(grep "minimizers:" prog_out.txt | sed -E 's/^ *minimizers: *([0-9.eE+-]+).*/\1/')
      bits_key=$(awk '/minimizers:/ {for(i=1;i<=NF;i++) if($i ~ /\[bits\/key\]/) print $(i-1)}' prog_out.txt | tr -d '(')

      parent_dir=$(dirname "$f")
      file_name=$(basename "$f")

      basename2=$(echo "$f" | sed 's/\_text.fasta$//')
      query=(${basename2}_query.fasta)

      echo $f >> $LOG 
      echo $query >> $LOG
      echo $query

      /usr/bin/time -l -o time.txt $PROGRAM query -i "${DIR}/${BASENAME}.index" -q $query > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG
      
      k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
      found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
      querytime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
      querymem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
      # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
      
      echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$space,$space_o,$space_m,$bits_key,$num_super_kmers,$querytime,$querymem,$k_mers",$found" >> "$CSV"
    done

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,buildtime[s],buildmem[B],indexsize[B], space[bits/kmer], space offsets[bits/kmer], space minimizers[bits/kmer], bits/key, no super kmers, querytime [ns/kmer],querymem[B],kmers,found" > "$CSV"
  run $data/
done
