#!/bin/bash

PROGRAM="../../sshash/build/sshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/unitigs" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="sshash-results-$today.csv"

m=15

run()
{
  FILES=$(ls $1*.fa.gz)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')
    length=$(python3 -c "import gzip; print(sum(len(line.strip()) for line in gzip.open('$f') if not line.startswith(b'>')))")

    ms=()
    for ((i=-1; i<=2; i++)); do
        m=$(echo "l($length)/l(4)+$i" | bc -l)
        m=$(printf "%.0f" "$m")
        ms+=("$m")
    done

    for m in "${ms[@]}"; do
      if [[ "$BASENAME" == *"bacterial"* || "$BASENAME" == *"human"* || "$BASENAME" == *"cod"* || "$BASENAME" == *"kestrel"* ]]; then
        if [ "$k" -le 32 ]; then
          echo $BASENAME

          echo $f >> $LOG

          /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -o "${BASENAME}.index" -k $k -m $m > prog_out.txt 2>&1

          cat prog_out.txt >> $LOG

          file_size=$(stat -f%z "${BASENAME}.index")
          buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
          buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
          num_super_kmers=$(grep "^num_super_kmers" prog_out.txt | awk '{print $2}')
          space=$(grep "total:" prog_out.txt | sed -E 's/^ *total: *([0-9.eE+-]+).*/\1/')
          space_o=$(grep "offsets:" prog_out.txt | sed -E 's/^ *offsets: *([0-9.eE+-]+).*/\1/')
          space_m=$(grep "minimizers:" prog_out.txt | sed -E 's/^ *minimizers: *([0-9.eE+-]+).*/\1/')
          bits_key=$(awk '/minimizers:/ {for(i=1;i<=NF;i++) if($i ~ /\[bits\/key\]/) print $(i-1)}' prog_out.txt | tr -d '(')

          parent_dir=$(dirname "$f")
          file_name=$(basename "$f")

          if [[ "$BASENAME" == *"bacterial"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5833294.K100.fastq.gz")
          elif [[ "$BASENAME" == *"cod"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR11449743_1.K100.fastq.gz")
          elif [[ "$BASENAME" == *"kestrel"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR12858649.K100.fastq.gz")
          elif [[ "$BASENAME" == *"human"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5901135_1.K100.fastq.gz")
          else
            queries=()
          fi

          for query in "${queries[@]}"; do
              echo $f >> $LOG
              echo $query >> $LOG
              echo $query
              /usr/bin/time -l -o time.txt $PROGRAM query -i "${BASENAME}.index" -q $query > prog_out.txt 2>&1

              cat prog_out.txt >> $LOG
              
              querytime=$(cat time.txt | grep "real" | awk '{print $1}')
              querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
              k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
              found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
              # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
              querytimekmer=$(grep 'ns/kmer' prog_out.txt | awk -F'/' '{print $4}' | awk '{print $1}')
              num_extensions=$(grep "num_extensions" prog_out.txt | sed -E 's/.*num_extensions = ([0-9]+).*/\1/')
              
              echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$space_o,$space_m,$bits_key,$num_super_kmers,$space,$querytimekmer,$querymem,$k_mers",$found",$num_extensions >> "$CSV"

          done

          # high hitrates
          if [[ "$BASENAME" == *"bacterial"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5901135_1.K100.fastq.gz")
          elif [[ "$BASENAME" == *"cod"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR12858649.K100.fastq.gz")
          elif [[ "$BASENAME" == *"kestrel"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR11449743_1.K100.fastq.gz")
          elif [[ "$BASENAME" == *"human"* ]]; then
            queries=("${parent_dir/unitigs/queries}/SRR5833294.K100.fastq.gz")
          else
            queries=()
          fi

          for query in "${queries[@]}"; do
              echo $f >> $LOG
              echo $query >> $LOG
              echo $query
              /usr/bin/time -l -o time.txt $PROGRAM query -i "${BASENAME}.index" -q $query > prog_out.txt 2>&1

              cat prog_out.txt >> $LOG
              
              querytime=$(cat time.txt | grep "real" | awk '{print $1}')
              querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
              k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
              found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
              # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
              querytimekmer=$(grep 'ns/kmer' prog_out.txt | awk -F'/' '{print $4}' | awk '{print $1}')
              num_extensions=$(grep "num_extensions" prog_out.txt | sed -E 's/.*num_extensions = ([0-9]+).*/\1/')
              
              echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$space_o,$space_m,$bits_key,$num_super_kmers,$space,$querytimekmer,$querymem,$k_mers",$found,$num_extensions" >> "$CSV"

          done


        fi
      fi
    done

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,buildtime[s],buildmem[B],indexsize[B], space offsets[bits/kmer], space minimizers[bits/kmer], bits/key, no super kmers, space[bits/kmer], querytime[ns/kmer],querymem[B],kmers,found,extensions" > "$CSV"
  run $data/
done
