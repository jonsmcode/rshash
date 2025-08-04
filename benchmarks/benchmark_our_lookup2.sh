#!/bin/bash

PROGRAM="../build/source/lookup2"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/unitigs" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="lookup2-results-$today.csv"

run()
{
  FILES=$(ls $1*.fa.gz)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')
    length=$(python3 -c "import gzip; print(sum(len(line.strip()) for line in gzip.open('$f') if not line.startswith(b'>')))")

    ms=()
    for ((i=0; i<=0; i++)); do
        m=$(echo "l($length)/l(4)+$i" | bc -l)
        m=$(printf "%.0f" "$m")
        ms+=("$m")
    done

    for m in "${ms[@]}"; do
      if [[ "$BASENAME" == *"bacterial"* || "$BASENAME" == *"human"* || "$BASENAME" == *"cod"* || "$BASENAME" == *"kestrel"* ]]; then
          if [ "$k" -le 32 ]; then
            if [ "$m" -le 32 ]; then
              echo $BASENAME
              echo $m

              echo $f >> $LOG

              /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -m $((m - 1)) -n $m > prog_out.txt 2>&1

              cat prog_out.txt >> $LOG

              file_size=$(stat -f%z "${BASENAME}.dict")
              buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
              buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
              textlength=$(grep "text length: " prog_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
              density_r1=$(awk '/density r1/ {print $3}' prog_out.txt | tr -d '%')
              density_r2=$(awk '/density r2/ {print $3}' prog_out.txt | tr -d '%')
              density_s1=$(awk '/density s1/ {print $3}' prog_out.txt | tr -d '%')
              density_s2=$(awk '/density s2/ {print $3}' prog_out.txt | tr -d '%')
              no_kmers=$(awk '/no kmers/ {print $3}' prog_out.txt)
              no_minimiser1=$(awk '/no minimiser1/ {print $3}' prog_out.txt)
              no_minimiser2=$(awk '/no minimiser2/ {print $3}' prog_out.txt)
              no_distinct_minimiser1=$(awk '/no distinct minimiser1/ {print $4}' prog_out.txt)
              no_distinct_minimiser2=$(awk '/no distinct minimiser2/ {print $4}' prog_out.txt)
              spaceoffsets1=$(grep '^offsets1:' prog_out.txt | cut -d':' -f2 | xargs)
              spaceoffsets2=$(grep '^offsets2:' prog_out.txt | cut -d':' -f2 | xargs)
              spacer1=$(grep '^R_1:' prog_out.txt | cut -d':' -f2 | xargs)
              spacer2=$(grep '^R_2:' prog_out.txt | cut -d':' -f2 | xargs)
              spaces1=$(grep '^S_1:' prog_out.txt | cut -d':' -f2 | xargs)
              spaces2=$(grep '^S_2:' prog_out.txt | cut -d':' -f2 | xargs)
              spacetotal=$(grep '^total:' prog_out.txt | cut -d':' -f2 | xargs)

              parent_dir=$(dirname "$f")
              file_name=$(basename "$f")

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
                  /usr/bin/time -l -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query > prog_out.txt 2>&1

                  cat prog_out.txt >> $LOG
                  
                  querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                  querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                  k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  
                  echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spacer1,$spacer2,$spaces1,$spaces2,$density_r1,$density_r2,$density_s1,$density_s2,$no_minimiser1,$no_minimiser2,$no_distinct_minimiser1,$no_distinct_minimiser2,$spacetotal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"
                  # echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$spaceoffsets,$spacer,$spaces,$density_r,$density_s,$no_minimiser,$spacetotal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

                  # rm -f time.txt
                  # rm -f prog_out.txt
              done

              # high hitrates
              if [[ "$BASENAME" == *"bacterial"* ]]; then
                queries=("${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz")
              elif [[ "$BASENAME" == *"cod"* ]]; then
                queries=("${parent_dir/unitigs/queries}/SRR12858649.fastq.gz")
              elif [[ "$BASENAME" == *"kestrel"* ]]; then
                queries=("${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz")
              elif [[ "$BASENAME" == *"human"* ]]; then
                queries=("${parent_dir/unitigs/queries}/SRR5833294.fastq.gz")
              else
                queries=()
              fi

              for query in "${queries[@]}"; do
                  echo $f >> $LOG
                  echo $query >> $LOG
                  echo $query
                  /usr/bin/time -l -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query > prog_out.txt 2>&1

                  cat prog_out.txt >> $LOG
                  
                  querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                  querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                  k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  
                  echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spacer1,$spacer2,$spaces1,$spaces2,$density_r1,$density_r2,$density_s1,$density_s2,$no_minimiser1,$no_minimiser2,$no_distinct_minimiser1,$no_distinct_minimiser2,$spacetotal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"
                  # echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$spaceoffsets,$spacer,$spaces,$density_r,$density_s,$no_minimiser,$spacetotal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

                  # rm -f time.txt
                  # rm -f prog_out.txt
              done
            fi
          fi
      fi
    done


  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],density_r1 [%],density_r2 [%],density_s1 [%],density_s2 [%],no minimizer1,no minimizer2, no distinct minimizer1,no distinct minimizer2,spacetotal [bits/kmer],querytime [ns/kmer],querymem [B],kmers,found" > "$CSV"
  # echo "textfile,queryfile,k,m,buildtime [s],buildmem [B],indexsize [B],spaceoffsets [bits/kmer],spaceR [bits/kmer],spaceS [bits/kmer],density_r [%],density_s [%], no minimizer, spacetotal [bits/kmer], querytime [ns/kmer],querymem [B],kmers,found" > "$CSV"
  run $data/
done
