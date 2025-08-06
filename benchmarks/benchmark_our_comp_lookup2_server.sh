#!/bin/bash

PROGRAM="../build/source/lookup2"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="comp_lookup2-results-$today.csv"

run()
{
  FILES=$(ls $1*.fa.gz)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')
    if [[ "$BASENAME" == *"bacterial"* || "$BASENAME" == *"human"* || "$BASENAME" == *"cod"* || "$BASENAME" == *"kestrel"* ]]; then
      if [ "$k" -le 32 ]; then
        length=$(python3 -c "import gzip; print(sum(len(line.strip()) for line in gzip.open('$f') if not line.startswith(b'>')))")

        m=$(echo "l($length)/l(4)" | bc -l)
        m=$(printf "%.0f" "$m")

        thresholds=(50 20 10)
        
        for ((i=-1; i<=0; i++)); do
          m1=($((m + i)))
          for ((j=1; j<=2; j++)); do
            m2=($((m + i + j)))
            for thres in "${thresholds[@]}"; do
              echo $BASENAME
              echo $m1
              echo $m2
              echo $thres

              echo $f >> $LOG

              /usr/bin/time -v -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -m $m1 -n $m2 -t $thres -c > prog_out.txt 2>&1

              cat prog_out.txt >> $LOG

              file_size=$(stat -c%s "${BASENAME}.dict")
              buildtime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
              buildmem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
              textlength=$(grep "text length: " prog_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
              density_r1=$(awk '/density r1/ {print $3}' prog_out.txt | tr -d '%')
              density_r2=$(awk '/density r2/ {print $3}' prog_out.txt | tr -d '%')
              density_s1=$(awk '/density s1/ {print $3}' prog_out.txt | tr -d '%')
              density_s2=$(awk '/density s2/ {print $3}' prog_out.txt | tr -d '%')
              density_ht=$(grep "no kmers HT:" prog_out.txt | awk '{print $(NF)}' | sed 's/%//')
              no_kmers=$(awk '/text kmers/ {print $3}' prog_out.txt)
              no_minimiser1=$(awk '/no minimiser1/ {print $3}' prog_out.txt)
              no_minimiser2=$(awk '/no minimiser2/ {print $3}' prog_out.txt)
              no_distinct_minimiser1=$(awk '/no distinct minimiser1/ {print $4}' prog_out.txt)
              no_distinct_minimiser2=$(awk '/no distinct minimiser2/ {print $4}' prog_out.txt)
              spaceoffsets1=$(grep '^offsets1:' prog_out.txt | cut -d':' -f2 | xargs)
              spaceoffsets2=$(grep '^offsets2:' prog_out.txt | cut -d':' -f2 | xargs)
              spaceht=$(grep '^Hashtable:' prog_out.txt | cut -d':' -f2 | xargs)
              spacer1=$(grep '^R_1:' prog_out.txt | cut -d':' -f2 | xargs)
              spacer2=$(grep '^R_2:' prog_out.txt | cut -d':' -f2 | xargs)
              spaces1=$(grep '^S_1:' prog_out.txt | cut -d':' -f2 | xargs)
              spaces2=$(grep '^S_2:' prog_out.txt | cut -d':' -f2 | xargs)
              spacetotal=$(grep '^total:' prog_out.txt | cut -d':' -f2 | xargs)
              spacetotalreal=$(echo "scale=2; $file_size / $no_kmers * 8" | bc)

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
                  /usr/bin/time -v -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query -c > prog_out.txt 2>&1

                  cat prog_out.txt >> $LOG
                  
                  querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                  querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                  k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  
                  echo "$f,$query,$k,$m1,$m2,$thres,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceht,$spacer1,$spacer2,$spaces1,$spaces2,$density_r1,$density_r2,$density_s1,$density_s2,$density_ht,$no_minimiser1,$no_minimiser2,$no_distinct_minimiser1,$no_distinct_minimiser2,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

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
                  /usr/bin/time -v -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query -c > prog_out.txt 2>&1

                  cat prog_out.txt >> $LOG
                  
                  querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                  querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                  k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  
                  echo "$f,$query,$k,$m1,$m2,$thres,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceht,$spacer1,$spacer2,$spaces1,$spaces2,$density_r1,$density_r2,$density_s1,$density_s2,$density_ht,$no_minimiser1,$no_minimiser2,$no_distinct_minimiser1,$no_distinct_minimiser2,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

                  # rm -f time.txt
                  # rm -f prog_out.txt
              done
              
            done
          done
        done

      fi
    fi

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m1,m2,thres,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceHT [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],density_r1 [%],density_r2 [%],density_s1 [%],density_s2 [%],kmers HT [%],no minimizer1,no minimizer2, no distinct minimizer1,no distinct minimizer2,space theo [bits/kmer],space real [bits/kmer],querytime [ns/kmer],querymem [B],kmers,found" > "$CSV"
  run $data/
done
