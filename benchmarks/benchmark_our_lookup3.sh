#!/bin/bash

PROGRAM="../build/source/lookup3"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/unitigs" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="comp_lookup3-results-$today.csv"

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

        t1=20
        t2=50
        t3=1000

        m1s=()
        for ((i=-1; i>=-1; i--)); do
            m1s+=($((m + i)))
        done
        m2s=()
        for ((i=1; i<=1; i++)); do
            m2s+=($((m + i)))
        done
        m3s=()
        for ((i=2; i<=2; i++)); do
            m3s+=($((m + i)))
        done
        for m1 in "${m1s[@]}"; do
          for m2 in "${m2s[@]}"; do
            for m3 in "${m3s[@]}"; do
              echo $BASENAME
              echo $m1
              echo $m2
              echo $m3

              echo $f >> $LOG

              /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -c --m1 $m1 --t1 $t1 --m2 $m2 --t2 $t2 --m3 $m3 --t3 $t3 > prog_out.txt 2>&1

              cat prog_out.txt >> $LOG

              file_size=$(stat -f%z "${BASENAME}.dict")
              buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
              buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
              textlength=$(grep "text length: " prog_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
              density_r1=$(awk '/density r1/ {print $3}' prog_out.txt | tr -d '%')
              density_r2=$(awk '/density r2/ {print $3}' prog_out.txt | tr -d '%')
              density_r3=$(awk '/density r3/ {print $3}' prog_out.txt | tr -d '%')
              density_s1=$(awk '/density s1/ {print $3}' prog_out.txt | tr -d '%')
              density_s2=$(awk '/density s2/ {print $3}' prog_out.txt | tr -d '%')
              density_s3=$(awk '/density s3/ {print $3}' prog_out.txt | tr -d '%')
              density_ht=$(grep "no kmers HT:" prog_out.txt | awk '{print $(NF)}' | sed 's/%//')
              no_kmers=$(awk '/text kmers/ {print $3}' prog_out.txt)
              no_minimiser1=$(awk '/no minimiser1/ {print $3}' prog_out.txt)
              no_minimiser2=$(awk '/no minimiser2/ {print $3}' prog_out.txt)
              no_minimiser3=$(awk '/no minimiser3/ {print $3}' prog_out.txt)
              no_distinct_minimiser1=$(awk '/no distinct minimiser1/ {print $4}' prog_out.txt)
              no_distinct_minimiser2=$(awk '/no distinct minimiser2/ {print $4}' prog_out.txt)
              no_distinct_minimiser3=$(awk '/no distinct minimiser3/ {print $4}' prog_out.txt)
              spaceoffsets1=$(grep '^offsets1:' prog_out.txt | cut -d':' -f2 | xargs)
              spaceoffsets2=$(grep '^offsets2:' prog_out.txt | cut -d':' -f2 | xargs)
              spaceoffsets3=$(grep '^offsets3:' prog_out.txt | cut -d':' -f2 | xargs)
              spaceht=$(grep '^Hashtable:' prog_out.txt | cut -d':' -f2 | xargs)
              spacer1=$(grep '^R_1:' prog_out.txt | cut -d':' -f2 | xargs)
              spacer2=$(grep '^R_2:' prog_out.txt | cut -d':' -f2 | xargs)
              spacer3=$(grep '^R_3:' prog_out.txt | cut -d':' -f2 | xargs)
              spaces1=$(grep '^S_1:' prog_out.txt | cut -d':' -f2 | xargs)
              spaces2=$(grep '^S_2:' prog_out.txt | cut -d':' -f2 | xargs)
              spaces3=$(grep '^S_3:' prog_out.txt | cut -d':' -f2 | xargs)
              spacetotal=$(grep '^total:' prog_out.txt | cut -d':' -f2 | xargs)
              spacetotalreal=$(echo "scale=2; $file_size/$no_kmers*8" | bc)

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
                  /usr/bin/time -l -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query  -c > prog_out.txt 2>&1

                  cat prog_out.txt >> $LOG
                  
                  querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                  querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                  k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  querytimekmer=$(grep 'ns_per_kmer' prog_out.txt | awk -F'=' '{print $2}' | awk '{print $1}')
                  
                  echo "$f,$query,$k,$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spaceht,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

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
                  /usr/bin/time -l -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query  -c > prog_out.txt 2>&1

                  cat prog_out.txt >> $LOG
                  
                  querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                  querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                  k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  querytimekmer=$(grep 'ns_per_kmer' prog_out.txt | awk -F'=' '{print $2}' | awk '{print $1}')
                  
                  echo "$f,$query,$k,$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spaceht,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

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
  echo "textfile,queryfile,k,m1,m2,m3,t1,t2,t3,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceoffsets3 [bits/kmer],spaceHT [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceR3 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],spaceS3 [bits/kmer],density_r1 [%],density_r2 [%],density_r3 [%],density_s1 [%],density_s2 [%],density_s3 [%],kmers HT [%],no minimizer1,no minimizer2,no minimizer3, no distinct minimizer1,no distinct minimizer2,no distinct minimizer3,space theo [bits/kmer],space real [bits/kmer],querytime [ns/kmer],querymem [B],kmers,found" > "$CSV"
  run $data/
done
