#!/bin/bash

PROGRAM="../build/source/lookup3compsimd"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets" >/dev/null 2>&1 && pwd )"
LOG="log3c.txt"
CSV="comp_lookup3compsimd-results-$today.csv"

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

        minimisers=( $((m + 2)) $((m + 4)) $((m + 4))  $((m + 3)) $((m + 4)) $((m + 4))  $((m + 4)) $((m + 4)) $((m + 4)) )
        thresholds=( 8 16 64  16 32 128  32 64 512 )
        spans=( $((k - m + 1)) )

        for ((i=0; i<${#minimisers[@]}; i+=3)); do
          m1=${minimisers[i]}
          m2=${minimisers[i+1]}
          m3=${minimisers[i+2]}

          for ((j=0; j<${#thresholds[@]}; j+=3)); do
            for span in "${spans[@]}"; do

              t1=${thresholds[j]}
              t2=${thresholds[j+1]}
              t3=${thresholds[j+2]}
              echo $BASENAME 
              echo $m1 $t1 $m2 $t2 $m3 $t3

              echo $f >> $LOG
              echo $m1 $m2 $m3 $t1 $t2 $t3 >> $LOG

              /usr/bin/time -v -o ctime3.txt $PROGRAM build -i "$f" -d "${BASENAME}.complookup3.dict" -k $k --m1 $m1 --t1 $t1 --m2 $m2 --t2 $t2 --m3 $m3 --t3 $t3 -s $span > comp3_out.txt 2>&1

              cat comp3_out.txt >> $LOG

              file_size=$(stat -c%s "${BASENAME}.complookup3.dict")
              buildtime=$(grep "User time" ctime3.txt | awk -F': ' '{print $2}')
              buildmem=$(grep "Maximum resident set size" ctime3.txt | awk -F': ' '{print $2}')
              textlength=$(grep "text length: " comp3_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
              density_r1=$(awk '/density r1/ {print $3}' comp3_out.txt | tr -d '%')
              density_r2=$(awk '/density r2/ {print $3}' comp3_out.txt | tr -d '%')
              density_r3=$(awk '/density r3/ {print $3}' comp3_out.txt | tr -d '%')
              density_s1=$(awk '/density s1/ {print $3}' comp3_out.txt | tr -d '%')
              density_s2=$(awk '/density s2/ {print $3}' comp3_out.txt | tr -d '%')
              density_s3=$(awk '/density s3/ {print $3}' comp3_out.txt | tr -d '%')
              density_ht=$(grep "no kmers HT:" comp3_out.txt | awk '{print $(NF)}' | sed 's/%//')
              no_kmers=$(awk '/text kmers/ {print $3}' comp3_out.txt)
              no_minimiser1=$(awk '/no minimiser1/ {print $3}' comp3_out.txt)
              no_minimiser2=$(awk '/no minimiser2/ {print $3}' comp3_out.txt)
              no_minimiser3=$(awk '/no minimiser3/ {print $3}' comp3_out.txt)
              no_distinct_minimiser1=$(awk '/no distinct minimiser1/ {print $4}' comp3_out.txt)
              no_distinct_minimiser2=$(awk '/no distinct minimiser2/ {print $4}' comp3_out.txt)
              no_distinct_minimiser3=$(awk '/no distinct minimiser3/ {print $4}' comp3_out.txt)
              spaceoffsets1=$(grep '^offsets1:' comp3_out.txt | cut -d':' -f2 | xargs)
              spaceoffsets2=$(grep '^offsets2:' comp3_out.txt | cut -d':' -f2 | xargs)
              spaceoffsets3=$(grep '^offsets3:' comp3_out.txt | cut -d':' -f2 | xargs)
              spaceht=$(grep '^Hashtable:' comp3_out.txt | cut -d':' -f2 | xargs)
              spacer1=$(grep '^R_1:' comp3_out.txt | cut -d':' -f2 | xargs)
              spacer2=$(grep '^R_2:' comp3_out.txt | cut -d':' -f2 | xargs)
              spacer3=$(grep '^R_3:' comp3_out.txt | cut -d':' -f2 | xargs)
              spaces1=$(grep '^S_1:' comp3_out.txt | cut -d':' -f2 | xargs)
              spaces2=$(grep '^S_2:' comp3_out.txt | cut -d':' -f2 | xargs)
              spaces3=$(grep '^S_3:' comp3_out.txt | cut -d':' -f2 | xargs)
              spacetotal=$(grep '^total:' comp3_out.txt | cut -d':' -f2 | xargs)
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
                  /usr/bin/time -v -o ctime3.txt $PROGRAM query -d "${BASENAME}.complookup3.dict" -q $query > comp3_out.txt 2>&1

                  cat comp3_out.txt >> $LOG
                  
                  querytime=$(grep "User time" ctime3.txt | awk -F': ' '{print $2}')
                  querymem=$(grep "Maximum resident set size" ctime3.txt | awk -F': ' '{print $2}')
                  k_mers=$(grep "num_kmers" comp3_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" comp3_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  querytimekmer=$(grep 'time_per_kmer' comp3_out.txt | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
                  extensions=$(grep 'extensions' comp3_out.txt | awk -F'=' '{print $2}' | awk '{print $1}')
                  
                  echo "$f,$query,$k,$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spaceht,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$extensions,$k_mers",$found" >> "$CSV"

                  # rm -f ctime3.txt
                  # rm -f comp3_out.txt
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
                  /usr/bin/time -v -o ctime3.txt $PROGRAM query -d "${BASENAME}.complookup3.dict" -q $query > comp3_out.txt 2>&1

                  cat comp3_out.txt >> $LOG
                  
                  querytime=$(grep "User time" ctime3.txt | awk -F': ' '{print $2}')
                  querymem=$(grep "Maximum resident set size" ctime3.txt | awk -F': ' '{print $2}')
                  k_mers=$(grep "num_kmers" comp3_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                  found=$(grep "num_positive_kmers" comp3_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                  # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                  querytimekmer=$(grep 'time_per_kmer' comp3_out.txt | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
                  extensions=$(grep 'extensions' comp3_out.txt | awk -F'=' '{print $2}' | awk '{print $1}')
                  
                  echo "$f,$query,$k,$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spaceht,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$extensions,$k_mers",$found" >> "$CSV"

                  # rm -f ctime3.txt
                  # rm -f comp3_out.txt
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
  echo "textfile,queryfile,k,m1,m2,m3,t1,t2,t3,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceoffsets3 [bits/kmer],spaceHT [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceR3 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],spaceS3 [bits/kmer],density_r1 [%],density_r2 [%],density_r3 [%],density_s1 [%],density_s2 [%],density_s3 [%],kmers HT [%],no minimizer1,no minimizer2,no minimizer3, no distinct minimizer1,no distinct minimizer2,no distinct minimizer3,space theo [bits/kmer],space real [bits/kmer],querytime [ns/kmer],querymem [B],extensions,kmers,found" > "$CSV"
  run $data/
done
