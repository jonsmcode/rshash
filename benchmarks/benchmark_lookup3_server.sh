#!/bin/bash

RSHASH="../build/source/lookup3compsimd"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="lookup-results-$today.csv"

run()
{
  FILES=$(ls $1*.fa.gz)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')

    # if [[  "$BASENAME" == *"human"* || "$BASENAME" == *"cod"* || "$BASENAME" == *"kestrel"* || "$BASENAME" == *"bacterial"* ]]; then
    if [[  "$BASENAME" == *"human"* || "$BASENAME" == *"cod"* || "$BASENAME" == *"kestrel"*  ]]; then
      if [ "$k" -le 32 ]; then
        length=$(python3 -c "import gzip; print(sum(len(line.strip()) for line in gzip.open('$f') if not line.startswith(b'>')))")

        m=$(echo "l($length)/l(4)" | bc -l)
        m=$(printf "%.0f" "$m")
        m=$((m + 1))

        # params=( $((m)) $((m + 2)) $((m + 3)) 64 128 256  $((m + 1)) $((m + 3)) $((m + 3)) 16 32 64 )
        # params=( $((m)) $((m + 2)) $((m + 3)) 64 64 64  $((m + 1)) $((m + 3)) $((m + 3)) 64 64 64 )
        params=( $((m + 3)) $((m + 3)) $((m + 3)) 64 64 64 )
        span=$((k - m - 2))
        
        for ((i=0; i<${#params[@]}; i+=6)); do
          m1=${params[i]}
          m2=${params[i+1]}
          m3=${params[i+2]}
          t1=${params[i+3]}
          t2=${params[i+4]}
          t3=${params[i+5]}

          echo $f >> $LOG
          echo $f
          echo $m1 $m2 $m3 $t1 $t2 $t3 $span

          /usr/bin/time -v -o time.txt $RSHASH build -i "$f" -d "${BASENAME}.dict" -k $k --m1 $m1 --t1 $t1 --m2 $m2 --t2 $t2 --m3 $m3 --t3 $t3 -s $span > out.txt 2>&1

          cat out.txt >> $LOG

          file_size=$(stat -c%s "${BASENAME}.dict")
          buildtime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
          buildmem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
          textlength=$(grep "text length: " out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
          density_r1=$(awk '/density r1/ {print $3}' out.txt | tr -d '%')
          density_r2=$(awk '/density r2/ {print $3}' out.txt | tr -d '%')
          density_r3=$(awk '/density r3/ {print $3}' out.txt | tr -d '%')
          density_s1=$(awk '/density s1/ {print $3}' out.txt | tr -d '%')
          density_s2=$(awk '/density s2/ {print $3}' out.txt | tr -d '%')
          density_s3=$(awk '/density s3/ {print $3}' out.txt | tr -d '%')
          density_ht=$(grep "no kmers HT:" out.txt | awk '{print $(NF)}' | sed 's/%//')
          no_kmers=$(awk '/text kmers/ {print $3}' out.txt)
          no_minimiser1=$(awk '/no minimiser1/ {print $3}' out.txt)
          no_minimiser2=$(awk '/no minimiser2/ {print $3}' out.txt)
          no_minimiser3=$(awk '/no minimiser3/ {print $3}' out.txt)
          no_distinct_minimiser1=$(awk '/no distinct minimiser1/ {print $4}' out.txt)
          no_distinct_minimiser2=$(awk '/no distinct minimiser2/ {print $4}' out.txt)
          no_distinct_minimiser3=$(awk '/no distinct minimiser3/ {print $4}' out.txt)
          spaceoffsets1=$(grep '^offsets1:' out.txt | cut -d':' -f2 | xargs)
          spaceoffsets2=$(grep '^offsets2:' out.txt | cut -d':' -f2 | xargs)
          spaceoffsets3=$(grep '^offsets3:' out.txt | cut -d':' -f2 | xargs)
          spaceht=$(grep '^Hashtable:' out.txt | cut -d':' -f2 | xargs)
          spacer1=$(grep '^R_1:' out.txt | cut -d':' -f2 | xargs)
          spacer2=$(grep '^R_2:' out.txt | cut -d':' -f2 | xargs)
          spacer3=$(grep '^R_3:' out.txt | cut -d':' -f2 | xargs)
          spaces1=$(grep '^S_1:' out.txt | cut -d':' -f2 | xargs)
          spaces2=$(grep '^S_2:' out.txt | cut -d':' -f2 | xargs)
          spaces3=$(grep '^S_3:' out.txt | cut -d':' -f2 | xargs)
          spacetotal=$(grep '^total:' out.txt | cut -d':' -f2 | xargs)
          spacetotalreal=$(echo "scale=2; $file_size / $no_kmers * 8" | bc)

          parent_dir=$(dirname "$f")
          file_name=$(basename "$f")

          /usr/bin/time -v -o time.txt $RSHASH lookup -d "${BASENAME}.dict" > out.txt 2>&1

          cat out.txt >> $LOG
          
          querymem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
          postimekmer=$(grep 'pos_time_per_kmer' out.txt | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
          negtimekmer=$(grep 'neg_time_per_kmer' out.txt | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
          posfound=$(grep 'num_positive_kmers' out.txt | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
          negfound=$(grep 'num_negative_kmers' out.txt | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
          
          echo "$f,$k,$m1,$m2,$m3,$t1,$t2,$t3,$span,$buildtime,$buildmem,$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spaceht,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$density_ht,$spacetotal,$spacetotalreal,$negtimekmer,$postimekmer,$querymem" >> "$CSV"
            
        done

      fi
    fi

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,k,m1,m2,m3,t1,t2,t3,span,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceoffsets3 [bits/kmer],spaceHT [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceR3 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],spaceS3 [bits/kmer],kmers HT [%],space theo [bits/kmer],space real [bits/kmer],neglookup [ns/kmer],poslookup [ns/kmer],querymem [B]" > "$CSV"
  run $data/
done
