#!/bin/bash

PROGRAM="../build/source/lookup"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="lookup-results-$today.csv"

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

        ms=()
        for ((i=-1; i<=0; i++)); do
            m=$(echo "l($length)/l(4)+$i" | bc -l)
            m=$(printf "%.0f" "$m")
            ms+=("$m")
        done

        thresholds=(100 50)

        for m in "${ms[@]}"; do
          for thres in "${thresholds[@]}"; do

            echo $BASENAME
            echo $m

            echo $f >> $LOG

            /usr/bin/time -v -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -m $m -t $thres > prog_out.txt 2>&1

            cat prog_out.txt >> $LOG

            file_size=$(stat -f%z "${BASENAME}.dict")
            buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
            buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
            textlength=$(grep "text length: " prog_out.txt | sed -E 's/.*text length: ([0-9]+).*/\1/')
            density_r=$(awk '/density r/ {print $3}' prog_out.txt | tr -d '%')
            density_s=$(awk '/density s/ {print $3}' prog_out.txt | tr -d '%')
            density_ht=$(grep "no kmers HT:" prog_out.txt | awk '{print $(NF)}' | sed 's/%//')
            no_kmers=$(awk '/text kmers/ {print $3}' prog_out.txt)
            no_minimiser=$(awk '/no minimiser/ {print $3}' prog_out.txt)
            no_distinct_minimiser=$(awk '/no distinct minimiser/ {print $4}' prog_out.txt)
            spaceoffsets=$(grep '^offsets:' prog_out.txt | cut -d':' -f2 | xargs)
            spacer=$(grep '^R:' prog_out.txt | cut -d':' -f2 | xargs)
            spaces=$(grep '^S:' prog_out.txt | cut -d':' -f2 | xargs)
            spaceht=$(grep '^Hashtable:' prog_out.txt | cut -d':' -f2 | xargs)
            spacetotal=$(grep '^total:' prog_out.txt | cut -d':' -f2 | xargs)
            spacetotalreal=$(echo "scale=2; $file_size/ $no_kmers * 8" | bc -l)

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
                /usr/bin/time -v -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query > prog_out.txt 2>&1

                cat prog_out.txt >> $LOG
                
                querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                
                echo "$f,$query,$k,$m,$thres,$buildtime,$buildmem",$file_size,$spaceoffsets,$spaceht,$spacer,$spaces,$density_r,$density_s,$density_ht,$no_minimiser,$no_distinct_minimiser,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

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
                /usr/bin/time -v -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $query > prog_out.txt 2>&1

                cat prog_out.txt >> $LOG
                
                querytime=$(cat time.txt | grep "real" | awk '{print $1}')
                querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
                k_mers=$(grep "num_kmers" prog_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                found=$(grep "num_positive_kmers" prog_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                
                echo "$f,$query,$k,$m,$thres,$buildtime,$buildmem",$file_size,$spaceoffsets,$spaceht,$spacer,$spaces,$density_r,$density_s,$density_ht,$no_minimiser,$no_distinct_minimiser,$spacetotal,$spacetotalreal,$querytimekmer,$querymem,$k_mers",$found" >> "$CSV"

                # rm -f time.txt
                # rm -f prog_out.txt
            done

          done
        done

      fi
    fi

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,thres,buildtime [s],buildmem [B],indexsize [B],spaceoffsets [bits/kmer],space HT [bits/kmer],spaceR [bits/kmer],spaceS [bits/kmer],density_r [%],density_s [%],kmers HT [%],no minimizer, no distinct minimizer,space theo [bits/kmer],space real [bits/kmer],querytime [ns/kmer],querymem [B],kmers,found" > "$CSV"
  run $data/
done
