#!/bin/bash

PROGRAM="../../sshash/build/sshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets/unitigs" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="pibiri-lookupresults-$today.csv"

m=15

run()
{
  FILES=$(ls $1*.fa.gz)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')

    if [[ "$BASENAME" == *"bacterial"* || "$BASENAME" == *"human"* || "$BASENAME" == *"cod"* || "$BASENAME" == *"kestrel"* ]]; then
      if [ "$k" -le 32 ]; then
        echo $BASENAME

        echo $f >> $LOG

        /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -o "${BASENAME}.index" -k $k -m $m > prog_out.txt 2>&1

        cat prog_out.txt >> $LOG

        # file_size=$(stat -f%z "${BASENAME}.dict")
        file_size=$(stat -f%z "${BASENAME}.index")
        buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
        buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')

        parent_dir=$(dirname "$f")
        file_name=$(basename "$f")

        # if [[ "$BASENAME" == *"celegans"* ]]; then
        #   query="${parent_dir/unitigs/queries}/SRR16288382_1.fastq.gz"
        # elif [[ "$BASENAME" == *"cod"* ]]; then
        #   query="${parent_dir/unitigs/queries}/SRR12858649.fastq.gz"
        # elif [[ "$BASENAME" == *"kestrel"* ]]; then
        #   query="${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz"
        # elif [[ "$BASENAME" == *"human"* ]]; then
        #   query="${parent_dir/unitigs/queries}/SRR5833294.fastq.gz"
        # fi
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

          /usr/bin/time -l -o time.txt $PROGRAM query -i "${BASENAME}.index" -q $query > prog_out.txt 2>&1

          cat prog_out.txt >> $LOG
          
          querytime=$(cat time.txt | grep "real" | awk '{print $1}')
          querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')

          k_mers=$(grep "num_kmers" prog_out.txt | awk -F'=' '{print $2}' | xargs)
          found=$(grep "num_positive_kmers" prog_out.txt | awk -F'=' '{print $2}' | awk '{print $1}')
          
          echo "$f,$query,$buildtime,$buildmem",$file_size,$querytime,$querymem,$k_mers",$found" >> "$CSV"

          # rm -f time.txt
          # rm -f prog_out.txt
        done
      fi
    fi
  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,buildtime,buildmem,space,querytime,querymem,kmers,found" > "$CSV"
  run $data/
done
