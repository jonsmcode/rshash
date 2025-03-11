#!/bin/bash

PROGRAM="../build/source/main"

today=$(date +%Y-%m-%d-%H-%M-%S)

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../test/datasets" >/dev/null 2>&1 && pwd )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../datasets/unitigs" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="results-$today.csv"

m=15

run()
{
  FILES=$(ls $1*.fa.gz)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')

    if [ "$k" -le 32 ]; then
      echo $BASENAME

      echo $f >> $LOG

      /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -m $m > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG

      file_size=$(stat -f%z "${BASENAME}.dict")
      
      buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
      buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')

      parent_dir=$(dirname "$f")
      file_name=$(basename "$f")

      if [[ "$BASENAME" == "$celegans"* ]]; then
        QUERY="${parent_dir/unitigs/queries}/SRR16288382.fastq"
      elif [[ "$BASENAME" == "$cod"* ]]; then
        QUERY="${parent_dir/unitigs/queries}/SRR12858649.fastq"
      elif [[ "$BASENAME" == "$kestrel"* ]]; then
        QUERY="${parent_dir/unitigs/queries}/SRR11449743.fastq"
      elif [[ "$BASENAME" == "$human"* ]]; then
        QUERY="${parent_dir/unitigs/queries}/SRR5833294.fastq"
      fi

      /usr/bin/time -l -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $QUERY -k $k -m $m > prog_out.txt 2>&1

      cat prog_out.txt >> $LOG
      
      querytime=$(cat time.txt | grep "real" | awk '{print $1}')
      querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
      
      
      echo "$f,$buildtime,$buildmem",$file_size",$querytime,$querymem" >> "$CSV"

      # rm -f time.txt
      # rm -f prog_out.txt
    fi
  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "file,buildtime,buildmem,space,querytime,querymem" > "$CSV"
  run $data/
done
