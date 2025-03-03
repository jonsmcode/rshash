#!/bin/bash

PROGRAM="../build/source/main"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../datasets" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="results-$today.csv"

k=20
m=10


run_program()
{
  FILES=$(ls $1*_text.fasta)

  for f in $FILES
  do
    BASENAME=$(echo "$f" | sed 's/_text\.fasta$//')
    echo $BASENAME
    QUERY="${BASENAME}_query.fasta"

    /usr/bin/time -l -o time.txt $PROGRAM bq -i "$f" -q "$QUERY" -k $k -m $m > prog_out.txt 2>&1

    if cmp -s "${BASENAME}.positions" prog_out.txt; then
      echo "correct."
    else
      echo "solution not correct."
    fi
    
    ELAPSED_TIME=$(cat time.txt | grep "real" | awk '{print $1}')
    USER_TIME=$(cat time.txt  | grep "user" | awk '{print $1}')
    SYSTEM_TIME=$(cat time.txt  | grep "sys" | awk '{print $1}')
    MAX_RSS=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
    
    echo "$f,$ELAPSED_TIME,$USER_TIME,$SYSTEM_TIME,$MAX_RSS" >> "$CSV"

    rm -f time.txt
    rm -f prog_out.txt
  done
}


for data in $(find $DIR -mindepth 1 -maxdepth 1 -type d); do
  FILENAME=$(basename $data)
  echo "File,Elapsed Time (s),User Time (s),System Time (s),Maximum Resident Set Size (KB)" > "$CSV"
  run_program $data/
done
