#!/bin/bash

PROGRAM="../build/source/main"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../datasets" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="results-$today.csv"

m=10

run_program()
{
  FILES=$(ls $1*_text.fasta)

  for f in $FILES
  do
    echo $f >> $LOG

    BASENAME=$(echo "$f" | sed 's/_text\.fasta$//')
    echo $BASENAME
    QUERY="${BASENAME}_query.fasta"
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')

    /usr/bin/time -l -o time.txt $PROGRAM bq -i "$f" -q "$QUERY" -k $k -m $m > prog_out.txt 2>&1

    cat prog_out.txt >> $LOG

    correct="";
    if cmp -s "${BASENAME}.positions" prog_out.txt; then
      echo "correct."
      correct="1";
    else
      echo "solution not correct."
      correct="0";
    fi
    
    ELAPSED_TIME=$(cat time.txt | grep "real" | awk '{print $1}')
    # USER_TIME=$(cat time.txt  | grep "user" | awk '{print $1}')
    # SYSTEM_TIME=$(cat time.txt  | grep "sys" | awk '{print $1}')
    MAX_RSS=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
    
    # echo "$f,$ELAPSED_TIME,$USER_TIME,$SYSTEM_TIME,$MAX_RSS" >> "$CSV"
    echo "$f,$ELAPSED_TIME,$MAX_RSS",$correct >> "$CSV"

    # rm -f time.txt
    # rm -f prog_out.txt
  done
}


for data in $(find $DIR -mindepth 1 -maxdepth 1 -type d); do
  FILENAME=$(basename $data)
  echo "file,time,memory,correct" > "$CSV"
  run_program $data/
done
