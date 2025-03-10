#!/bin/bash

PROGRAM="../build/source/main"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/datasets" >/dev/null 2>&1 && pwd )"
LOG="log.txt"

m=15

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

    $PROGRAM bq -i "$f" -q "$QUERY" -k $k -m $m -p > prog_out.txt 2>&1

    cat prog_out.txt >> $LOG

    correct="";
    if cmp -s "${BASENAME}.positions" prog_out.txt; then
      echo "correct."
      correct="1";
    else
      echo "solution not correct."
      correct="0";
    fi
    
  done
}


for data in $(find $DIR -mindepth 1 -maxdepth 1 -type d); do
  FILENAME=$(basename $data)
  run_program $data/
done
