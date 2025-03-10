#!/bin/bash

PROGRAM="../build/source/main"

today=$(date +%Y-%m-%d-%H-%M-%S)

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../datasets" >/dev/null 2>&1 && pwd )"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../test/datasets" >/dev/null 2>&1 && pwd )"
LOG="log.txt"
CSV="results-$today.csv"

m=15

run()
{
  FILES=$(ls $1*_text.fasta)

  for f in $FILES
  do
    echo $f >> $LOG

    BASENAME=$(echo "$f" | sed 's/_text\.fasta$//')
    echo $BASENAME
    QUERY="${BASENAME}_query.fasta"
    k=$(echo "${BASENAME##*k}" | grep -o '[0-9]*')

    /usr/bin/time -l -o time.txt $PROGRAM build -i "$f" -d "${BASENAME}.dict" -k $k -m $m > prog_out.txt 2>&1

    cat prog_out.txt >> $LOG

    file_size=$(stat -f%z "${BASENAME}.dict")
    
    buildtime=$(cat time.txt | grep "real" | awk '{print $1}')
    buildmem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')

    /usr/bin/time -l -o time.txt $PROGRAM query -d "${BASENAME}.dict" -q $QUERY -k $k -m $m > prog_out.txt 2>&1

    cat prog_out.txt >> $LOG
    
    querytime=$(cat time.txt | grep "real" | awk '{print $1}')
    querymem=$(cat time.txt  | grep "maximum resident set size" | awk '{print $1}')
    
    
    echo "$f,$buildtime,$buildmem",$file_size",$querytime,$querymem" >> "$CSV"

    # rm -f time.txt
    # rm -f prog_out.txt
  done
}


for data in $(find $DIR -mindepth 1 -maxdepth 1 -type d); do
  FILENAME=$(basename $data)
  echo "file,buildtime,buildmem,space,querytime,querymem" > "$CSV"
  run $data/
done
