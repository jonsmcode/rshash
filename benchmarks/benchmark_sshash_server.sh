#!/bin/bash

PROGRAM="../../../sshash/build/sshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets" >/dev/null 2>&1 && pwd )"
LOG="logsshash.txt"
CSV="sshash-results-$today.csv"

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
          for ((i=1; i<=5; i++)); do
              m=$(echo "l($length)/l(4)+$i" | bc -l)
              m=$(printf "%.0f" "$m")
              ms+=("$m")
          done

          for m in "${ms[@]}"; do
                echo $BASENAME

                echo $f >> $LOG

                /usr/bin/time -v -o time.txt $PROGRAM build -i "$f" -o "${BASENAME}.index" -k $k -m $m --canonical-parsing > sshash_out.txt 2>&1

                cat sshash_out.txt >> $LOG

                file_size=$(stat -c%z "${BASENAME}.index")
                buildtime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
                buildmem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
                num_super_kmers=$(grep "^num_super_kmers" sshash_out.txt | awk '{print $2}')
                space=$(grep "total:" sshash_out.txt | sed -E 's/^ *total: *([0-9.eE+-]+).*/\1/')
                space_o=$(grep "offsets:" sshash_out.txt | sed -E 's/^ *offsets: *([0-9.eE+-]+).*/\1/')
                space_m=$(grep "minimizers:" sshash_out.txt | sed -E 's/^ *minimizers: *([0-9.eE+-]+).*/\1/')
                bits_key=$(awk '/minimizers:/ {for(i=1;i<=NF;i++) if($i ~ /\[bits\/key\]/) print $(i-1)}' sshash_out.txt | tr -d '(')

                parent_dir=$(dirname "$f")
                file_name=$(basename "$f")

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
                    /usr/bin/time -v -o time.txt $PROGRAM query -i "${BASENAME}.index" -q $query > sshash_out.txt 2>&1

                    cat sshash_out.txt >> $LOG
                    
                    querytime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
                    querymem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
                    k_mers=$(grep "num_kmers" sshash_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                    found=$(grep "num_positive_kmers" sshash_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                    # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                    querytimekmer=$(grep 'ns/kmer' sshash_out.txt | awk -F'/' '{print $4}' | awk '{print $1}')
                    num_extensions=$(grep "num_extensions" sshash_out.txt | sed -E 's/.*num_extensions = ([0-9]+).*/\1/')
                    
                    echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$space_o,$space_m,$bits_key,$num_super_kmers,$space,$querytimekmer,$querymem,$k_mers",$found",$num_extensions >> "$CSV"

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
                    /usr/bin/time -v -o time.txt $PROGRAM query -i "${BASENAME}.index" -q $query > sshash_out.txt 2>&1

                    cat sshash_out.txt >> $LOG
                    
                    querytime=$(grep "User time" time.txt | awk -F': ' '{print $2}')
                    querymem=$(grep "Maximum resident set size" time.txt | awk -F': ' '{print $2}')
                    k_mers=$(grep "num_kmers" sshash_out.txt | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                    found=$(grep "num_positive_kmers" sshash_out.txt | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                    # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                    querytimekmer=$(grep 'ns/kmer' sshash_out.txt | awk -F'/' '{print $4}' | awk '{print $1}')
                    num_extensions=$(grep "num_extensions" sshash_out.txt | sed -E 's/.*num_extensions = ([0-9]+).*/\1/')
                    
                    echo "$f,$query,$k,$m,$buildtime,$buildmem",$file_size,$space_o,$space_m,$bits_key,$num_super_kmers,$space,$querytimekmer,$querymem,$k_mers",$found,$num_extensions" >> "$CSV"

                done

          done

        fi
      fi

  done
}


for data in $(find $DIR -mindepth 0 -maxdepth 0 -type d); do
  FILENAME=$(basename $data)
  echo "textfile,queryfile,k,m,buildtime[s],buildmem[B],indexsize[B], space offsets[bits/kmer], space minimizers[bits/kmer], bits/key, no super kmers, space[bits/kmer], querytime[ns/kmer],querymem[B],kmers,found,extensions" > "$CSV"
  run $data/
done
