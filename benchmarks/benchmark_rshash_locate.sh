#!/bin/bash

PROGRAM="../build/source/rshash"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets" >/dev/null 2>&1 && pwd )"
LOG="locate_log.txt"
OUT="locateout.txt"
CSV="rshash-locate-results-$today.csv"

# shapes (binary, decimal) per configuration "<no shapes>,<core length> -- <hit rate>"
# 1,19 -- 45.8487%
# 111111111111111111100001110111011010111 549754826455

# 2,19 -- 63.5086%
# 111111111111111111100101101011100110111 549754951479
# 111100111011011010001111111111111111111 523369971711

# 3,19 -- 69.892%
# 111111111111111111100110010101110101111 549754973103
# 111101010011111111111111111110010110111 526670363831
# 111101101101010110001111111111111111111 530072469503

# 4,19 -- 73.6257%
# 111110011010100111111111111111111100111 536149491687
# 111011010111011100001111111111111111111 509952393215
# 111111111111111111100100111001101011111 549754925919
# 110101111111111111111111000101101110111 463856438135

# 1,21 -- 43.3296%
# 111111111111111111111000110101011011011 549755579099

# 2,21 -- 59.973%
# 111111111111111111111000110101011011011 549755579099
# 110110110101011000111111111111111111111 471022436351

# 3,21 -- 64.0283%
# 110110110101011000111111111111111111111 471022436351
# 111111111111111111111000101001101101111 549755573103
# 111001001111111111111111111110010101111 491773754543

# 4,21 -- 67.0643%
# 111011100010110100111111111111111111111 511480692735
# 01101011111001100111111111111111111111  115857162239
# 111111111111111111111001011101101111    68719449967
# 111111111111111111111001101010011010111 549755606231

# 1,23 -- 39.535%
# 111111111111111111111110001101010011011 549755755163

# 2,23 -- 53.0609%
# 111111111111111111111110001100101011011 549755754843
# 110101100110100011111111111111111111111 460442304511

# 3,23 -- 55.6012%
# 111111111111111111111110001100101011011 549755754843
# 111001001111111111111111111111100100111 491773755175
# 110110101001100011111111111111111111111 469434892287

# 4,23 -- 57.9174%
# 11111111111111111111111001011110111     34359735031
# 111111111111111111111110001010011001111 549755753679
# 011111010110011111111111111111111111    33663483903
# 111100100101010011111111111111111111111 520404074495

# shape configurations "core length:shape values" (from shape optimisation; % = hit rate);
# no shape values: contiguous k-mers with k = core length
shape_configs=(
  "31:"                                                     # 31-mers
  "19:"                                                     # 19-mers
  "19:549754826455"                                         # 1,19 -- 45.8487%
  "19:549754951479 523369971711"                            # 2,19 -- 63.5086%
  "19:549754973103 526670363831 530072469503"               # 3,19 -- 69.892%
  "19:536149491687 509952393215 549754925919 463856438135"  # 4,19 -- 73.6257%
  "21:"                                                     # 21-mers
  "21:549755579099"                                         # 1,21 -- 43.3296%
  "21:549755579099 471022436351"                            # 2,21 -- 59.973%
  "21:471022436351 549755573103 491773754543"               # 3,21 -- 64.0283%
  "21:511480692735 115857162239 68719449967 549755606231"   # 4,21 -- 67.0643%
  "23:"                                                     # 23-mers
  "23:549755755163"                                         # 1,23 -- 39.535%
  "23:549755754843 460442304511"                            # 2,23 -- 53.0609%
  "23:549755754843 491773755175 469434892287"               # 3,23 -- 55.6012%
  "23:34359735031 549755753679 33663483903 520404074495"    # 4,23 -- 57.9174%
)

run()
{
    f="$1"
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    k=31
        # params per shape core length: longer cores use longer minimisers
        if [[ "$BASENAME" == *"human"* ]]; then
          params19=( 15 17 0 64 128 0 2  17 0 0 128 0 0 1 )
          params21=( 17 19 0 64 128 0 2  19 0 0 128 0 0 1 )
          params23=( 17 21 0 64 128 0 2  21 0 0 128 0 0 1 )
          params31=( 17 23 25 64 64 64 3  19 0 0 128 0 0 1 )
        elif [[ "$BASENAME" == *"cod"* ]]; then
          params19=( 15 17 0 64 128 0 2  17 0 0 128 0 0 1 )
          params21=( 17 19 0 64 128 0 2  19 0 0 128 0 0 1 )
          params23=( 17 21 0 64 128 0 2  21 0 0 128 0 0 1 )
          params31=( 17 23 25 64 64 64 3  19 0 0 128 0 0 1 )
        elif [[ "$BASENAME" == *"kestrel"* ]]; then
          params19=( 15 17 0 64 128 0 2  17 0 0 128 0 0 1 )
          params21=( 17 19 0 64 128 0 2  19 0 0 128 0 0 1 )
          params23=( 17 21 0 64 128 0 2  21 0 0 128 0 0 1 )
          params31=( 17 23 25 64 64 64 3  19 0 0 128 0 0 1 )
        fi

        for config in "${shape_configs[@]}"; do
          IFS=':' read -r core_length shape_values <<< "$config"
          read -ra shapes <<< "$shape_values"
          no_shapes=${#shapes[@]}
          shape_args=()
          for s in "${shapes[@]}"; do
            shape_args+=( --shape "$s" )
          done
          if [[ $no_shapes -eq 0 ]]; then
            shape_args=( -k "$core_length" )
          fi
          p="params${core_length}[@]"
          params=( "${!p}" )

          for ((i=0; i<${#params[@]}; i+=7)); do
          m1=${params[i]}
          m2=${params[i+1]}
          m3=${params[i+2]}
          t1=${params[i+3]}
          t2=${params[i+4]}
          t3=${params[i+5]}
          l=${params[i+6]}

          echo $BASENAME
          echo $m1 $t1 $m2 $t2 $m3 $t3
          echo "shapes: $no_shapes, core length: $core_length, ${shapes[*]}"
          echo $f >> $LOG
          echo $m1 $t1 $m2 $t2 $m3 $t3 >> $LOG
          echo "shapes: $no_shapes, core length: $core_length, ${shapes[*]}" >> $LOG

          rm -f "${BASENAME}.loc.rshash"
          /usr/bin/time -v -o time_loc.txt $PROGRAM build -i "$f" -d "${BASENAME}.loc.rshash" "${shape_args[@]}" --level $l --m1 $m1 --m2 $m2 --m3 $m3 --t1 $t1 --t2 $t2 --t3 $t3 --loc > OUT 2>&1
          build_status=$?

          cat OUT >> $LOG

          if [[ $build_status -ne 0 || ! -f "${BASENAME}.loc.rshash" ]]; then
            echo "build failed (shapes: $no_shapes, core length: $core_length), skipping" | tee -a $LOG
            continue
          fi

            file_size=$(stat -c%s "${BASENAME}.loc.rshash")
            # wall clock time in seconds (the last level is built in parallel), given as h:mm:ss or m:ss.ss
            buildtime=$(grep "Elapsed (wall clock)" time_loc.txt | awk '{n = split($NF, t, ":"); s = 0; for (i = 1; i <= n; i++) s = s * 60 + t[i]; print s}')
            buildmem=$(grep "Maximum resident set size" time_loc.txt | awk -F': ' '{print $2}')
            textlength=$(grep "text length: " OUT | sed -E 's/.*text length: ([0-9]+).*/\1/')
            density_r1=$(awk '/density r1/ {print $3}' OUT | tr -d '%')
            density_r2=$(awk '/density r2/ {print $3}' OUT | tr -d '%')
            density_r3=$(awk '/density r3/ {print $3}' OUT | tr -d '%')
            density_s1=$(awk '/density s1/ {print $3}' OUT | tr -d '%')
            density_s2=$(awk '/density s2/ {print $3}' OUT | tr -d '%')
            density_s3=$(awk '/density s3/ {print $3}' OUT | tr -d '%')
            density_ht=$(grep "no freq kmers:" OUT | awk '{print $(NF)}' | sed 's/%//')
            no_kmers=$(awk '/textkmers/ {print $2}' OUT)
            no_minimiser1=$(awk '/no minimiser1/ {print $3}' OUT)
            no_minimiser2=$(awk '/no minimiser2/ {print $3}' OUT)
            no_minimiser3=$(awk '/no minimiser3/ {print $3}' OUT)
            no_distinct_minimiser1=$(awk '/no distinct minimiser1/ {print $4}' OUT)
            no_distinct_minimiser2=$(awk '/no distinct minimiser2/ {print $4}' OUT)
            no_distinct_minimiser3=$(awk '/no distinct minimiser3/ {print $4}' OUT)
            spaceoffsets1=$(grep '^offsets1:' OUT | cut -d':' -f2 | xargs)
            spaceoffsets2=$(grep '^offsets2:' OUT | cut -d':' -f2 | xargs)
            spaceoffsets3=$(grep '^offsets3:' OUT | cut -d':' -f2 | xargs)
            spaceoffsets4=$(grep '^offsets4:' OUT | cut -d':' -f2 | xargs)
            spaceht=$(grep '^Last_level:' OUT | cut -d':' -f2 | xargs)
            spacer1=$(grep '^R_1:' OUT | cut -d':' -f2 | xargs)
            spacer2=$(grep '^R_2:' OUT | cut -d':' -f2 | xargs)
            spacer3=$(grep '^R_3:' OUT | cut -d':' -f2 | xargs)
            spaces1=$(grep '^S_1:' OUT | cut -d':' -f2 | xargs)
            spaces2=$(grep '^S_2:' OUT | cut -d':' -f2 | xargs)
            spaces3=$(grep '^S_3:' OUT | cut -d':' -f2 | xargs)
            spaces4=$(grep '^S_4:' OUT | cut -d':' -f2 | xargs)
            spacetotal=$(grep '^total:' OUT | cut -d':' -f2 | xargs)
            memtotal=$(grep '^total:' OUT | cut -d':' -f2 | xargs)
            spacetotal=$(echo "scale=2; $file_size / $no_kmers * 8" | bc)

            # single k-mers/shapes: lookup and locate of random text k-mers (pos) and random k-mers (neg)
            $PROGRAM bench -d "${BASENAME}.loc.rshash" > BENCH 2>&1
            cat BENCH >> $LOG
            single_lookup_pos=$(grep '^pos_time_per_kmer =' BENCH | cut -d= -f2 | xargs)
            single_lookup_neg=$(grep '^neg_time_per_kmer =' BENCH | cut -d= -f2 | xargs)
            single_locate_pos=$(grep '^pos_locate_time_per_kmer =' BENCH | cut -d= -f2 | xargs)
            single_locate_neg=$(grep '^neg_locate_time_per_kmer =' BENCH | cut -d= -f2 | xargs)
            single_locate_kmers=$(grep -A3 '^==== positive locate:' BENCH | grep '^num_kmers =' | cut -d= -f2 | xargs)
            single_locate_positions=$(grep -A3 '^==== positive locate:' BENCH | grep '^num_positions =' | cut -d= -f2 | xargs)
            single="$single_lookup_pos,$single_lookup_neg,$single_locate_pos,$single_locate_neg,$single_locate_kmers,$single_locate_positions"

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
                /usr/bin/time -v -o time_loc.txt $PROGRAM locate -d "${BASENAME}.loc.rshash" -q $query > OUT 2>&1

                cat OUT >> $LOG
              
                querytime=$(grep "User time" time_loc.txt | awk -F': ' '{print $2}')
                querymem=$(grep "Maximum resident set size" time_loc.txt | awk -F': ' '{print $2}')
                k_mers=$(grep "num_kmers" OUT | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                found=$(grep "num_positive_kmers" OUT | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                positions=$(grep "num_positions" OUT | sed -E 's/.*num_positions = ([0-9]+).*/\1/')
                sensitivity=$(grep "num_found_queries" OUT | sed -E 's/.*num_found_queries = [0-9]+ \(([0-9.e+-]+)%\).*/\1/')
                # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                querytimekmer=$(grep 'time_per_kmer' OUT | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
                extensions=$(grep 'extensions' OUT | awk -F'=' '{print $2}' | awk '{print $1}')

                echo "$f,$query,$k,$no_shapes,$core_length,${shapes[*]},$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spaceoffsets4,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$spaces4,$spaceht,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$memtotal,$querytimekmer,$querymem,$extensions,$k_mers",$found,$positions,$sensitivity,$single" >> "$CSV"

                # rm -f time_loc.txt
                # rm -f OUT
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
                /usr/bin/time -v -o time_loc.txt $PROGRAM locate -d "${BASENAME}.loc.rshash" -q $query > OUT 2>&1

                cat OUT >> $LOG
              
                querytime=$(grep "User time" time_loc.txt | awk -F': ' '{print $2}')
                querymem=$(grep "Maximum resident set size" time_loc.txt | awk -F': ' '{print $2}')
                k_mers=$(grep "num_kmers" OUT | sed -E 's/.*num_kmers = ([0-9]+).*/\1/')
                found=$(grep "num_positive_kmers" OUT | sed -E 's/.*num_positive_kmers = ([0-9]+).*/\1/')
                positions=$(grep "num_positions" OUT | sed -E 's/.*num_positions = ([0-9]+).*/\1/')
                sensitivity=$(grep "num_found_queries" OUT | sed -E 's/.*num_found_queries = [0-9]+ \(([0-9.e+-]+)%\).*/\1/')
                # querytimekmer=$(echo "scale=10; $querytime / $k_mers * 1000000000" | bc)
                querytimekmer=$(grep 'time_per_kmer' OUT | awk -F'=' '{print $2}' | awk '{print $1}' | sed 's/ns//')
                extensions=$(grep 'extensions' OUT | awk -F'=' '{print $2}' | awk '{print $1}')
              
                echo "$f,$query,$k,$no_shapes,$core_length,${shapes[*]},$m1,$m2,$m3,$t1,$t2,$t3,$buildtime,$buildmem",$file_size,$spaceoffsets1,$spaceoffsets2,$spaceoffsets3,$spaceoffsets4,$spacer1,$spacer2,$spacer3,$spaces1,$spaces2,$spaces3,$spaces4,$spaceht,$density_r1,$density_r2,$density_r3,$density_s1,$density_s2,$density_s3,$density_ht,$no_minimiser1,$no_minimiser2,$no_minimiser3,$no_distinct_minimiser1,$no_distinct_minimiser2,$no_distinct_minimiser3,$spacetotal,$memtotal,$querytimekmer,$querymem,$extensions,$k_mers",$found,$positions,$sensitivity,$single" >> "$CSV"

                # rm -f time_loc.txt
                # rm -f OUT
            done

        
          done
        done

}

# FILES=( "../DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs.fa.ust.fa.gz" "../DNA_datasets/Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs.fa.ust.fa.gz" "../DNA_datasets/Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs.fa.ust.fa.gz" "../DNA_datasets/bacterial.genome.fixed.fa.unitigs.fa.ust.fa.gz" )
FILES=( "../../datasets/cod.genome.preproc.fa" "../../datasets/kestrel.genome.preproc.fa.gz" "../../datasets/human.genome.preproc.fa.gz" )
echo "textfile,queryfile,k,no shapes,core length,shapes,m1,m2,m3,t1,t2,t3,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceoffsets3 [bits/kmer],spaceoffsets4 [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceR3 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],spaceS3 [bits/kmer],spaceS4 [bits/kmer],last level [bits/kmer],density_r1 [%],density_r2 [%],density_r3 [%],density_s1 [%],density_s2 [%],density_s3 [%],freq kmers [%],no minimizer1,no minimizer2,no minimizer3, no distinct minimizer1,no distinct minimizer2,no distinct minimizer3,space disk [bits/kmer],mem theo [bits/kmer],querytime [ns/kmer],querymem [B],extensions,kmers,found,positions,sensitivity [%],single lookup pos [ns/kmer],single lookup neg [ns/kmer],single locate pos [ns/kmer],single locate neg [ns/kmer],single locate kmers,single locate positions" > "$CSV"
for f in "${FILES[@]}"; do
  run $f
done
