#!/bin/bash

# hash set (lookup) and hash map (locate) baselines for the rshash benchmarks
# (benchmark_rshash_lookup.sh, benchmark_rshash_locate.sh): same datasets, queries, shape configurations and
# csv columns; columns without a hash table equivalent (minimisers, levels, index file) stay empty.
# the table is built once per configuration and queried with both query files in the same process.

PROGRAM="../build/test/hashtable"

today=$(date +%Y-%m-%d-%H-%M-%S)

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../datasets" >/dev/null 2>&1 && pwd )"
# own file names, so that this can run in parallel to the rshash benchmarks
LOG="ht_log.txt"
OUT="ht_out.txt"
TIME="ht_time.txt"
CSV="ht-results-$today.csv"

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
  "19:549754826455"                                         # 1,19 -- 45.8487%
  "19:549754951479 523369971711"                            # 2,19 -- 63.5086%
  "19:549754973103 526670363831 530072469503"               # 3,19 -- 69.892%
  "19:536149491687 509952393215 549754925919 463856438135"  # 4,19 -- 73.6257%
  "21:549755579099"                                         # 1,21 -- 43.3296%
  "21:549755579099 471022436351"                            # 2,21 -- 59.973%
  "21:471022436351 549755573103 491773754543"               # 3,21 -- 64.0283%
  "21:511480692735 115857162239 68719449967 549755606231"   # 4,21 -- 67.0643%
  "23:549755755163"                                         # 1,23 -- 39.535%
  "23:549755754843 460442304511"                            # 2,23 -- 53.0609%
  "23:549755754843 491773755175 469434892287"               # 3,23 -- 55.6012%
  "23:34359735031 549755753679 33663483903 520404074495"    # 4,23 -- 57.9174%
)

# value of "<name> = <value>" in a report
get()
{
    grep "^$1 = " | head -1 | sed -E "s/^$1 = ([^ ]+).*/\1/"
}

run()
{
    f="$1"
    BASENAME=$(echo "$f" | sed 's/\.fa.gz$//')
    parent_dir=$(dirname "$f")
    k=31

    # low and high hitrate query of the dataset (as in the rshash benchmarks)
    if [[ "$BASENAME" == *"bacterial"* ]]; then
      queries=("${parent_dir/unitigs/queries}/SRR5833294.fastq.gz" "${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz")
    elif [[ "$BASENAME" == *"cod"* ]]; then
      queries=("${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz" "${parent_dir/unitigs/queries}/SRR12858649.fastq.gz")
    elif [[ "$BASENAME" == *"kestrel"* ]]; then
      queries=("${parent_dir/unitigs/queries}/SRR12858649.fastq.gz" "${parent_dir/unitigs/queries}/SRR11449743_1.fastq.gz")
    elif [[ "$BASENAME" == *"human"* ]]; then
      queries=("${parent_dir/unitigs/queries}/SRR5901135_1.fastq.gz" "${parent_dir/unitigs/queries}/SRR5833294.fastq.gz")
    else
      queries=()
    fi
    query_args=()
    for query in "${queries[@]}"; do
      query_args+=( -q "$query" )
    done

    for structure in hashset hashmap; do
      mode_args=()
      if [[ "$structure" == "hashmap" ]]; then
        mode_args=( --locate )
      fi

      for config in "${shape_configs[@]}"; do
        IFS=':' read -r core_length shape_values <<< "$config"
        read -ra shapes <<< "$shape_values"
        no_shapes=${#shapes[@]}
        shape_args=()
        for s in "${shapes[@]}"; do
          shape_args+=( --shapes "$s" )
        done
        if [[ $no_shapes -eq 0 ]]; then
          shape_args=( -k "$core_length" )
        fi

        echo $BASENAME $structure
        echo "shapes: $no_shapes, core length: $core_length, ${shapes[*]}"
        echo $f $structure >> $LOG
        echo "shapes: $no_shapes, core length: $core_length, ${shapes[*]}" >> $LOG

        /usr/bin/time -v -o $TIME $PROGRAM -i "$f" "${shape_args[@]}" "${mode_args[@]}" --bench "${query_args[@]}" > $OUT 2>&1
        status=$?

        cat $OUT >> $LOG
        cat $TIME >> $LOG

        if [[ $status -ne 0 ]]; then
          echo "$structure failed (shapes: $no_shapes, core length: $core_length), skipping" | tee -a $LOG
          continue
        fi

        # build: wall clock time [s] and peak memory [kB] after the build, measured in the program
        buildtime=$(get build_time < $OUT)
        buildmem=$(get build_mem < $OUT)
        spaceht=$(grep '^Last_level:' $OUT | cut -d':' -f2 | xargs)
        spacetotal=$(grep '^total:' $OUT | cut -d':' -f2 | xargs)
        memtotal=$spacetotal

        # single k-mers/shapes: random text k-mers (pos) and random k-mers (neg)
        single_lookup_pos=$(get pos_time_per_kmer < $OUT)
        single_lookup_neg=$(get neg_time_per_kmer < $OUT)
        single_locate_pos=$(get pos_locate_time_per_kmer < $OUT)
        single_locate_neg=$(get neg_locate_time_per_kmer < $OUT)
        single_locate_kmers=$(grep -A3 '^==== positive locate:' $OUT | get num_kmers)
        single_locate_positions=$(grep -A3 '^==== positive locate:' $OUT | get num_positions)
        single="$single_lookup_pos,$single_lookup_neg,$single_locate_pos,$single_locate_neg,$single_locate_kmers,$single_locate_positions"

        for query in "${queries[@]}"; do
          report=$(awk -v q="==== query report: $query" '$0 == q {p = 1} p {print} /^==== end query report$/ {p = 0}' $OUT)

          querytimekmer=$(get time_per_kmer <<< "$report")
          querymem=$(get query_mem <<< "$report")
          k_mers=$(get num_kmers <<< "$report")
          found=$(get num_positive_kmers <<< "$report")
          positions=$(get num_positions <<< "$report")
          sensitivity=$(grep "^num_found_queries" <<< "$report" | sed -E 's/.*num_found_queries = [0-9]+ \(([0-9.e+-]+)%\).*/\1/')

          # columns as in benchmark_rshash_locate.sh (rshash specific ones empty) + structure
          echo "$f,$query,$k,$no_shapes,$core_length,${shapes[*]},,,,,,,$buildtime,$buildmem,,,,,,,,,,,,,$spaceht,,,,,,,,,,,,,,,$memtotal,$querytimekmer,$querymem,,$k_mers,$found,$positions,$sensitivity,$single,$structure" >> "$CSV"
        done
      done
    done
}

FILES=( "../../datasets/cod.genome.preproc.fa" "../../datasets/kestrel.genome.preproc.fa.gz" "../../datasets/human.genome.preproc.fa.gz" )
echo "textfile,queryfile,k,no shapes,core length,shapes,m1,m2,m3,t1,t2,t3,buildtime [s],buildmem [B],indexsize [B],spaceoffsets1 [bits/kmer],spaceoffsets2 [bits/kmer],spaceoffsets3 [bits/kmer],spaceoffsets4 [bits/kmer],spaceR1 [bits/kmer],spaceR2 [bits/kmer],spaceR3 [bits/kmer],spaceS1 [bits/kmer],spaceS2 [bits/kmer],spaceS3 [bits/kmer],spaceS4 [bits/kmer],last level [bits/kmer],density_r1 [%],density_r2 [%],density_r3 [%],density_s1 [%],density_s2 [%],density_s3 [%],freq kmers [%],no minimizer1,no minimizer2,no minimizer3, no distinct minimizer1,no distinct minimizer2,no distinct minimizer3,space disk [bits/kmer],mem theo [bits/kmer],querytime [ns/kmer],querymem [B],extensions,kmers,found,positions,sensitivity [%],single lookup pos [ns/kmer],single lookup neg [ns/kmer],single locate pos [ns/kmer],single locate neg [ns/kmer],single locate kmers,single locate positions,structure" > "$CSV"
for f in "${FILES[@]}"; do
  run $f
done
