#!/bin/bash
# Pivot a headerless long table (sample, feature, value) into a wide
# feature-by-sample matrix, filling absent pairs with 0.
#
# Equivalent to
#   mlr --tsv --implicit-tsv-header label hash,feature,value \
#     then reshape -s hash,value then unsparsify --fill-with 0
# except that rows are sorted by feature and columns follow the samples' first
# appearance in the input. Miller holds every record in memory; this streams
# after an external sort, so memory use stays small for tables of any size.
#
# Usage: long_to_presence_matrix.sh INPUT OUTPUT TMPDIR THREADS SORT_MEMORY
set -euo pipefail

input=$1 output=$2 tmp_root=$3 threads=$4 sort_memory=$5
export LC_ALL=C

work=$(mktemp -d "$tmp_root/presence-matrix.XXXXXX")
trap 'rm -rf "$work"' EXIT

awk -F '\t' '!seen[$1]++ { print $1 }' "$input" > "$work/samples"

sort -t "$(printf '\t')" -k2,2 -s -S "$sort_memory" --parallel="$threads" -T "$work" "$input" |
awk -F '\t' -v samples="$work/samples" '
    function flush(    i, line) {
        line = feature
        for (i = 1; i <= n; i++) line = line "\t" ((i in row) ? row[i] : 0)
        print line
    }
    BEGIN {
        while ((getline name < samples) > 0) { n++; column[name] = n; header = header "\t" name }
        print "feature" header
    }
    NR > 1 && $2 != feature { flush(); split("", row) }
    { feature = $2; row[column[$1]] = $3 }
    END { if (NR > 0) flush() }
' > "$output"
