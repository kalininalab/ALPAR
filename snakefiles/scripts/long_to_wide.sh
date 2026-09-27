#!/bin/bash
# Pivot a headerless long table (sample, feature, value) into a wide table.
#
# ROW_FIELD and COLUMN_FIELD select which of the first two fields become rows
# and columns; the third field supplies the cell values, and absent pairs get
# FILL. Duplicate pairs keep their last value. This is equivalent to
#   mlr --tsv --implicit-tsv-header label hash,feature,value \
#     then reshape -s <column>,value then unsparsify --fill-with <FILL>
# except that rows are sorted by key and columns follow their first appearance
# in the input. Miller holds every record in memory; this streams after an
# external sort, so memory grows only with the number of distinct columns.
#
# Usage: long_to_wide.sh INPUT OUTPUT TMPDIR THREADS SORT_MEMORY ROW_FIELD COLUMN_FIELD ROW_HEADER [FILL]
set -euo pipefail

input=$1 output=$2 tmp_root=$3 threads=$4 sort_memory=$5
row_field=$6 column_field=$7 row_header=$8 fill=${9:-}
export LC_ALL=C

work=$(mktemp -d "$tmp_root/long-to-wide.XXXXXX")
trap 'rm -rf "$work"' EXIT

awk -F '\t' -v key="$column_field" '!seen[$key]++ { print $key }' "$input" > "$work/columns"

sort -t "$(printf '\t')" -k"$row_field","$row_field" -s -S "$sort_memory" \
    --parallel="$threads" -T "$work" "$input" |
awk -F '\t' -v columns="$work/columns" -v row_key="$row_field" -v column_key="$column_field" \
    -v row_header="$row_header" -v fill="$fill" '
    # Print field by field: concatenating rows with millions of columns is quadratic.
    function flush(    i) {
        printf "%s", current
        for (i = 1; i <= n; i++) printf "\t%s", ((i in row) ? row[i] : fill)
        printf "\n"
    }
    BEGIN {
        printf "%s", row_header
        while ((getline name < columns) > 0) { n++; index_of[name] = n; printf "\t%s", name }
        printf "\n"
    }
    NR > 1 && $row_key != current { flush(); split("", row) }
    { current = $row_key; row[index_of[$column_key]] = $3 }
    END { if (NR > 0) flush() }
' > "$output"
