#!/bin/bash
# Run a bash snippet for every cluster in a shard manifest, THREADS at a time.
#
# Usage: for_each_cluster.sh MANIFEST THREADS LOG SNIPPET
#
# SNIPPET runs as `bash -e -o pipefail -c SNIPPET _ CLUSTER TOOL_LOG`. $1 is the
# cluster filename and $2 a scratch path the snippet may pass as a tool's log
# file. The snippet's output and the tool log are written to LOG per cluster,
# in manifest order. Every cluster runs even if others fail; the script then
# exits non-zero if any failed.
set -euo pipefail

manifest=$1 threads=$2 log=$3 snippet=$4
logs=$(mktemp -d)
trap 'rm -rf "$logs"' EXIT
export ALPAR_CLUSTER_SNIPPET=$snippet ALPAR_CLUSTER_LOGS=$logs

status=0
# A worker exit status of 1 lets xargs continue with other clusters; it then exits 123.
xargs -r -d '\n' -P "$threads" -n 1 bash -c '
    cluster=$1 out="$ALPAR_CLUSTER_LOGS/$1.log" tool="$ALPAR_CLUSTER_LOGS/$1.tool"
    bash -e -o pipefail -c "$ALPAR_CLUSTER_SNIPPET" _ "$cluster" "$tool" > "$out" 2>&1 || failed=1
    if [ -f "$tool" ]; then cat "$tool" >> "$out"; fi
    if [ -n "${failed:-}" ]; then echo "FAILED: $cluster" >> "$out"; exit 1; fi
' _ < "$manifest" || status=$?

: > "$log"
while IFS= read -r cluster; do
    if [ -f "$logs/$cluster.log" ]; then
        printf '>> %s\n' "$cluster" >> "$log"
        cat "$logs/$cluster.log" >> "$log"
    fi
done < "$manifest"
if [ "$status" -ne 0 ]; then
    echo "$(grep -c '^FAILED: ' "$log") cluster(s) failed; see FAILED lines above." >> "$log"
fi
exit "$status"
