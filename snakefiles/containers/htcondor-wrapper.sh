#!/bin/bash
set -euo pipefail

# Controller scratch variables can be forwarded to child jobs by the executor.
# Allocate cache and temp paths inside this worker's own container instead.
export TMPDIR=/tmp
export XDG_CACHE_HOME="$(mktemp -d /tmp/alpar-worker-cache.XXXXXXXX)"

export USER="${USER:-joca00004}"
export LOGNAME="${LOGNAME:-joca00004}"

# The executor may retain the Python launcher before the Snakemake arguments.
if [[ "${2:-}" == "-m" && "${3:-}" == "snakemake" ]]; then
    shift 3
elif [[ "${1:-}" == "-m" && "${2:-}" == "snakemake" ]]; then
    shift 2
fi

# Prefer a rule-env Snakemake when present. Its shebang selects the Python used
# for script: rules; forcing /opt/snakemake/bin/snakemake would lose tool packages.
exec snakemake "$@"
