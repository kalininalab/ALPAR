#!/bin/bash
set -euo pipefail

# Keep worker source caches in job-owned storage, not a shared /tmp/.cache.
export XDG_CACHE_HOME="${_CONDOR_SCRATCH_DIR:-$PWD}/.cache"

# The executor may retain the Python launcher before the Snakemake arguments.
if [[ "${2:-}" == "-m" && "${3:-}" == "snakemake" ]]; then
    shift 3
elif [[ "${1:-}" == "-m" && "${2:-}" == "snakemake" ]]; then
    shift 2
fi

# Prefer a rule-env Snakemake when present. Its shebang selects the Python used
# for script: rules; forcing /opt/snakemake/bin/snakemake would lose tool packages.
exec snakemake "$@"
