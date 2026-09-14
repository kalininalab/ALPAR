#!/bin/bash
set -euo pipefail

# Keep worker source caches in job-owned storage, not a shared /tmp/.cache.
export XDG_CACHE_HOME="${_CONDOR_SCRATCH_DIR:-$PWD}/.cache"

# Prefer a rule-env Snakemake when present. Its shebang selects the Python used
# for script: rules; forcing /opt/snakemake/bin/snakemake would lose tool packages.
exec snakemake "$@"
