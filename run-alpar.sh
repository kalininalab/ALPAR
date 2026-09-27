#!/usr/bin/env bash
# Submit a foreground Snakemake controller from an exact workflow commit.
set -Eeuo pipefail

die() { printf 'error: %s\n' "$*" >&2; exit 1; }

if [[ ${1:-} == --inside-condor ]]; then
    run_dir=$(readlink -f -- "${2:?run directory required}")
    export HOME=/home/joca00004 USER=joca00004 LOGNAME=joca00004
    export XDG_CACHE_HOME="${_CONDOR_SCRATCH_DIR:?not running inside HTCondor}/.cache"
    export CONDOR_CONFIG="$run_dir/condor.conf"
    export _CONDOR_SEC_TOKEN_DIRECTORY="$run_dir/auth"
    # The plugin creates Schedd() clients for submission, history and cancellation.
    # All of them must use the submission host, rather than the execution host.
    unset _CONDOR_SCHEDD_ADDRESS_FILE
    cd "$run_dir/src"
    exec > >(tee -a "$run_dir/controller.log") 2>&1
    trap 'result=$?; printf "%s\n" "$result" > "$run_dir/exit-status"' EXIT
    printf 'Controller host: %s\nWorkflow commit: %s\n' "$(hostname)" "$(cat "$run_dir/commit")"
    controller_python=$(cat "$run_dir/source-python")
    if "$controller_python" -c 'import snakemake, htcondor2' 2>/dev/null; then
        printf 'Using existing environment: %s\n' "$controller_python"
    else
        # A venv's interpreter may be absent or incompatible in the container.
        # Recreate its pinned packages on shared storage if needed.
        controller_python="$run_dir/controller-env/bin/python"
        if [[ ! -f "$run_dir/controller-env/ready" ]]; then
            python3 -m venv "$run_dir/controller-env"
            "$controller_python" -m pip install --disable-pip-version-check \
                -r "$run_dir/requirements.txt"
            touch "$run_dir/controller-env/ready"
        fi
    fi
    export PATH="$(dirname "$controller_python"):$PATH"
    "$controller_python" -u - "$run_dir" <<'PY'
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys

import htcondor2 as htcondor

run = Path(sys.argv[1])
schedd = htcondor.Schedd()
constraint = 'Owner == "joca00004"'
schedd.query(constraint=constraint, projection=["ClusterId"], limit=1)
list(schedd.history(constraint, ["ClusterId"], 1))
# A held probe verifies WRITE and cancellation without consuming a worker slot.
probe = htcondor.Submit({
    "universe": "docker",
    "docker_image": "docker.io/cambouu/alpar-smk-python313:1.0.0",
    "executable": "/bin/true",
    "transfer_executable": "False",
    "hold": "True",
    "request_cpus": "1",
    "request_memory": "128MB",
    "+ALPARControllerProbe": "True",
})
probe_id = schedd.submit(probe).cluster()
(run / "probe-job-id").write_text(f"{probe_id}.0\n")
schedd.act(htcondor.JobAction.Remove, [f"{probe_id}.0"])
print("Scheduler query, history, submission and cancellation checks passed.", flush=True)

args = json.loads((run / "snakemake-args.json").read_text())
source_cache = run / "runtime-source-cache"
source_cache.mkdir(exist_ok=True)
deployment_cache = run / "apptainer-cache"
deployment_cache.mkdir(exist_ok=True)
bootstrap = """
import faulthandler
import signal
import sys
from snakemake.cli import main

faulthandler.enable()
faulthandler.register(signal.SIGUSR1, all_threads=True)
faulthandler.dump_traceback_later(180, repeat=True)
sys.argv[0] = "snakemake"
sys.exit(main())
"""
command = [sys.executable, "-u", "-c", bootstrap,
           "--profile", str(run / "src/snakefiles/profiles/htcondor-containers/profile.v9+.yaml"),
           "--runtime-source-cache-path", str(source_cache),
           "--apptainer-prefix", str(deployment_cache),
           "--rerun-incomplete", *args]
print("Starting:", " ".join(command), flush=True)
child = subprocess.Popen(command)

def forward_signal(signum, frame):
    child.send_signal(signal.SIGINT)

signal.signal(signal.SIGTERM, forward_signal)
signal.signal(signal.SIGINT, forward_signal)
result = child.wait()
peak_mib = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
(run / "controller-usage.json").write_text(json.dumps({
    "snakemake_exit_code": result, "snakemake_peak_rss_mib": peak_mib,
}, indent=2) + "\n")
print(f"Snakemake exit={result}; peak RSS={peak_mib:.1f} MiB", flush=True)
sys.exit(result if result >= 0 else 128 - result)
PY
    exit 0
fi

repo=${ALPAR_REPO:-/home/joca00004/ALPAR}
runs_dir=${ALPAR_RUNS_DIR:-/home/joca00004/runs}
snakemake_python=${ALPAR_SNAKEMAKE_PYTHON:-/home/joca00004/.venvs/snakemake-htcondor/bin/python}
controller_image=${ALPAR_CONTROLLER_IMAGE:-docker.io/python:3.11-bookworm}
controller_memory=${ALPAR_CONTROLLER_MEMORY:-16GB}
# Snakemake state (metadata, incomplete markers, locks) shared by every run
# that writes to the same output directory.
state_dir=${ALPAR_STATE_DIR:-/home/joca00004/snakemake-state/out}
commit_ref=HEAD

usage() {
    cat <<'EOF'
Usage: ./run-alpar.sh [--commit COMMIT] [Snakemake target and options ...]
Example: ./run-alpar.sh --commit 7e1d1464 automatix --latency-wait 120

Exports the selected commit and submits one HTCondor controller (2 CPUs, 16 GiB).
The default target is automatix. Uncommitted and untracked files are excluded.
Uses the selected commit's config, including its output paths.
Snakemake state is shared across runs through ALPAR_STATE_DIR, so metadata and
incomplete markers survive a restart and the lock prevents concurrent controllers.
ALPAR_REPO, ALPAR_RUNS_DIR, ALPAR_STATE_DIR, ALPAR_CONTROLLER_MEMORY and
ALPAR_CONTROLLER_IMAGE can override the defaults. The controller token is valid for two days.
EOF
}

case ${1:-} in
    --help|-h) usage; exit 0 ;;
    --commit) commit_ref=${2:?a commit is required}; shift 2 ;;
esac
(( $# )) || set -- automatix
for program in git condor_submit condor_token_fetch condor_config_val timeout tar; do
    command -v "$program" >/dev/null || die "required command not found: $program"
done
[[ -x $snakemake_python ]] || die "controller Python not found: $snakemake_python"
commit=$(git -C "$repo" rev-parse --verify "${commit_ref}^{commit}")
if [[ -n $(git -C "$repo" status --porcelain) ]]; then
    printf 'Using commit %s; local changes are excluded.\n' "$commit" >&2
fi
mkdir -p "$runs_dir"
# ALPAR_PREPARED_RUN is for completing a run whose token was prepared separately.
if [[ -n ${ALPAR_PREPARED_RUN:-} ]]; then
    run_dir=$(readlink -f -- "$ALPAR_PREPARED_RUN")
    [[ -d $run_dir && ! -e $run_dir/src ]] || die "prepared run must exist and have no src directory"
else
    run_dir=$(mktemp -d "$runs_dir/$(date -u +%Y%m%dT%H%M%SZ)-${commit:0:12}-XXXXXX")
fi
mkdir "$run_dir/src" "$run_dir/logs"
mkdir -p -m 700 "$run_dir/auth"
git -C "$repo" archive "$commit" | tar -xf - -C "$run_dir/src"
# Without shared state, each run starts with empty metadata: Snakemake 9.27 then
# reruns every script rule ("Code has changed") and cannot see partial outputs
# left by a killed controller. Adopt metadata from older private-state runs,
# oldest first, keeping the newest record for each output.
mkdir -p "$state_dir/metadata"
for legacy in "$runs_dir"/*/src/.snakemake/metadata; do
    [[ -d $legacy && ! -L ${legacy%/metadata} ]] || continue
    find "$legacy" -maxdepth 1 -type f -exec cp -pu -t "$state_dir/metadata" {} +
done
ln -s "$state_dir" "$run_dir/src/.snakemake"
printf '%s\n' "$commit" > "$run_dir/commit"
cp -- "$(readlink -f -- "$0")" "$run_dir/run-alpar.sh"
sha256sum "$run_dir/run-alpar.sh" > "$run_dir/launcher.sha256"
printf '%s\n' "$snakemake_python" > "$run_dir/source-python"
"$snakemake_python" -m pip freeze > "$run_dir/requirements.txt"
"$snakemake_python" -c 'import json,sys; from pathlib import Path; Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:]) + "\n")' \
    "$run_dir/snakemake-args.json" "$@"
if [[ ! -s "$run_dir/auth/controller.token" ]]; then
    (umask 077; timeout 30 condor_token_fetch -authz READ -authz WRITE \
        -lifetime 172800 -file "$run_dir/auth/controller.token")
fi
chmod 600 "$run_dir/auth/controller.token"
collector=$(condor_config_val COLLECTOR_HOST)
schedd_host=$(hostname -f)
cat > "$run_dir/condor.conf" <<EOF
COLLECTOR_HOST = $collector
SCHEDD_HOST = $schedd_host
UID_DOMAIN = cs.uni-saarland.de
SEC_CLIENT_AUTHENTICATION_METHODS = IDTOKENS
SEC_DEFAULT_AUTHENTICATION = REQUIRED
SEC_DEFAULT_ENCRYPTION = REQUIRED
SEC_TOKEN_DIRECTORY = $run_dir/auth
EOF
cat > "$run_dir/controller.sub" <<EOF
universe = docker
docker_image = $controller_image
executable = /bin/bash
transfer_executable = False
arguments = "$run_dir/run-alpar.sh --inside-condor $run_dir"
initialdir = $run_dir
output = $run_dir/logs/controller.\$(ClusterId).out
error = $run_dir/logs/controller.\$(ClusterId).err
log = $run_dir/logs/controller.\$(ClusterId).log
stream_output = True
stream_error = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = ""
request_cpus = 2
request_memory = $controller_memory
request_disk = 4GB
requirements = UidDomain == "cs.uni-saarland.de"
+WantGPUHomeMounted = true
+ALPARCommit = "$commit"
+ALPARController = true
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
queue 1
EOF
printf 'Run: %s\nCommit: %s\n' "$run_dir" "$commit"
condor_submit -dry-run "$run_dir/controller.classad" "$run_dir/controller.sub"
condor_submit -terse "$run_dir/controller.sub" | tee "$run_dir/job-id"
printf 'Controller log: %s/controller.log\nExit status: %s/exit-status\n' "$run_dir" "$run_dir"
