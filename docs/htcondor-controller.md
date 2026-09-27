# Run the Snakemake controller as an HTCondor job

Run `run-alpar.sh` on `conduit2.hpc.uni-saarland.de`:

```bash
bash ./run-alpar.sh --commit <commit-sha> automatix --latency-wait 120
```

The launcher exports the selected commit with `git archive`. Uncommitted changes
and untracked files are excluded. It saves a copy and checksum of the launcher,
the full workflow commit, Snakemake arguments, and the current environment's
pinned package versions in a new directory under `/home/joca00004/runs`.
It does not start tmux or modify the source checkout.

The controller requests 2 CPUs and 16 GiB RAM. The existing Snakemake profile
still controls resources for its child jobs. Local rules run inside the
controller's allocation. Set `ALPAR_CONTROLLER_MEMORY=32GB`, for example, to
request more controller memory on the next run.

The controller uses `/home/joca00004/.venvs/snakemake-htcondor` when its Python
interpreter and imports work in the container. Otherwise, it recreates the
pinned packages in `controller-env` under the run directory. This fallback needs
package download access. The default controller image is
`docker.io/python:3.11-bookworm`; set `ALPAR_CONTROLLER_IMAGE` to an image digest
to fix the container version as well as the workflow commit.

The worker wrapper creates its cache under the worker's `/tmp`. HTCondor can
forward the controller's scratch variables to child jobs, so workers must not
reuse the controller's cache path.
The controller passes a source cache under the private run directory to
Snakemake. Snakemake forwards that path to child jobs, which can then read and
write it from the shared home filesystem.
It also stores Snakemake's Apptainer prefix there because remote Snakemake
workers initialize that directory even for HTCondor Docker jobs.

## Authentication and startup

The launcher requests a two-day token limited to HTCondor READ and WRITE access.
The token stays in the private run directory on the shared home filesystem; it
is not written to logs or copied to this repository. The generated HTCondor
configuration points the Python plugin to the submission host and its collector.

Before starting ALPAR, the controller checks queue and history access, submits
one held probe, and removes it. The probe never runs on a worker. A failed check
stops startup with a nonzero exit status. If removal fails, `probe-job-id` records
the probe to inspect and remove from the submission host.

The token must remain valid while Snakemake runs. To renew it from conduit2
without displaying its contents (replace the run path):

```bash
run=/home/joca00004/runs/<run-directory>
(umask 077; condor_token_fetch -authz READ -authz WRITE -lifetime 172800 \
    -file "$run/controller.token.new") &&
    mv -- "$run/controller.token.new" "$run/auth/controller.token"
```

## Check a run

The launcher prints the controller's job ID and run directory. On conduit2:

```bash
condor_q <job-id>
tail -f /home/joca00004/runs/<run-directory>/controller.log
```

`logs/controller.<cluster>.log` records HTCondor events. Nonzero exits leave the
controller held for inspection. `exit-status` records the launcher's exit code
when it can finish, and `controller-usage.json` records Snakemake's peak RSS and
exit code. A forced kill may prevent those files from being written; use the
HTCondor event log in that case.

The controller prints a Python stack every three minutes. To request a stack
immediately, enter the running job with `condor_ssh_to_job <job-id>` and send
`SIGUSR1` to its Snakemake Python process. The signal prints a stack without
terminating Snakemake.

Pangenome clusters are processed in `pangenome_shards` fixed shards (config,
default 256), one job per shard and rule, with no checkpoint. The earlier
design created one job per cluster behind a checkpoint. With 34,148 clusters
and two antibiotics that is about 300,000 jobs. After the checkpoint, the
controller spent its time in Snakemake's `DAG.validate_group`, which searches
the upstream graph once per group, and its RSS grew past 6 GiB. The input
callbacks themselves took under a second. Shard counts above the number of
clusters are allowed; empty shards produce empty output directories. Raise
`--jobs` on the command line if more shards should run concurrently than the
profile's limit allows.

Snakemake keeps job metadata, including the markers for incomplete outputs, in
`<run-directory>/src/.snakemake`, which is private to each run. A new controller
therefore cannot tell that a killed job left a partial output behind, and will
treat that output as complete. Before relaunching after removing a controller,
wait until its children have left the queue. Then delete every path listed
under `<old-run>/src/.snakemake/incomplete`; the file names are base64-encoded
output paths (`for m in $(ls <old-run>/src/.snakemake/incomplete); do echo "$m" | base64 -d; echo; done`). When a child job itself is the only one
still running, as `binary_mutation_table` was on 2026-09-27, it is safer to let
it finish.

The selected commit's configuration determines input, output, and temporary
paths. The current configuration reuses `/home/joca00004/out`, so this resumes
existing results. Use distinct output and temporary paths in a committed config
for an independent experiment. Separate source exports alone do not isolate
results. Run only one controller against a given output directory, and inspect
remaining child jobs before restarting after a killed or evicted controller.

The previous launcher creates tmux sessions and detached Git worktrees. Existing
snapshots and results remain usable; this launcher does not delete them.
