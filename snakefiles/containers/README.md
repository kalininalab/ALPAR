# ALPAR rule-environment images

Each `snakefiles/envs/alpar-smk-*.yaml` file is build-time input for one image.
The Conda environment is installed at `/opt/rule-env`; an independent
Snakemake worker is installed at `/opt/snakemake`. Runtime activation and nested
Conda/Apptainer execution are intentionally unnecessary.

The current inventory contains 17 images. The workflow currently references 16
of them; `alpar-smk-panpa.yaml` is retained and built as a legacy environment.

## Build and test locally

Run commands from the repository root in WSL:

```bash
cd ~/ALPAR

# Start with one small image.
snakefiles/containers/build-env-images.sh bubblegun
snakefiles/containers/test-env-images.sh bubblegun

# Build and test every alpar-smk-*.yaml environment.
snakefiles/containers/build-env-images.sh
snakefiles/containers/test-env-images.sh
```

Images default to `linux/amd64`, Snakemake `9.26.1`, and tag `1.0.0`, for
example `alpar-smk-bubblegun:1.0.0`. Override these settings with environment
variables:

```bash
PLATFORM=linux/amd64 \
SNAKEMAKE_VERSION=9.26.1 \
IMAGE_TAG=1.0.0 \
snakefiles/containers/build-env-images.sh
```

Pass one or more environment keys to build or test a subset:

```bash
snakefiles/containers/build-env-images.sh mafft prokka snippy
snakefiles/containers/test-env-images.sh mafft prokka snippy
```

## Publish to Docker Hub

Log in, choose the Docker Hub namespace, and push only after local tests pass:

```bash
docker login

IMAGE_NAMESPACE=docker.io/DOCKERHUB_USERNAME \
IMAGE_TAG=1.0.0 \
PUSH=1 \
snakefiles/containers/build-env-images.sh
```

The namespace is part of the build tag, so names become, for example,
`docker.io/DOCKERHUB_USERNAME/alpar-smk-bubblegun:1.0.0`.

Record the immutable digest after each push:

```bash
docker buildx imagetools inspect \
  docker.io/DOCKERHUB_USERNAME/alpar-smk-bubblegun:1.0.0
```

Use `docker://docker.io/DOCKERHUB_USERNAME/alpar-smk-bubblegun@sha256:...`
in the source profile to pin an immutable image rather than a mutable tag.


## Pangenome job groups

Run the workflow from the repository root with the checked-in profile:

```bash
snakemake pangenome --profile snakefiles/profiles/htcondor-containers
```

Pangenome alignment, graph construction, BubbleGun, training features, test
alignment, and test features each have their own group. The profile bundles up
to 2000 independent jobs of each rule per submission, like SNP and Prokka.
`cores` and `mem_mb` bound concurrency inside the group; `group-components`
controls how many jobs share a submission. Each group uses one environment
image. Local execution uses ordinary per-file jobs.

The cluster-splitting checkpoint is the only dynamic discovery boundary.
Cluster filenames are cached after one directory scan and refreshed if the
checkpoint directory changes. All subsequent jobs use those filenames as
wildcards, including PanPA graph construction. PanPA's index reads an explicit
alignment inventory, so unrelated or stale files in the alignment folder are
excluded.

Per-cluster feature files now live under
`pangenome/bubble_features/{train,test}/{antibiotic}/`. LOR lookups and test GAFs
are written directly to their final locations. The merged train/test TSV paths
and the `pangenome` completion target are unchanged. Table merging streams the
exact declared inputs and replaces previous contents on rerun.

Existing batch folders are no longer inputs and are not deleted automatically.
A resumed run builds missing per-file outputs; old batch directories are not
imported. Custom profiles referring to `batch_*` or `batched_*` pangenome rules
must use the new rule names in the checked-in profile. `panpa_build_gfa` now
uses one thread per graph instead of parallelizing a whole directory internally.

## HTCondor clusters requiring local SIF images

Some clusters reject registry URLs in `container_image` and require a full path
to a local `.sif` file. The checked-in `htcondor-containers` profile is the source
inventory of registry images. Prepare a local profile before submitting on these
clusters:

```bash
cd ~/ALPAR
# Use the Python environment containing Snakemake and PyYAML.
# Apptainer or Singularity must also be on PATH for image pulls.
python snakefiles/containers/prepare-htcondor-sif.py --image-dir "$HOME/alpar-images"
snakemake pangenome --profile .snakemake/htcondor-sif-profile/
```

Choose an image directory allowed by your cluster and accessible to the scheduler
(for example, your group's scratch directory). Apptainer's `pull` converts each
Docker image to SIF. The preparation command downloads each distinct image once,
reuses existing nonempty files, and writes the profile only after all pulls succeed.
It preserves scheduling settings, sets absolute image and wrapper paths, and
resolves the configfile relative to the source profile. Rerun preparation after
changing image references or other settings in the source profile.

Use `--dry-run` to preview pulls without downloading or writing files, or
`--profile-dir PATH` to choose where to write the generated profile. The default
output is ignored by Git under `.snakemake/`. HTCondor still owns container
execution; keep Snakemake's Conda/Apptainer deployment disabled with this profile.
