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
in production HTCondor configuration rather than a mutable tag.
