# syntax=docker/dockerfile:1.7

# Miniforge 26.7.2-0, pinned to its multi-platform manifest digest.
FROM condaforge/miniforge3:26.7.2-0@sha256:eeb947cc87d61d46820b123bd7c26e1cbdc4b182ff7d0331e501a32a936b82e3

ARG ENV_FILE
ARG ENV_NAME
ARG SNAKEMAKE_VERSION=9.26.1

# VCS-based pip dependencies in the PanPA environments need Git while building.
RUN apt-get update \
 && apt-get install --yes --no-install-recommends ca-certificates git \
 && rm -rf /var/lib/apt/lists/*

COPY ${ENV_FILE} /opt/environment.yaml

# The rule tools and the remote Snakemake worker live in separate environments.
RUN mamba env create \
      --yes \
      --prefix /opt/rule-env \
      --file /opt/environment.yaml \
 && mamba create \
      --yes \
      --prefix /opt/snakemake \
      --channel conda-forge \
      --channel bioconda \
      "snakemake-minimal=${SNAKEMAKE_VERSION}" \
 && mamba clean --all --yes

LABEL org.opencontainers.image.title="ALPAR ${ENV_NAME} rule environment" \
      org.opencontainers.image.source="https://github.com/kalininalab/ALPAR" \
      org.opencontainers.image.description="ALPAR rule environment from ${ENV_FILE}, with a separate Snakemake worker" \
      io.alpar.conda-environment="${ENV_NAME}" \
      io.alpar.snakemake-version="${SNAKEMAKE_VERSION}"

# No activation is needed at runtime. Snakemake's launcher retains its own Python
# shebang, while ordinary rule commands resolve from /opt/rule-env first.
ENV PATH="/opt/rule-env/bin:/opt/snakemake/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    HOME="/tmp" \
    TMPDIR="/tmp" \
    PYTHONNOUSERSITE="1"

WORKDIR /work

# HTCondor supplies the executable; an inherited entrypoint must not wrap it.
ENTRYPOINT []
CMD ["/bin/bash"]
