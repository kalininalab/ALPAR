#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ENV_DIR="${REPO_ROOT}/snakefiles/envs"
DOCKERFILE="${SCRIPT_DIR}/conda-env.Dockerfile"

IMAGE_NAMESPACE="${IMAGE_NAMESPACE:-}"
IMAGE_TAG="${IMAGE_TAG:-1.0.0}"
PLATFORM="${PLATFORM:-linux/amd64}"
SNAKEMAKE_VERSION="${SNAKEMAKE_VERSION:-9.26.1}"
PUSH="${PUSH:-0}"

usage() {
    cat <<'EOF'
Build one Docker image for each ALPAR rule environment.

Usage:
  snakefiles/containers/build-env-images.sh [ENV ...]

ENV may be a short key such as "mafft", an image/environment name such as
"alpar-smk-mafft", or an environment YAML path. With no arguments, all files
matching snakefiles/envs/alpar-smk-*.yaml are built.

Configuration is supplied through environment variables:
  IMAGE_NAMESPACE     Optional registry namespace, e.g. docker.io/my-user
  IMAGE_TAG           Image tag (default: 1.0.0)
  PLATFORM            Target platform (default: linux/amd64)
  SNAKEMAKE_VERSION   Worker version (default: 9.26.1)
  PUSH                Set to 1 to push each successful build (default: 0)

Examples:
  snakefiles/containers/build-env-images.sh bubblegun
  IMAGE_NAMESPACE=docker.io/my-user snakefiles/containers/build-env-images.sh
  IMAGE_NAMESPACE=docker.io/my-user PUSH=1 snakefiles/containers/build-env-images.sh
EOF
}

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

resolve_env_file() {
    local requested="$1"
    local candidate

    if [[ -f "${requested}" ]]; then
        candidate="$(cd "$(dirname "${requested}")" && pwd)/$(basename "${requested}")"
    else
        requested="${requested%.yaml}"
        requested="${requested#alpar-smk-}"
        candidate="${ENV_DIR}/alpar-smk-${requested}.yaml"
    fi

    [[ -f "${candidate}" ]] || die "environment not found: ${1}"
    case "${candidate}" in
        "${ENV_DIR}"/alpar-smk-*.yaml) printf '%s\n' "${candidate}" ;;
        *) die "environment must match ${ENV_DIR}/alpar-smk-*.yaml: ${candidate}" ;;
    esac
}

build_one() {
    local env_file="$1"
    local env_basename env_name image_repository image

    env_basename="$(basename "${env_file}")"
    env_name="${env_basename%.yaml}"
    image_repository="${env_name}"
    if [[ -n "${IMAGE_NAMESPACE}" ]]; then
        image_repository="${IMAGE_NAMESPACE%/}/${image_repository}"
    fi
    image="${image_repository}:${IMAGE_TAG}"

    printf '\nBuilding %s from %s\n' "${image}" "${env_file#"${REPO_ROOT}/"}"
    docker build \
        --platform "${PLATFORM}" \
        --build-arg "ENV_FILE=${env_file#"${REPO_ROOT}/"}" \
        --build-arg "ENV_NAME=${env_name}" \
        --build-arg "SNAKEMAKE_VERSION=${SNAKEMAKE_VERSION}" \
        --file "${DOCKERFILE}" \
        --tag "${image}" \
        "${REPO_ROOT}"

    if [[ "${PUSH}" == "1" ]]; then
        [[ -n "${IMAGE_NAMESPACE}" ]] || die "IMAGE_NAMESPACE is required when PUSH=1"
        docker push "${image}"
    fi
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    usage
    exit 0
fi

command -v docker >/dev/null 2>&1 || die "docker is not installed or not on PATH"
docker info >/dev/null 2>&1 || die "the Docker daemon is not available"
[[ -f "${DOCKERFILE}" ]] || die "Dockerfile not found: ${DOCKERFILE}"
[[ "${PUSH}" == "0" || "${PUSH}" == "1" ]] || die "PUSH must be 0 or 1"

env_files=()
if (( $# > 0 )); then
    for requested in "$@"; do
        env_files+=("$(resolve_env_file "${requested}")")
    done
else
    shopt -s nullglob
    env_files=("${ENV_DIR}"/alpar-smk-*.yaml)
    shopt -u nullglob
fi

(( ${#env_files[@]} > 0 )) || die "no alpar-smk-*.yaml files found in ${ENV_DIR}"

printf 'Building %d image(s) for %s with Snakemake %s\n' \
    "${#env_files[@]}" "${PLATFORM}" "${SNAKEMAKE_VERSION}"
for env_file in "${env_files[@]}"; do
    build_one "${env_file}"
done

printf '\nAll %d image(s) built successfully.\n' "${#env_files[@]}"
