#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ENV_DIR="${REPO_ROOT}/snakefiles/envs"

IMAGE_NAMESPACE="${IMAGE_NAMESPACE:-}"
IMAGE_TAG="${IMAGE_TAG:-1.0.0}"
SNAKEMAKE_VERSION="${SNAKEMAKE_VERSION:-9.26.1}"

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

resolve_env_name() {
    local requested="$1"
    local env_name

    requested="$(basename "${requested}")"
    requested="${requested%.yaml}"
    if [[ "${requested}" == alpar-smk-* ]]; then
        env_name="${requested}"
    else
        env_name="alpar-smk-${requested}"
    fi

    [[ -f "${ENV_DIR}/${env_name}.yaml" ]] || die "environment not found: ${1}"
    printf '%s\n' "${env_name}"
}

test_one() {
    local env_name="$1"
    local image_repository image

    image_repository="${env_name}"
    if [[ -n "${IMAGE_NAMESPACE}" ]]; then
        image_repository="${IMAGE_NAMESPACE%/}/${image_repository}"
    fi
    image="${image_repository}:${IMAGE_TAG}"

    printf 'Testing %s as uid:gid 12345:12345\n' "${image}"
    docker image inspect "${image}" >/dev/null 2>&1 || die "image not found: ${image}"
    docker run --rm \
        --user 12345:12345 \
        --env HOME=/tmp \
        --env "EXPECTED_SNAKEMAKE_VERSION=${SNAKEMAKE_VERSION}" \
        "${image}" \
        bash -euc '
            test -r /opt/environment.yaml
            test -d /opt/rule-env
            test -x /opt/snakemake/bin/snakemake
            test "$(/opt/snakemake/bin/snakemake --version)" = "${EXPECTED_SNAKEMAKE_VERSION}"
            test "$(command -v snakemake)" = "/opt/snakemake/bin/snakemake"
            probe="${HOME}/alpar-container-write-test-${RANDOM}"
            touch "${probe}"
            rm "${probe}"
        '
}

command -v docker >/dev/null 2>&1 || die "docker is not installed or not on PATH"
docker info >/dev/null 2>&1 || die "the Docker daemon is not available"

env_names=()
if (( $# > 0 )); then
    for requested in "$@"; do
        env_names+=("$(resolve_env_name "${requested}")")
    done
else
    shopt -s nullglob
    for env_file in "${ENV_DIR}"/alpar-smk-*.yaml; do
        env_names+=("$(basename "${env_file}" .yaml)")
    done
    shopt -u nullglob
fi

(( ${#env_names[@]} > 0 )) || die "no environments selected"

for env_name in "${env_names[@]}"; do
    test_one "${env_name}"
done

printf 'All %d image(s) passed non-root runtime checks.\n' "${#env_names[@]}"
