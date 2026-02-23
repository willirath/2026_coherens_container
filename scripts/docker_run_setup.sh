#!/usr/bin/env bash
set -euo pipefail

SETUP="seamount_smoke"
if [ $# -gt 0 ]; then
  SETUP="$1"
  shift
fi
SETUP="${SETUP#setups/}"
IMAGE="${COHERENS_IMAGE:-coherens-local}"
RUN_ID="${COHERENS_RUN_ID:-$SETUP}"
ROOT="/opt/2026_coherens_container"
HOST_RUNS_DIR="${HOST_RUNS_DIR:-$(pwd)/runs}"
mkdir -p "$HOST_RUNS_DIR"

HOST_SETUP_DIR="$(pwd)/setups/$SETUP"

DOCKER_CMD=(
  docker run --rm
  -e "COHERENS_SETUP=$SETUP"
  -e "COHERENS_RUN_ID=$RUN_ID"
  -e "COHERENS_RESET=${COHERENS_RESET:-1}"
  -e "COHERENS_LAUNCH=${COHERENS_LAUNCH:-serial}"
  -e "COHERENS_MPI_PROCS=${COHERENS_MPI_PROCS:-4}"
  -e "COHERENS_MAKE_JOBS=${COHERENS_MAKE_JOBS:-1}"
  -v "$HOST_RUNS_DIR:/workspace/runs"
)

if [ -d "$HOST_SETUP_DIR" ]; then
  DOCKER_CMD+=(-v "$HOST_SETUP_DIR:$ROOT/coherens/setups/$SETUP")
fi

DOCKER_CMD+=("$IMAGE" "$SETUP")

exec "${DOCKER_CMD[@]}" "$@"
