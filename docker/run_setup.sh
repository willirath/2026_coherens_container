#!/usr/bin/env bash
set -euo pipefail

SETUP="${COHERENS_SETUP:-tutorial/river}"
if [ $# -gt 0 ]; then
  SETUP="$1"
  shift
fi
CUSTOM_CMD=("$@")

COHERENS_ROOT="${COHERENS_ROOT:-/opt/2026_coherens_container/coherens}"
FLAGS_FILE="${COHERENS_FLAGS_FILE:-/opt/2026_coherens_container/docker/coherensflags_docker.cmp}"
DEFAULT_RUN_ID="$(echo "$SETUP" | tr '/ ' '__')"
RUN_ID="${COHERENS_RUN_ID:-$DEFAULT_RUN_ID}"
WORKDIR_BASE="${COHERENS_WORKDIR_BASE:-/workspace/runs}"
WORKDIR="${COHERENS_WORKDIR:-$WORKDIR_BASE/$RUN_ID}"

if [ "${COHERENS_RESET:-0}" = "1" ] && [ -d "$WORKDIR" ]; then
  rm -rf "$WORKDIR"
fi
mkdir -p "$WORKDIR"
cd "$WORKDIR"

ln -snf "$COHERENS_ROOT" COHERENS

"$COHERENS_ROOT"/install_test -t "$SETUP" -o "$FLAGS_FILE"

if [ -n "${COHERENS_DEFRUNS_CONTENT:-}" ]; then
  printf "%s\n" "$COHERENS_DEFRUNS_CONTENT" > defruns
elif [ -n "${COHERENS_DEFRUNS_FILE:-}" ]; then
  cp "$COHERENS_DEFRUNS_FILE" defruns
fi

MAKE_TARGET="${COHERENS_MAKE_TARGET:-linux-gfort}"
MAKE_JOBS="${COHERENS_MAKE_JOBS:-$(nproc)}"
if [ "${COHERENS_SKIP_BUILD:-0}" != "1" ]; then
  make -j"$MAKE_JOBS" "$MAKE_TARGET"
fi

if [ ${#CUSTOM_CMD[@]} -gt 0 ]; then
  exec "${CUSTOM_CMD[@]}"
fi

if [ "${COHERENS_SKIP_RUN:-0}" = "1" ]; then
  exec "${SHELL:-/bin/bash}"
fi

if [ "${COHERENS_LAUNCH:-serial}" = "mpi" ]; then
  MPI_PROCS="${COHERENS_MPI_PROCS:-4}"
  exec mpirun -np "$MPI_PROCS" ./coherens
else
  exec ./Run
fi
