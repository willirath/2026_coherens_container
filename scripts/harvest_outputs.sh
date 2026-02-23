#!/usr/bin/env bash
set -euo pipefail

RUNS_DIR="${RUNS_DIR:-$(pwd)/runs}"
OUTPUT_DIR="${OUTPUT_DIR:-$(pwd)/outputs}"
mkdir -p "$OUTPUT_DIR"

patterns=("*.nc" "*.tsout*" "*.log" "*.tst" "*.runlog" "*.errlog" "*.warlog" "*.inilog" "*.timing")
shopt -s nullglob
for dir in "$RUNS_DIR"/*; do
  [ -d "$dir" ] || continue
  setup=$(basename "$dir")
  dest="$OUTPUT_DIR/$setup"
  mkdir -p "$dest"
  moved=false
  for pat in "${patterns[@]}"; do
    for file in "$dir"/$pat; do
      mv "$file" "$dest/"
      moved=true
    done
  done
  if [ "$moved" = true ]; then
    echo "Saved outputs to $dest"
  fi
  moved=false
  unset moved
  if [ -z "$(ls -A "$dest")" ]; then
    rmdir "$dest" 2>/dev/null || true
  fi
  if [ -z "$(ls -A "$dir")" ]; then
    rmdir "$dir" 2>/dev/null || true
  fi
done
shopt -u nullglob
