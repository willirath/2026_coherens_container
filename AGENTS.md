# Repository Guidelines

## Project Structure & Module Organization
- `coherens/`: git submodule containing the full COHERENS v3 source tree (physics, sediment, biology, tracers, utilities, and the official `setups/`). Treat this directory as upstream; never edit it directly without going through their repo.
- `docker/`: helper artifacts for containerized development (`docker_run_setup.sh` entrypoint and `coherensflags_docker.cmp`). Extend this area with additional tooling scripts or config snippets.
- `Dockerfile` & `.dockerignore`: build context for the runtime image. Local experiment outputs belong in `runs/` (ignored by Docker and Git); mount it into the container to persist builds.

## Build, Test, and Development Commands
- `docker build -t coherens-local .` — builds the Debian/gfortran image that ships the submodule and tooling.
- `docker run --rm -v $(pwd)/runs:/workspace/runs coherens-local tutorial/river` — stages and runs a bundled case inside the container, reusing cached artifacts from `runs/`.
- `docker run --rm -it -v $(pwd):/opt/2026_coherens_container coherens-local sedrough bash` — drop into a shell with your live checkout mounted; handy for debugging custom setups.
- Inside the container, the helper script runs `make linux-gfort` by default. Override with `-e COHERENS_MAKE_TARGET=linux-gfort-g` or `-e COHERENS_LAUNCH=mpi` as needed.

## Coding Style & Naming Conventions
- Keep shell helpers POSIX-compliant; use `set -euo pipefail` and provide clear environment variable knobs.
- When authoring supplemental Fortran modules, match upstream COHERENS style: two-space indentation, uppercase subroutine names in comments, and descriptive `USE` blocks grouped at the top.
- Follow existing naming patterns (`Usrdef_*` for user hooks, `COHERENS_*` paths for env vars) to avoid confusing the automation in `docker/docker_run_setup.sh`.

## Testing Guidelines
- Prefer running the smallest relevant setup (e.g., `tutorial/river`) before/after changes to verify the toolchain and physics compile.
- For MPI paths, exercise `docker run --rm -e COHERENS_LAUNCH=mpi -e COHERENS_MPI_PROCS=4 coherens-local <setup>` to ensure OpenMPI wiring remains intact.
- Capture runtime logs from `/workspace/runs/<setup>/` and attach excerpts to merge requests when behavior changes.

## Commit & Pull Request Guidelines
- Use short imperative commit subjects (`Add docker entrypoint override`, `Document custom setups`) and include context in the body when touching the submodule pointer or Docker image.
- Reference upstream issues/Zenodo docs when relevant, and list the tested command (docker build/run) in the PR description.
- For changes affecting container behavior, include before/after build or run output snippets so reviewers can confirm the workflow.
