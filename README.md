# COHERENS Container

This repository packages the COHERENS v3 ocean model in a Docker container so you can build and run setups without installing Fortran compilers on your host. The container handles compilation; your outputs land on the host via a mounted volume.

## Requirements
- Docker
- Git with submodule support

## Setup
```bash
git submodule update --init --recursive
docker build -t coherens-local .
```

## Run a built-in setup
```bash
docker run --rm \
  -v "$(pwd)/runs:/workspace/runs" \
  coherens-local seamount_smoke
```

Outputs land in runs/seamount_smoke/ on the host.

## Run again with a different ID (if the first run is still there)
```bash
docker run --rm \
  -v "$(pwd)/runs:/workspace/runs" \
  -e COHERENS_RUN_ID=seamount_smoke_2 \
  coherens-local seamount_smoke
```

## Bring your own setup
Copy a setup from the coherens submodule or use one from setups/:

```bash
docker run --rm \
  -v "$(pwd)/runs:/workspace/runs" \
  -v "$(pwd)/setups/my_case:/opt/coherens/setups/my_case" \
  coherens-local my_case
```

## Analyse outputs
The notebooks/ directory contains Jupyter notebooks for analysing results.
You need: xarray, matplotlib, numpy, netcdf4, notebook
```bash
pip install xarray matplotlib numpy netcdf4 notebook
jupyter notebook notebooks/
```

## Possible next steps
- Run with MPI: `-e COHERENS_LAUNCH=mpi -e COHERENS_MPI_PROCS=4`
