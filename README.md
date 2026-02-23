# COHERENS Quick Runner

This repository wraps the official `coherens_v3` source tree in a reproducible workflow:

1. Build a Linux container with the full Fortran/MPI toolchain.
2. Stage and run COHERENS setups (stock or custom copies) inside the container while storing artifacts under `runs/` on the host.
3. Analyse NetCDF outputs locally using a Pixi-managed Python environment (xarray, matplotlib, Jupyter).

## Requirements
- Docker Desktop (macOS/Windows) or Docker Engine (Linux).
- [Pixi](https://pixi.sh) for the lightweight analysis environment.
- Git with submodule support (`git submodule update --init --recursive`).

## Repository Layout
```
.
├── coherens/            # upstream submodule (leave untouched)
├── docker/              # entrypoint script plus Dockerfile
├── notebooks/           # Jupyter analysis (e.g., seamount_demo.ipynb)
├── runs/                # host-mounted build directories (ignored by Git)
├── outputs/             # auto-harvested outputs (via pixi run harvest)
├── scripts/             # helper shell wrappers
├── setups/              # curated copies of upstream test cases
│   ├── seamount_smoke/  # 1-hour validation case
│   └── seamount_demo/   # 12-hour demo with richer diagnostics
└── pixi.toml            # analysis environment + task definitions
```

## Quick Start
```bash
# Install Pixi dependencies (xarray, matplotlib, notebook, netCDF4)
pixi install

# Build the Docker image containing the COHERENS stack
pixi run build-image

# Run any setup directory directly (this example writes to runs/seamount_demo)
pixi run experiment setups/seamount_demo

# Collect NetCDF/log files under outputs/<setup>/
pixi run harvest
```

All runs share the same [runs/](runs) directory on the host. Every `pixi run experiment …` invocation starts from a clean slate (`COHERENS_RESET=1`). If you need to keep an existing run directory for inspection, override with `COHERENS_RESET=0 pixi run experiment …`.

## Pixi Tasks
| Task | Description |
| --- | --- |
| `build-image` | `docker build -t coherens-local .` |
| `experiment` | Accepts a setup path (`setups/<case>`) and runs it inside the container (default `setups/seamount_smoke`). |
| `harvest` | Moves `*.nc`, `*.log`, etc. from `runs/<setup>/` into `outputs/<setup>/`. |

Override values at call time via env vars when needed, for example:

```bash
COHERENS_MAKE_JOBS=4 pixi run experiment setups/my_case
```

## Example Workflows
### Smoke Test (fast validation)
- Definition: copy of seamount case trimmed to a single hour, only `seamountA` stage ([setups/seamount_smoke](setups/seamount_smoke)).
- Command:

  ```bash
  pixi run experiment setups/seamount_smoke
  ```

- Outputs: `seamountA_1.tsout3.nc` (3-D fields), `seamountA_2.tsout0.nc` (time series), `.log/.tst` files under [runs/seamount_smoke](runs/seamount_smoke).

### Seamount Demo (richer plotting)
- Definition: 12-hour run with the same grid and physics but more time samples ([setups/seamount_demo](setups/seamount_demo)).
- Command:

  ```bash
  pixi run experiment setups/seamount_demo
  ```

- Outputs: same file structure under [runs/seamount_demo](runs/seamount_demo); suitable for time animations or vertical profiles.

After a run, collect the NetCDF/log files with `pixi run harvest`. Each setup gets its own folder under [outputs/](outputs). Runs reset automatically; set `COHERENS_RESET=0` if you need to keep intermediate build folders.

## Output Files
After harvesting, you will typically see several NetCDF products (for example, under [outputs/seamount_demo](outputs/seamount_demo)):

- `*_1.tsout3.nc` – full 3-D time series (x, y, z, time) of fields such as temperature and velocity.
- `*_1.tsout2.nc` – depth-averaged 2-D fields from the same run.
- `*_2.tsout0.nc` – point or domain-integrated diagnostics (kinetic energy, mean speeds, etc.).

## Plotting & Notebooks
Once a run completes (and you've harvested into, e.g., [outputs/seamount_demo](outputs/seamount_demo)), launch Jupyter manually inside the Pixi environment:
```bash
pixi run python -m notebook notebooks/seamount_demo.ipynb
```
The provided notebook reads data from [outputs/seamount_demo](outputs/seamount_demo) and plots simple diagnostics. Feel free to duplicate it for other setups.

## Custom Setups
- Copy any upstream case:

  ```bash
  cp -R coherens/setups/<name> setups/<my_case>
  ```
- Adjust `defruns`, `cifruns`, and `Usrdef_*` routines as needed (e.g., shorten `CEndDateTime`).
- Run it with `pixi run experiment setups/<my_case>`; the wrapper will mount that folder into the container automatically.
