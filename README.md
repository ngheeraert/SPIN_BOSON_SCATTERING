# Multi-polaron waveguide-QED dynamics (Fortran)

This repository contains a Fortran implementation of time-dependent variational dynamics for a two-level system coupled to a 1D bosonic continuum (spin–boson / waveguide-QED setting). The field is represented using a **superposition of multimode coherent states** (often called a *multi-polaron* / *MCS* ansatz). The code evolves the variational parameters in time, produces text data files (`.d`) for observables, and supports parameter sweeps via the provided cluster submission scripts.

## Contents

- `main.f90` — program entry point: parses CLI arguments, runs one trajectory.
- `systm.f90` — core types (`param`, `state`, `traj`) + physics/numerics routines (EOM, overlaps, observables).
- `output.f90` — time-integration driver + output writers (time traces, spectra, correlation functions).
- `typedefs.f90` — numerical kinds (single/double precision selection via preprocessing).
- `consts.f90` — constants and small helpers.
- `lapackmodule.f90`, `inverse.f90` — LAPACK/BLAS helpers used by the linear algebra in the EOM.
- `bash_scat.sh`, `bash_1scat.sh`, `gs.sh` — example sweep / submission scripts for `qsub`-style clusters.

> **Data format:** outputs are plain text with whitespace-separated columns, suitable for post-processing in Python/Matlab/gnuplot.

---

## Requirements

### Compiler
- A Fortran compiler supporting free-form Fortran 90/95 (e.g. `gfortran`, `ifort`).

### Libraries
- BLAS + LAPACK (system libraries, OpenBLAS, MKL, etc.).

### (Optional) Cluster utilities
The provided bash scripts assume:
- `qsub` is available (PBS/Torque-like scheduler),
- `bc` and `python` exist on the submission node,
- a cluster-specific `jobfile` template exists (not included here).

---

## Building

### Recommended: double precision
Precision is selected in `typedefs.f90` via preprocessing:
- `-DDP` enables double precision (`r_type = dp`).

**Example (gfortran + system BLAS/LAPACK):**
```bash
gfortran -O3 -cpp -DDP   typedefs.f90 consts.f90 lapackmodule.f90 inverse.f90 systm.f90 output.f90 main.f90   -llapack -lblas -o mpol
```

**Example (ifort + MKL, illustrative):**
```bash
ifort -O3 -fpp -DDP   typedefs.f90 consts.f90 lapackmodule.f90 inverse.f90 systm.f90 output.f90 main.f90   -mkl -o mpol
```

---

## Running

`main.f90` calls `getParameters(sys)` which parses command-line flags. The most commonly used flags are:

### Physical / model parameters
- `-al <alpha>`   : dimensionless coupling strength.
- `-del <Delta>`  : bare qubit splitting.
- `-wc <wc>`      : cutoff frequency.
- `-wmax <wmax>`  : maximum retained frequency (band edge / smooth cutoff).
- `-nm <nmode>`   : number of discretized bosonic modes.

### Wavepacket / drive parameters
- `-k0 <k0>`      : central wavevector / frequency of the incident packet.
- `-k0_2 <k0_2>`  : secondary frequency (used in a two-tone preparation mode).
- `-x0 <x0>`      : initial packet center (often negative to start left of the scatterer).
- `-sigma <sigma>`: spectral width of the packet.
- `-n <n_wp>`     : intensity parameter (used for coherent-state preparation modes).
- `-xmin <xmin>`  : spatial cutoff used in some Fourier/field post-processing.
- `-xg2 <xg2>`    : reference position used in some correlation outputs.

### Time stepping
- `-dt <dt>`      : integration time step.
- `-tmax <tmax>`  : final time.

### Variational basis size (coherent components)
- `-npini <np>`   : initial number of coherent components (“polarons”).
- `-npadd <n>`    : number of additional components to add during evolution (0 disables adding).
- `-p0 <p0>`      : initial amplitude assigned to newly added components.
- `-tref <tref>`  : reference time used in the adding schedule.
- `-me <merr>`    : error threshold used to decide when to add components.
- `-mde <max_deltaE>` : energy-drift threshold used in some stability checks.
- `-bt <0/1>`     : enable/disable “back-tracking” logic in the adaptive procedure (when available).

### Preparation mode (`-prep`)
`-prep` selects how the initial state is prepared inside `output:printTrajectory_DL`:
- `50` : coherent-state packet + dressed ground-state background (requires `gs_data/`).
- `51` : single-photon-like packet + dressed ground-state background (requires `gs_data/`).
- `52` : single-photon-like packet without ground-state subtraction.
- `54` : two-frequency coherent packet + dressed ground-state background (requires `gs_data/`).
- `60` : restart mode (loads a saved state from file after setting up reference packet metadata).

*(Other values exist in the historical help string; the modes above are the ones used by the current driver logic.)*

### Output directory selector
The code writes most outputs under `sys%file_path` (default: `data`).
- `-path d` : use `data/` (default)
- `-path .` : write into the current directory (`.`)

> Note: this flag is implemented as a single-character selector (`d` or `.`).

### Quick example (single run)
```bash
mkdir -p data
./mpol   -prep 50 -al 0.1 -del 0.1   -nm 1200 -wc 1 -wmax 3   -k0 0.0783 -x0 -700 -sigma 0.005 -n 0.5 -xmin 400   -dt 0.05 -tmax 1500   -npini 6 -npadd 20 -p0 0.000005 -tref 3 -me 1e-7
```

---

## Ground-state data (`gs_data/`)

Several preparation modes (`prep = 50, 51, 54, 60`) load a dressed ground-state file. By default the code looks in the `gs_data/` directory for a filename of the form:

```
gs_data/FinalFksP_al<alpha>_del<Delta>_pols<np>_nk<nmode>_wm<wmax>_wc<wc>.d
```

The provided `gs.sh` shows how ground-state files were generated on the authors’ cluster (it assumes a separate ground-state executable and cluster paths). If you do not have `gs_data/`, use a preparation mode that does not require it (e.g. `prep = 52`) or generate the expected file via your own ground-state solver.

---

## Outputs

### Key scalar outputs (scattering)
At the end of each trajectory, the code writes a short file in the **current working directory**:
- `file_<k0>.d` — first line contains `k0  R  T` (reflection and transmission).

### Time-series diagnostics (written under `data/` by default)
Filenames are tagged by a parameter string produced by `systm:parameterchar(sys)`.

Common time-trace files include:
- `data/ErEN_<...>.d` — time, (internal) error metric, energy, norm.
- `data/spinXYZ_<...>.d` — time, ⟨σx⟩, ⟨σy⟩, ⟨σz⟩.
- `data/np_<...>.d` — time, number of coherent components in the basis.
- `data/trefs_<...>.d` — records times / thresholds used during adaptive basis growth.

Field / spectrum / correlation outputs commonly include:
- `data/nk_<...>.d`, `data/nk_av_<...>.d`, `data/nphk_<...>.d` — k-space photon observables.
- `data/nx_<...>.d`, `data/nphx_<...>.d`, `data/nbx_<...>.d` — real-space photon observables.
- `data/g1_<...>.d`, `data/g2_<...>.d`, `data/g2_0_<...>.d` — first/second order correlation outputs (when enabled).

Failure / debugging:
- `data/FAILED_<...>.d` — emitted when the evolution hits a numerical failure mode and writes a “fix” snapshot.

---

## Parameter sweeps on a cluster

`bash_scat.sh` and `bash_1scat.sh` demonstrate how the executable was launched over a grid of `k0` values using `qsub`, with each run producing its own folder of outputs. These scripts are cluster-specific:
- they assume a `jobfile` template exists,
- they use absolute scratch paths,
- they create per-run directories and copy `gs_data/` into each run folder.

If you adapt them:
1. Set `MPOL_DIR_SCRATCH` (or equivalent) to your working directory.
2. Ensure `mpol` is built and available at the expected location.
3. Provide a scheduler-compatible `jobfile` template.
4. Confirm `gs_data/` is available for the chosen `prep` mode.

---

## Notes on numerics and performance

- Memory and runtime scale roughly with `nmode × np` (and with additional overhead from the linear solves).
- Results depend on discretization (`nmode`, `wmax`, cutoff), timestep (`dt`), and basis size (`npini`, `npadd`).
- If you see `FAILED_...` outputs, try smaller `dt`, larger `npini/npadd`, or adjust adaptive thresholds (`-me`, `-mde`, `-tref`).

---

## Citation

If you use this code or data derived from it, please cite:

N. Gheeraert *et al.*, “Particle production in ultrastrong-coupling waveguide QED”, **Phys. Rev. A 98**, 043816 (2018).
