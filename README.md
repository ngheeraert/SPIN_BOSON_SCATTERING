# Multi-polaron (MCS) waveguide-QED dynamics (Fortran)

This repository contains a Fortran implementation of **time-dependent variational dynamics** for a two-level system (“qubit”) coupled to a 1D bosonic continuum (spin–boson / waveguide-QED setting).  
It is the simulation engine used for the core results of:

> N. Gheeraert *et al.*, **Particle production in ultrastrong-coupling waveguide QED**, *Phys. Rev. A* **98**, 043816 (2018).

The code evolves a many-body wavefunction written as a **superposition of multimode coherent states** (often called *multi-polaron* or *MCS* ansatz), and outputs plain-text `.d` files for spectra, photon densities, and correlation functions.

---

## Quick start

```bash
# 1) build
make            # gfortran + BLAS/LAPACK (see makefile)
# or:
# make -f makefile_intel

# 2) run (example close to the paper’s Fig. 3 parameter set)
mkdir -p data
./mpol \
  -prep 50 -al 0.1 -del 0.1 -wc 1 -wmax 3 -nm 1200 \
  -k0 0.16 -sigma 0.005 -n 0.5 -x0 -700 -xmin 400 \
  -dt 0.05 -tmax 1500 \
  -npini 6 -npadd 24 -p0 1e-6 -tref 3 -me 1e-7
```

Outputs land in `data/` by default (plus a small `file_<k0>.d` summary in the run directory).

---

## Physics in one page

### Model: spin–boson / waveguide QED
The simulation targets a standard waveguide-QED model (spin–boson Hamiltonian): a two-level system coupled to a continuum of modes of a 1D waveguide. In the paper’s convention one works in a basis where the coupling is diagonal in the qubit degree of freedom (σz-coupling), while the bare splitting appears as a transverse term (σx). The waveguide dispersion is taken linear, ωk = |k| (c = 1), and the coupling uses a smooth high-frequency cutoff.

**What this buys you physically:** counter-rotating terms are fully retained (no RWA), which is essential in the ultrastrong-coupling regime. This is precisely the regime where *off-resonant* nonlinear emission can be dominated by non-RWA “particle-production” pathways (frequency down-conversion into multiple lower-energy photons).

### Ultrastrong coupling & dressing cloud
At α ~ 0.1–0.3 the qubit is **nonperturbatively dressed** by a cloud of virtual waveguide photons, and its resonance is renormalized (Δ → ΔR). Scattering must therefore be performed on top of a dressed ground state, not a bare vacuum.

### Variational many-body wavefunction: multi-polaron / MCS ansatz
The core approximation is a truncated superposition of multimode coherent states:
- the total wavefunction is expanded into **Ncs coherent components** (“polarons”),
- each component carries **mode displacements** (one set for each qubit state) and a complex weight.

The coherent states are genuinely multimode (hundreds to thousands of discretized modes), which is what allows the method to capture a dressed impurity *and* outgoing wavepackets in a single representation.

### Equations of motion (TDVP)
The time evolution is obtained from the **time-dependent variational principle** (TDVP): inserting the MCS ansatz into a real Lagrangian and minimizing with respect to all variational parameters yields coupled ODEs for the weights and displacements. The integrator used in this code is RK4 (see `output.f90`).

### Scattering protocol implemented here (conceptually)
1. **Prepare** the dressed qubit+waveguide ground state (many-body “cloud”).
2. **Add** an incoming wavepacket (typically a coherent Gaussian pulse) far from the impurity.
3. **Evolve** in time until the outgoing packets separate from the dressing region.
4. **Analyze** the outgoing field: elastic peaks at ±k0 plus inelastic emission (frequency conversion and particle production).
5. **Optionally** extract number-resolved spectra and second-order correlations g².

---

## Discretization, cutoff, and coupling conventions (as implemented)

### Frequency grid
The continuum is discretized on positive frequencies:
- `nmode` modes between 0 and `wmax`
- uniform spacing `dk = wmax / nmode`
- mode frequencies `ω_i = dk * (i - 1/2)` for i = 1..nmode

### Coupling function
The code uses an Ohmic-like scaling with a smooth exponential cutoff:
- `g_i ∝ sqrt( α * ω_i * dk * exp(-ω_i/wc) )`

This corresponds to an Ohmic spectral density with exponential cutoff, consistent with the discussion in the paper (Appendix notes also compare to a hard cutoff for some correlation computations).

### Units
Typical runs set the cutoff scale to 1 (`wc = 1`) and express all frequencies in units of ωc (i.e. `del`, `k0`, `sigma`, `wmax` are dimensionless).

---

## Convergence and adaptive basis growth (important)

The MCS basis must be large enough to:
- represent the dressed ground state **and**
- represent the outgoing many-body radiation (including weak inelastic tails).

Two practical mechanisms are used:
1. **Start** with a small number `npini` of coherent components.
2. **Grow** the basis by adding up to `npadd` extra components when an error monitor exceeds a threshold `me` (see paper discussion of monitoring `Err(t) = ||(i∂t - H)|Ψ(t)||²` and adding states with near-zero weights).

Heuristics that usually help:
- Decrease `dt` if energy/norm drift appears.
- Increase total coherent states (`npini + npadd`) when inelastic spectra look noisy (inelastic signals are often orders of magnitude smaller than elastic peaks).
- For strong power saturation, larger Ncs is essential.

---

## Repository contents

- `main.f90` — program entry point: parses CLI args, runs one trajectory; includes a post-processing mode for g².
- `systm.f90` — core types (`param`, `state`, `traj`) and physics routines (mode grid, couplings, overlaps, observables, EOM).
- `output.f90` — time-integration driver (RK4) and file writers (time traces, spectra, correlations).
- `typedefs.f90` — numeric kinds; compile with `-DDP` for double precision.
- `consts.f90` — constants and small helpers.
- `lapackmodule.f90`, `inverse.f90` — LAPACK helpers / linear algebra used in the EOM.
- `makefile`, `makefile_intel` — build recipes (gfortran + BLAS/LAPACK or ifort + MKL).
- `bash_scat.sh`, `bash_1scat.sh`, `gs.sh` — examples for cluster sweeps (qsub-style; may need adaptation).

---

## Building

### Recommended: double precision
Precision is selected in `typedefs.f90` via preprocessing:
- `-DDP` enables double precision (`r_type = dp`)

**Using the Makefile (gfortran):**
```bash
make clean && make
```

**Manual build (gfortran + BLAS/LAPACK):**
```bash
gfortran -O3 -cpp -DDP \
  typedefs.f90 consts.f90 lapackmodule.f90 inverse.f90 systm.f90 output.f90 main.f90 \
  -llapack -lblas -o mpol
```

**Intel example (ifort + MKL):**
```bash
make -f makefile_intel clean && make -f makefile_intel
```

---

## Running: key flags

`getParameters(sys)` parses command-line arguments. Common flags:

### Physical / model parameters
- `-al <alpha>`   : dimensionless coupling strength α
- `-del <Delta>`  : bare qubit splitting Δ (in units of ωc when wc=1)
- `-wc <wc>`      : cutoff frequency scale (appears in exp(-ω/wc))
- `-wmax <wmax>`  : maximum retained frequency (band edge of discretization)
- `-nm <nmode>`   : number of discretized modes

### Wavepacket parameters
- `-k0 <k0>`      : central momentum/frequency of incoming packet
- `-k0_2 <k0_2>`  : secondary k0 (used in two-tone preparation)
- `-x0 <x0>`      : initial packet center in real space (choose negative to start left of impurity)
- `-sigma <sigma>`: spectral width (smaller = more monochromatic)
- `-n <n_wp>`     : mean photon number / intensity parameter of coherent input (n̄)

### Time stepping
- `-dt <dt>`      : time step
- `-tmax <tmax>`  : final time (choose long enough that packets exit the interaction region)

### Variational basis (coherent components)
- `-npini <np>`   : initial number of coherent components (Ncs)
- `-npadd <n>`    : number of components to add during evolution
- `-p0 <p0>`      : small initial weight assigned to new components
- `-tref <tref>`  : reference time used in the addition logic
- `-me <merr>`    : error threshold used to trigger basis growth
- `-mde <max_deltaE>` : energy-drift threshold for stability checks
- `-bt <0/1>`     : enable/disable “back-tracking” (adaptive stability)

### Preparation mode (`-prep`)
The prep mode selects how the initial state is built (inside the trajectory driver). Commonly used modes in this repo:
- `50` : coherent-state packet + dressed ground-state background (**requires `gs_data/`**)
- `51` : single-photon-like packet + dressed ground-state background (**requires `gs_data/`**)
- `52` : single-photon-like packet without ground-state subtraction
- `54` : two-frequency coherent packet + dressed ground-state background (**requires `gs_data/`**)
- `60` : restart mode (loads a saved state after setting up metadata)
- `2601..2699` : post-processing mode: reload a saved state and compute g² (see `main.f90`)

### Output directory selector
- `-path d` : write to `data/` (default)
- `-path .` : write to current directory

---

## Ground-state data (`gs_data/`)

Several preparation modes load a dressed ground-state file. By default the code expects:

```
gs_data/FinalFksP_al<alpha>_del<Delta>_pols<np>_nk<nmode>_wm<wmax>_wc<wc>.d
```

If you do not have `gs_data/`, use a prep mode that does not require it (e.g. `prep = 52`) or generate the expected file via your own ground-state solver.

---

## Outputs

### Key scalar outputs (scattering)
At the end of each trajectory:
- `file_<k0>.d` — first line contains `k0  R  T` (reflection and transmission)

### Time-series diagnostics (under `data/` by default)
Filenames are tagged by a parameter string.

Common files:
- `data/ErEN_<...>.d` — time, error metric, energy, norm
- `data/spinXYZ_<...>.d` — time, ⟨σx⟩, ⟨σy⟩, ⟨σz⟩
- `data/np_<...>.d` — time, number of coherent components in the basis
- `data/trefs_<...>.d` — times/thresholds used during adaptive basis growth

Field / spectrum / correlation outputs:
- `data/nk_<...>.d`, `data/nk_av_<...>.d`, `data/nphk_<...>.d` — k-space photon observables / spectra
- `data/nx_<...>.d`, `data/nphx_<...>.d`, `data/nbx_<...>.d` — real-space photon observables / densities
- `data/g1_<...>.d`, `data/g2_<...>.d`, `data/g2_0_<...>.d` — correlation outputs (when enabled)

Failure / debugging:
- `data/FAILED_<...>.d` — emitted when the evolution hits a numerical failure mode and writes a “fix” snapshot

---

## Reproducing the paper’s headline results (guidelines)

The paper emphasizes that inelastic features can be **tiny** compared to elastic peaks, and require careful convergence (enough modes, enough coherent states, stable timestep).

### (A) Reflection curve (paper Fig. 2 style)
- Sweep `k0` across the (renormalized) resonance.
- Use small `sigma` (narrowband pulse) and a long enough `tmax` that packets exit the interaction region.
- For stronger input power (larger `n`), increase total coherent states (paper discussion indicates Ncs ~ 16 for convergence of strong-power reflection).

Practical recipe:
1) Fix α, Δ, wc, wmax, nmode, sigma, n, etc.  
2) Run a bash loop over `k0`, collecting `file_<k0>.d` into a single curve.

### (B) Off-resonant frequency conversion / particle production (paper Fig. 3 style)
A representative parameter set discussed in the paper is:
- α = 0.1, Δ = 0.1 ωc
- k0 = 0.16 ωc, σ = 0.005 ωc, n̄ = 0.5
- Nmodes ≈ 1200 and Ncs ≈ 30
- integrate until the pulse is far away from the impurity region

The example command in **Quick start** is designed to be close to this regime.

### (C) g² correlations
Second-order correlations are typically harder to converge than average photon numbers. The paper notes that for some g² computations it is useful to simplify the cutoff (hard cutoff) to reach a larger number of coherent states.

---

## Notes on performance
Rough scaling:
- memory/time grows with `nmode × Ncs**6`, with additional overhead from dense linear algebra.

If you get failures:
- reduce `dt`
- increase `npini` / `npadd`
- tighten or retune `-me`, `-mde`, `-tref`

---

## Citation
If you use this code or data derived from it, please cite:

N. Gheeraert *et al.*, “Particle production in ultrastrong-coupling waveguide QED”, **Phys. Rev. A 98**, 043816 (2018).
