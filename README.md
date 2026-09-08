# Drag-Sail Attitude Dynamics & Orientation-Averaged Ballistic Coefficient

MATLAB code for a dissertation study on passive drag-sail de-orbit performance. A
single degree-of-freedom pitch model is propagated under aerodynamic and
gravity-gradient torque, and the resulting residence-time distribution is used to
weight the orientation-dependent drag coefficient into an **orientation-averaged
drag coefficient** $\bar{C}_D$ and **ballistic coefficient** $\bar{\beta}$ across a
sweep of sail geometry, surface accommodation coefficient, and orbital altitude.

## Pipeline

| Stage | What it does | In this repo? |
|-------|--------------|:---:|
| **ADM** — Aerodynamic Database Model | Generates $C_D(\theta)$ and $C_M(\theta)$ over the full pitch range for each geometry / accommodation coefficient / altitude, using the CLL gas–surface interaction model in **ADBSat**. | **No** — see below |
| **ROADM** — Reduced-Order Attitude Dynamics Model | Reads $C_M(\theta)$ from the ADM and propagates the 1-DOF pitch equation of motion. Aerodynamic torque at each timestep interpolates $C_M(\theta)$ at the instantaneous pitch angle and scales by orbit-varying dynamic pressure from NRLMSISE-00. Outputs $\theta(t)$, $\dot\theta(t)$. | Yes |
| **SAM** — Statistical Averaging Model | Builds a residence-time PDF $p(\theta)$ from the ROADM time history and uses it as weights in a Riemann sum over $C_D(\theta)$, giving $\bar{C}_D$ and $\bar{\beta}$ over the full $(\phi,\alpha,h)$ space. | Yes |

## Scripts

Five scripts in the repo root, run in this order:

| Script | Purpose | Output |
|--------|---------|--------|
| **`setup_databases.m`** | Builds `inertia_tensors.mat` and `ref_geometry.mat`; `setup_databases rename '<path>'` also reshapes raw ADBSat output. No argument = build both. | `*.mat` inputs |
| **`ROADM.m`** | Propagates the 1-DOF pitch equation of motion for every $(\phi,\alpha,h)$ case. Aerodynamic + gravity-gradient torque, NRLMSISE-00 dynamic pressure, `ode45`. | `roadm_results.mat` — $\theta(t)$, $\dot\theta(t)$ per case |
| **`SAM.m`** | Loads `roadm_results.mat`, builds the residence-time PDF $p(\theta)$, Riemann-sums it against $C_D(\theta)$ for $\bar C_D$, and forms $\bar\beta = m/(\bar C_D A_{\mathrm{ref}})$. | `sam_results.mat` — `Cd_bar_3D`, `beta_bar_3D` |
| **`montecarlo.m`** | Monte Carlo sensitivity of $\bar C_D$ to $\theta_0$, $\dot\theta_0$ and a density scale factor. All settings in the `CFG` block at the top. | `mc_results.mat`, `mc_workspace.mat` |
| **`make_figures.m`** | Dispatcher for the plotting scripts in `lib/figures/`. `make_figures` with no argument lists the sets; `make_figures beta_cd_suite` runs one. | `.fig` / `.png` |

`ROADM.m` runs **one** orbit by default (the canonical `MASTERCODE.m` setting).
Set `N_ORBITS = 15` near the top for the dissertation headline runs.

### Typical run order

```matlab
setup_databases          % only to regenerate the .mat inputs
ROADM                    % -> roadm_results.mat   (slow; the full sweep)
SAM                      % -> sam_results.mat
make_figures beta_cd_suite
montecarlo               % -> mc_results.mat, mc_workspace.mat
```

### `validation/`

Supporting one-off studies and checks, kept as standalone scripts (run each
directly, not wired to a dispatcher):

`accom_pdf`, `accom_transient` (accommodation-coefficient effects) ·
`tumbling_multiday`, `tumbling_sphere` (multi-day propagation / Roberts-sphere
validation) · `error_convergence` (step-count / bin-width / density-scaling
convergence) · `sensitivity_theta0`, `sensitivity_isolated` (isolated
initial-condition sweeps) · `decay_map`, `lifetime_kinghele` (orbital
decay-time and King-Hele lifetime figures) · `montecarlo_multiparam`,
`montecarlo_sensitivity`, `montecarlo_barchart` (the granular Monte Carlo
scripts `montecarlo.m` was distilled from) · `patch_one_case` (recompute and
splice one bad cell in `sam_results.mat`).

`montecarlo_barchart`, `sensitivity_theta0`, `sensitivity_isolated` load
`mc_workspace.mat` — run `montecarlo` first. `decay_map` needs `decay_times.mat`.
`sensitivity_isolated` uses `parfor` (Parallel Computing Toolbox).

## The ADM / ADBSat data is not included

The aerodynamic database is produced with **ADBSat** (Sinpetru et al., 2022),
which belongs to its authors and is not redistributed here. The scripts expect
ADBSat output at:

```
adbsat_processed/<alt>km/<phi>deg_CLL_accom_<a>.mat   e.g. adbsat_processed/450km/75deg_CLL_accom_0p75.mat
```

To reproduce: obtain ADBSat, run it on the `*.obj` sail geometries in this repo
for each altitude / accommodation-coefficient combination, then
`setup_databases rename '<path to ADBSat inou/results>'`. Each per-panel `.mat`
must contain `aedb.aero` with fields `Cm_BY` and `Cf_wX`.

## Requirements

- MATLAB (developed on R2023+)
- **Aerospace Toolbox** — `atmosnrlmsise00`, `eci2lla`
- **Parallel Computing Toolbox** — only for `validation/sensitivity_isolated.m`
- ADBSat output (see above) for `ROADM.m`, `montecarlo.m` and most of `validation/`

## Layout

```
setup_databases.m   ROADM.m   SAM.m   montecarlo.m   make_figures.m   <- run these

obj files/               Sail geometries (apex half-angle 45°–90°) — ADM input
ref_geometry.mat         Reference geometry per sail (regenerate: setup_databases geometry)
inertia_tensors.mat      Pitch inertia per sail (regenerate: setup_databases inertia)
decay_times.mat          Precomputed decay lifetimes for validation/decay_map.m

lib/
  project_root.m         Resolves the repo root (replaces the old absolute paths)
  get_ref_geometry.m     Reads a .obj -> reference area / length
  figures/*.m            One script per make_figures option
  variants/*.m           Alternative full-sweep formulations (NEWMASTER.m, mastergit.m)
validation/*.m           Supporting studies and checks (see above)
```

Generated output (`*.png`, `*.fig`) and large result `*.mat` (`roadm_results.mat`,
`sam_results.mat`, …) are git-ignored — the repo is code plus the small
hand-built inputs only.

## Notes

- `lib/figures/` and `validation/` scripts are the original per-task bodies,
  moved almost verbatim (only the hard-coded absolute paths were replaced by
  `project_root()`). They each keep their own copy of the attitude ODE /
  density-torque code where it differs intentionally (density scaling,
  `wrapTo180` binning, `ode45` vs `ode113`, equatorial vs Sun-synchronous).
- `lib/variants/`: `NEWMASTER.m` (January epoch, equatorial, `ode113`, 15 orbits)
  and `mastergit.m` (`run_drag_sail_sweep(config)` function form) — alternative
  formulations of the `ROADM.m` + `SAM.m` sweep. Run directly; `mastergit.m`
  expects the working directory to be the repo root.
- `make_figures beta_cd_suite` is the headline $\bar C_D$ / $\bar\beta$ figure
  set. `sam_sweep` / `sam_extra` need the 15-orbit output files and will error
  without them.

## References

- **ADBSat** — <https://github.com/nhcrisp/ADBSat>
- Sinpetru, L. A., Crisp, N. H., et al. *ADBSat: methodology of a novel panel
  method tool for aerodynamic analysis of satellites.* Computer Physics
  Communications 275 (2022) 108326.
- Sinpetru, L. A., Crisp, N. H., et al. *ADBSat: verification and validation of a
  novel panel method for quick aerodynamic analysis of satellites.* Computer
  Physics Communications 272 (2022) 108234.

## License

MIT — see [LICENSE](LICENSE).
