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

## Four entry points

Everything is driven from four scripts in the repo root. Helper code lives in
`lib/` — the original per-task script bodies, moved there almost verbatim.

| Script | Purpose |
|--------|---------|
| **`setup_databases.m`** | Builds `inertia_tensors.mat` and `ref_geometry.mat`; can also rename raw ADBSat output into the expected layout. `setup_databases` with no argument does the first two. |
| **`run_sweep.m`** | The ROADM + SAM parameter sweep. Canonical settings from the former `MASTERCODE.m` (June epoch, 98° inclination, `ode45`, one orbit). Writes `drag_sail_attitude_results.mat`. |
| **`run_analyses.m`** | Dispatcher for every secondary analysis. `run_analyses` with no argument lists them; `run_analyses montecarlo_multiparam` runs one. |
| **`make_figures.m`** | Dispatcher for every figure-generation script. `make_figures` with no argument lists them; `make_figures beta_cd_suite` runs one. |

`run_analyses` options: `montecarlo_multiparam`, `montecarlo_sensitivity`,
`montecarlo_barchart`, `accom_pdf`, `accom_transient`, `sensitivity_theta0`,
`sensitivity_isolated`, `tumbling_multiday`, `tumbling_sphere`,
`error_convergence`, `decay_map`, `lifetime_kinghele`.

`make_figures` options: `beta_cd_suite`, `sam_sweep`, `sam_extra`, `adm`, `pdf`,
`sam_vs_uniform`, `cd_feedback`, `angular_velocity`, `density_vs_time`,
`geometry_patch`.

### Typical run order

```matlab
setup_databases                       % (only needed to regenerate the .mat inputs)
run_sweep                             % -> drag_sail_attitude_results.mat
make_figures beta_cd_suite            % headline Cd_bar / beta figures
run_analyses montecarlo_multiparam    % -> mc_workspace.mat, mc_results.mat
run_analyses lifetime_kinghele        % orbital-lifetime figures
```

Some analyses depend on earlier outputs: `montecarlo_barchart`,
`sensitivity_theta0` and `sensitivity_isolated` need `mc_workspace.mat` from
`montecarlo_multiparam`; `decay_map` needs `decay_times.mat`.
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
- **Parallel Computing Toolbox** — only for `run_analyses sensitivity_isolated`
- ADBSat output (see above) for `run_sweep` and most analyses

## Layout

```
setup_databases.m   run_sweep.m   run_analyses.m   make_figures.m   <- entry points

obj files/               Sail geometries (apex half-angle 45°–90°) — ADM input
ref_geometry.mat         Reference geometry per sail (regenerate: setup_databases geometry)
inertia_tensors.mat      Pitch inertia per sail (regenerate: setup_databases inertia)
decay_times.mat          Precomputed decay lifetimes for run_analyses decay_map

lib/
  project_root.m         Resolves the repo root (replaces the old absolute paths)
  get_ref_geometry.m     Reads a .obj -> reference area / length
  analyses/*.m           One self-contained script per run_analyses option
  figures/*.m            One script per make_figures option
  variants/*.m           Alternative sweep formulations (not wired to run_sweep)
```

Generated output (`*.png`, `*.fig`) and large result `*.mat` files are
git-ignored — the repo is code plus the small hand-built inputs only.
Regenerate figures with `make_figures`.

## Notes

- The `lib/` scripts are the original bodies with only the hard-coded absolute
  paths replaced by `project_root()`. Each analysis keeps
  its own copy of the attitude ODE / density-torque / SAM code, because the
  differences between them (density scaling, `wrapTo180` binning, `ode45` vs
  `ode113`, equatorial vs Sun-synchronous) are intentional.
- Alternative sweep formulations are in `lib/variants/`: `NEWMASTER.m` (January
  epoch, equatorial, `ode113`, 15 orbits) and `mastergit.m`
  (`run_drag_sail_sweep(config)` function form). Run these directly; they are
  not called by `run_sweep`. `mastergit.m` expects the working directory to be
  the repo root.
- `.fig` / `.png` figures and large `.mat` results are git-ignored; regenerate
  them locally.

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
