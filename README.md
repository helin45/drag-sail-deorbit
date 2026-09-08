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
| **ROADM** — Reduced-Order Attitude Dynamics Model | Reads $C_M(\theta)$ from the ADM and propagates the 1-DOF pitch equation of motion over 15 orbits. Aerodynamic torque at each timestep interpolates $C_M(\theta)$ at the instantaneous pitch angle and scales by orbit-varying dynamic pressure from NRLMSISE-00. Outputs $\theta(t)$, $\dot\theta(t)$. | Yes |
| **SAM** — Statistical Averaging Model | Builds a residence-time PDF $p(\theta)$ from the ROADM time history and uses it as weights in a Riemann sum over $C_D(\theta)$, giving $\bar{C}_D$ and $\bar{\beta}$ over the full $(\phi,\alpha,h)$ space. | Yes |

## The ADM / ADBSat data is not included

The aerodynamic database is produced with **ADBSat** (Sinpetru et al., 2022),
which belongs to its authors and is not redistributed here. The scripts expect
ADBSat output at:

```
adbsat_processed/<alt>km/<phi>deg_CLL_accom_<a>.mat   e.g. adbsat_processed/450km/75deg_CLL_accom_0p75.mat
```

To reproduce: obtain ADBSat, run it on the `*.obj` sail geometries in this repo
for each altitude / accommodation-coefficient combination, and place the results
in that folder layout. The per-panel `.mat` files must contain `aedb.aero` with
fields `Cm_BY` and `Cf_wX`.

## Requirements

- MATLAB (developed on R2023+)
- **Aerospace Toolbox** — `atmosnrlmsise00`, `eci2lla`
- ADBSat output (see above) for anything that runs the sweep

## Layout

```
*.obj                         Sail geometries (apex half-angle 45°–90°) — ADM input
drag_sail_inertias.m          Builds inertia_tensors.mat
get_ref_geometry.m            Reads .obj -> reference area / length; builds ref_geometry.mat
ref_geometry.mat              Reference geometry per sail (tracked)
inertia_tensors.mat           Pitch inertia per sail (tracked)
decay_times.mat               Precomputed decay lifetimes for the decay plots (tracked)

MASTERCODE.m                  Full (phi, accom, alt) ROADM+SAM sweep  -> drag_sail_attitude_results.mat
NEWMASTER.m                   Sweep with Roberts & Harkness (2007) epoch/IC -> attitude_SAM_sweep_*.mat
mastergit.m                   Same sweep refactored as run_drag_sail_sweep(config)
accomMASTECODE.m              Accommodation-coefficient focused sweep
accomtransient.m              Accommodation effect on the pitch transient

montecarli.m                  Monte Carlo over phi0, dphi0, rho scaling
MonteCarlo/mc_sensitivity_new.m   Monte Carlo sensitivity (theta0, dphi0, rho_scale)
barCHARTMC.m                  Monte Carlo bar-chart summary
error_analysis.m             n_steps convergence, bin-width, density-scaling checks

ORBITALDECAY.m               Decay-time colour maps over (accom, phi)
Publishing/KingHeleOrbital.m King–Hele orbital lifetime vs altitude

plotting.m, attitude_SAM_plots_all_v2.m, newfigs.m, admplot.m, pdfplots.m,
CDCOMPARE.m, fefdbackplot.m, newplothelin.m, angularvel.m   Figure generation
patchold.m, patchsingle.m, patchsingleem.m                  Geometry rendering
spheretest.m, test_tumbling.m, test_tehta0.m, isolated_tests.m   Validation tests
```

## Known issues

- Several scripts contain hard-coded absolute paths
  (`/Users/helintaha/Library/CloudStorage/...`). Edit `adbsat_base` /
  `addpath` lines near the top of each script to point at your local
  `adbsat_processed/`.
- `.fig` files and large `.mat` result files are git-ignored; PNG exports of the
  figures are tracked.

## Reference

Sinpetru, L. A., et al. *ADBSat: methodology of a novel panel method tool for
aerodynamic analysis of satellites.* Computer Physics Communications, 2022.

## License

MIT — see [LICENSE](LICENSE).
