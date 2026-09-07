# PIMS — Polarimetric Imaging-based Mirror Soiling detection

Simulation code for the DoLP-vs-reflectance model in:

> M. Tian, M. Z. E. Rafique, K. Chidambaranathan, R. Brost, D. Small, D. Novick,
> J. Yellowhair, Y. Yao, *A Polarimetry-based Field-deployable Non-interruptive
> Mirror Soiling Detection Method*, Solar Energy (2026).

Section and equation numbers in the code comments refer to that paper.


## Start here

1. Open MATLAB in this folder.
2. Run `run\_pims\_simulation.m` to execute the cleaned version of the Sahara-soil simulation workflow.
3. Use `functions/process\_polarization\_image.m` for raw polarization-camera image processing.

## Directory structure

* `run\_pims\_simulation.m` - main student-facing simulation entry point.
* `functions/` - standalone Mie, Rayleigh-skylight, geometry, Mueller-matrix, and utility functions extracted from the original script.
* `data/` - uploaded particle/skylight data.
* `examples/` - short usage examples.

## Scientific model

The simulation combines two contributions in Stokes space:

1. **Reflected skylight** from the mirror, with polarization determined by the Rayleigh-scattered skylight pattern.
2. **Sunlight scattered by soil particles**, modeled with single-particle Mie scattering and weighted by the particle-size distribution.

Particle concentration is varied to change the modeled soiled area fraction. The combined Stokes vector is converted to DoLP, producing the model relationship between DoLP and relative reflectance.

## Important reproducibility note

This package is a conservative refactor. Numerical constants, empirical scale factors, and measurement arrays in the selected source script are intentionally retained. They should not be interpreted as universal physical constants. See `docs/CODE\_AUDIT.md`.

## MATLAB dependencies

The uploaded code is largely self-contained. It uses standard MATLAB numerical, plotting, image, and interpolation functions. `rms` may require Signal Processing Toolbox in some MATLAB versions. If unavailable, replace `rms(error)` with `sqrt(mean(error.^2))`.



