# Simulation Workflow

1. Load the particle-size distribution.
2. Define optical parameters, source direction, and camera angular coverage.
3. Calculate the Rayleigh skylight Stokes map.
4. Apply the skylight intensity weighting.
5. Calculate Mie-scattered sunlight maps for each particle size.
6. Sweep particle concentration to represent different surface-soiling levels.
7. Convert particle concentration to surface coverage and relative reflectance.
8. Weight particle-size-resolved scattering using Mie attenuation.
9. Combine scattered sunlight and reflected skylight in Stokes space.
10. Calculate camera-averaged DoLP.
11. Scale the simulated DoLP curve to the measured maximum DoLP.
12. Interpolate the model at measured reflectance values and report RMS DoLP error.
