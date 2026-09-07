function results = run_pims_spectral_400_700()
%% PIMS broadband simulation: 400-700 nm
% Broadband extension of the PIMS DoLP-to-reflectance model.
%
% Place this file in the repository root next to:
%   functions/
%   data/particle_size_distribution.txt
%
% The calculation follows the same model structure as the monochromatic
% simulation, but repeats the Mie calculation over a visible wavelength
% grid and combines the wavelength-resolved Stokes vectors before DoLP is
% evaluated.
%
% IMPORTANT:
%   Do not average DoLP over wavelength. The detector integrates optical
%   signal first, so S0, S1, S2 and S3 are spectrally combined first, and
%   DoLP is calculated from the broadband Stokes vector afterward.
%
% This script uses the existing helper functions in functions/ and does not
% modify them.

close all;
clc;
t_start = tic;

project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'functions'));

%% ------------------------------------------------------------------------
% 1. Broadband model settings
% -------------------------------------------------------------------------

% Wavelength grid. A 20 nm step is a practical starting point with the
% current Mie implementation. For a final convergence test, repeat with
% 10 nm or smaller steps and verify that the output curve changes little.
wavelength_nm = 400:20:700;
wavelength_m = wavelength_nm * 1e-9;

medium_refractive_index = 1.0;

% Current repository model uses a wavelength-independent particle index.
% Replace this vector with measured n(lambda)+i*k(lambda), if available.
particle_refractive_index = (1.57 + 0i) * ones(size(wavelength_nm));

% Broadband sunlight/skylight ratio retained from the fitted model.
% The spectral weights within each component are normalized separately, so
% these remain the integrated broadband fractions.
sun_fraction = 0.82;
sky_fraction = 1 - sun_fraction;

% Measured maximum DoLP used by the existing model as a final scale factor.
% Set to 1 to inspect the unscaled broadband model shape.
max_measured_dolp = 0.6428;

% Fixed source direction used by the current repository example.
sun_azimuth_deg = 205.3034;
sun_zenith_deg = 69.6823;

% Camera angular coverage.
camera_azimuth_center_deg = 180;
camera_zenith_center_deg = 45;
camera_azimuth_span_deg = 360;
camera_zenith_span_deg = 90;

camera_azimuth_limits = camera_azimuth_center_deg + ...
    [-0.5, 0.5] * camera_azimuth_span_deg;
camera_zenith_limits = camera_zenith_center_deg + ...
    [-0.5, 0.5] * camera_zenith_span_deg;

% Angular grid. The broadband run is much more expensive than the
% monochromatic run because the Mie map is recalculated at every wavelength.
angle_resolution_deg = 1;
mirror_azimuth_deg = 0:angle_resolution_deg:360;
mirror_zenith_deg = 0:angle_resolution_deg:90;

% Fraction retained from the existing model for scattered light entering
% the near-specular camera acceptance region.
camera_scatter_acceptance = 0.06;

% Retain attenuation/cross-section weighting of particle sizes.
use_cross_section_weighting = true;

% Area used to convert the particle count distribution to surface coverage.
simulation_area_m2 = (280e-6)^2;

%% ------------------------------------------------------------------------
% 2. Spectral weighting
% -------------------------------------------------------------------------
% The manuscript specifies a 400-700 nm operating range but does not define
% a unique wavelength binning/weighting rule. Here we use a physically
% transparent broadband weighting that can be replaced with measured data:
%
%   sunlight weight ~ solar spectrum x detector response x mirror response
%   skylight weight ~ sunlight weight x Rayleigh(lambda^-4)
%
% The two weight vectors are normalized independently so sun_fraction and
% sky_fraction continue to represent the fitted broadband mixture.

solar_temperature_K = 5782;
solar_spectrum = planck_spectral_shape(wavelength_m, solar_temperature_K);

% Relative detector responsivity. Keep at 1 for a flat detector model.
% If you have camera QE(lambda), a photon-counting response is proportional
% to QE(lambda)*lambda (the constant 1/(h*c) cancels after normalization).
detector_response = ones(size(wavelength_nm));

% Relative mirror throughput. Replace with measured spectral reflectance if
% desired. A scalar-like flat response reproduces the current simplification.
mirror_response = ones(size(wavelength_nm));

sun_weight = solar_spectrum .* detector_response .* mirror_response;

% Relative spectral shape of Rayleigh-scattered skylight.
rayleigh_weight = (550 ./ wavelength_nm).^4;
sky_weight = sun_weight .* rayleigh_weight;

sun_weight = sun_weight / sum(sun_weight);
sky_weight = sky_weight / sum(sky_weight);

%% ------------------------------------------------------------------------
% 3. Load particle-size distribution
% -------------------------------------------------------------------------
particle_table = readmatrix( ...
    fullfile(project_root, 'data', 'particle_size_distribution.txt'), ...
    'NumHeaderLines', 1);

% Preserve the current repository convention: column 1 is the size passed
% to the Mie calculation and is treated as particle diameter in the model.
% Verify the radius/diameter convention before substituting another PSD.
particle_distribution = [particle_table(:,1), particle_table(:,3) * 250];
particle_distribution = sortrows(particle_distribution, 1);

particle_diameter_nm = particle_distribution(:,1)';
base_particle_count = particle_distribution(:,2)';
particle_radius_m = particle_diameter_nm * 1e-9 / 2;
particle_geometric_area_m2 = pi * particle_radius_m.^2;

%% ------------------------------------------------------------------------
% 4. Rayleigh skylight angular polarization map
% -------------------------------------------------------------------------
% The existing helper is called once to obtain the angular Rayleigh Stokes
% pattern. In the ideal single-scattering Rayleigh model, the normalized
% polarization state versus scattering angle is wavelength independent.
% Wavelength dependence of skylight intensity is therefore handled by the
% sky_weight vector above.

[~,~,~,~,~,~,skylight_zenith_deg,skylight_azimuth_deg,~,~,~,~, ...
    skylight_I,skylight_Q,skylight_U,skylight_V] = ...
    get_skylight_Rayleigh_Scat( ...
        sun_zenith_deg, sun_azimuth_deg, ...
        mirror_azimuth_deg, mirror_zenith_deg);

skylight_data = [skylight_azimuth_deg, skylight_zenith_deg, ...
    skylight_I, skylight_Q, skylight_U, skylight_V];
skylight_data = sortrows(skylight_data, 1);

skylight_camera_data = select_camera_region( ...
    skylight_data, camera_azimuth_limits, camera_zenith_limits);

camera_angles = skylight_camera_data(:,1:2);
normalized_skylight_stokes = normalize_stokes_rows( ...
    skylight_camera_data(:,3:6));

% Since the normalized Rayleigh Stokes pattern is the same at each
% wavelength in this model, spectral integration returns the same state.
% This explicit loop keeps the broadband bookkeeping clear and makes it
% easy to add wavelength-dependent mirror Mueller matrices later.
broadband_skylight_stokes = zeros(size(normalized_skylight_stokes));
for wavelength_index = 1:numel(wavelength_nm)
    broadband_skylight_stokes = broadband_skylight_stokes + ...
        sky_weight(wavelength_index) * normalized_skylight_stokes;
end

%% ------------------------------------------------------------------------
% 5. Wavelength-resolved Mie-scattered sunlight
% -------------------------------------------------------------------------
% For every wavelength:
%   a) update k0 = 2*pi/lambda and the Mie size parameter x,
%   b) calculate the scattering map for every particle size,
%   c) combine particle sizes using the same attenuation-weighting rule as
%      the current monochromatic model,
%   d) normalize the scattered Stokes vector at each camera direction,
%   e) add it to the broadband Stokes vector with the spectral weight.

incident_sun_stokes = [1; 0; 0; 0];
broadband_scattered_stokes_raw = zeros(size(normalized_skylight_stokes));

mean_scattered_dolp_by_wavelength = zeros(size(wavelength_nm));

for wavelength_index = 1:numel(wavelength_nm)
    lambda_m = wavelength_m(wavelength_index);
    k0 = 2*pi/lambda_m;
    particle_index = particle_refractive_index(wavelength_index);

    size_parameter = k0 * medium_refractive_index * particle_radius_m;

    fprintf('Wavelength %d/%d: %.0f nm\n', ...
        wavelength_index, numel(wavelength_nm), wavelength_nm(wavelength_index));

    [I_forward,Q_forward,U_forward,V_forward] = ...
        get_mirror_stokes_map_size_distribution( ...
            mirror_azimuth_deg, mirror_zenith_deg, ...
            sun_azimuth_deg, sun_zenith_deg, ...
            particle_index, medium_refractive_index, ...
            size_parameter, k0, incident_sun_stokes);

    [I_mirrored,Q_mirrored,U_mirrored,V_mirrored] = ...
        get_mirror_stokes_map_size_distribution( ...
            mirror_azimuth_deg, mirror_zenith_deg, ...
            sun_azimuth_deg, 180-sun_zenith_deg, ...
            particle_index, medium_refractive_index, ...
            size_parameter, k0, incident_sun_stokes);

    % Particle-size probability at this wavelength.
    % A uniform concentration multiplier cancels after normalization, so
    % this only needs to be calculated once per wavelength rather than once
    % for every soiling level.
    relative_index = particle_index / medium_refractive_index;
    attenuation_coefficient = zeros(size(particle_diameter_nm));

    for particle_index_id = 1:numel(particle_diameter_nm)
        x = size_parameter(particle_index_id);
        mie_result = mie(relative_index, x);
        q_extinction = mie_result(5) + mie_result(6); % Qsca + Qabs

        surface_number_density = ...
            base_particle_count(particle_index_id) / simulation_area_m2;

        attenuation_coefficient(particle_index_id) = ...
            surface_number_density * ...
            particle_geometric_area_m2(particle_index_id) * q_extinction;
    end

    if use_cross_section_weighting
        particle_probability = ...
            attenuation_coefficient .* base_particle_count;
    else
        particle_probability = base_particle_count;
    end
    particle_probability = particle_probability / sum(particle_probability);

    scattered_I = 0;
    scattered_Q = 0;
    scattered_U = 0;
    scattered_V = 0;

    for particle_index_id = 1:numel(particle_probability)
        w_particle = particle_probability(particle_index_id);

        scattered_I = scattered_I + ...
            (I_forward{particle_index_id,1} + ...
             I_mirrored{particle_index_id,1}) * w_particle;
        scattered_Q = scattered_Q + ...
            (Q_forward{particle_index_id,1} + ...
             Q_mirrored{particle_index_id,1}) * w_particle;
        scattered_U = scattered_U + ...
            (U_forward{particle_index_id,1} + ...
             U_mirrored{particle_index_id,1}) * w_particle;
        scattered_V = scattered_V + ...
            (V_forward{particle_index_id,1} + ...
             V_mirrored{particle_index_id,1}) * w_particle;
    end

    scattered_data = [skylight_data(:,1:2), ...
        real(scattered_I(:)), real(scattered_Q(:)), ...
        real(scattered_U(:)), real(scattered_V(:))];

    scattered_camera_data = select_camera_region( ...
        scattered_data, camera_azimuth_limits, camera_zenith_limits);

    scattered_camera_stokes = real(scattered_camera_data(:,3:6));

    % S1/S2 from Mie theory are dimensionless amplitude functions. The
    % differential scattering cross section is proportional to |S|^2/k^2,
    % so 1/k0^2 is retained here when comparing different wavelengths.
    % This factor is irrelevant in a single-wavelength normalized model but
    % becomes important for broadband spectral integration.
    scattered_camera_stokes = scattered_camera_stokes / k0^2;

    broadband_scattered_stokes_raw = ...
        broadband_scattered_stokes_raw + ...
        sun_weight(wavelength_index) * scattered_camera_stokes;

    normalized_scattered_stokes = normalize_stokes_rows( ...
        scattered_camera_stokes);
    spectral_dolp = sqrt(normalized_scattered_stokes(:,2).^2 + ...
        normalized_scattered_stokes(:,3).^2) ./ ...
        normalized_scattered_stokes(:,1);
    mean_scattered_dolp_by_wavelength(wavelength_index) = ...
        mean(spectral_dolp, 'omitnan');
end

% Normalize only after all wavelengths have been combined. This is the
% key broadband step: the detector integrates Stokes-resolved signal across
% wavelength before DoLP is evaluated.
broadband_scattered_stokes = normalize_stokes_rows( ...
    broadband_scattered_stokes_raw);

%% ------------------------------------------------------------------------
% 6. Sweep soiling level / particle concentration
% -------------------------------------------------------------------------
concentration_modifier = ...
    [0.0225 0.025 0.0272 0.03 0.04 0.045 0.05 0.06 0.07 0.08 ...
     0.09 0.10 0.12 0.15 0.18 0.20 0.22 0.25 0.28 0.30 ...
     0.40 0.50 0.65 1.00 1.50 2.00 3.00 8.00] / 10;

reflectance = zeros(size(concentration_modifier));
average_dolp = zeros(size(concentration_modifier));

for concentration_index = 1:numel(concentration_modifier)
    number_concentration = base_particle_count ./ ...
        concentration_modifier(concentration_index);

    covered_area_m2 = sum( ...
        particle_geometric_area_m2 .* number_concentration);
    covered_fraction = covered_area_m2 / simulation_area_m2;
    uncovered_fraction = 1 - covered_fraction;

    % Preserve the current relative-reflectance approximation.
    reflectance(concentration_index) = uncovered_fraction;

    % Same broadband sunlight/skylight mixing structure used by the
    % monochromatic model, now using spectrally integrated Stokes states.
    total_stokes = ...
        broadband_scattered_stokes * sun_fraction * covered_fraction + ...
        broadband_skylight_stokes * sky_fraction * ...
        (uncovered_fraction + ...
         camera_scatter_acceptance * covered_fraction);

    dolp = sqrt(total_stokes(:,2).^2 + total_stokes(:,3).^2) ./ ...
        total_stokes(:,1);

    average_dolp(concentration_index) = mean(dolp, 'omitnan');
end

%% ------------------------------------------------------------------------
% 7. Scale and plot broadband DoLP-reflectance curve
% -------------------------------------------------------------------------
model_dolp_unscaled = real(average_dolp);
model_dolp = model_dolp_unscaled / max(model_dolp_unscaled) * ...
    max_measured_dolp;

figure('Name','PIMS Broadband 400-700 nm');
plot(reflectance*100, model_dolp, '-o', 'LineWidth', 2);
xlabel('Relative Reflectance (%)');
ylabel('Broadband DoLP');
title('PIMS spectral simulation, 400-700 nm');
grid on;
box on;

figure('Name','Spectral weights');
plot(wavelength_nm, sun_weight, '-o', 'LineWidth', 1.5);
hold on;
plot(wavelength_nm, sky_weight, '-s', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Normalized spectral weight');
legend('Sunlight','Skylight','Location','best');
grid on;
box on;

figure('Name','Mie scattered DoLP vs wavelength');
plot(wavelength_nm, mean_scattered_dolp_by_wavelength, '-o', ...
    'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Mean normalized scattered-light DoLP');
grid on;
box on;

%% ------------------------------------------------------------------------
% 8. Return useful outputs
% -------------------------------------------------------------------------
results = struct();
results.wavelength_nm = wavelength_nm;
results.sun_weight = sun_weight;
results.sky_weight = sky_weight;
results.camera_angles = camera_angles;
results.reflectance = reflectance;
results.model_dolp_unscaled = model_dolp_unscaled;
results.model_dolp = model_dolp;
results.mean_scattered_dolp_by_wavelength = ...
    mean_scattered_dolp_by_wavelength;
results.broadband_scattered_stokes = broadband_scattered_stokes;
results.broadband_skylight_stokes = broadband_skylight_stokes;
results.elapsed_seconds = toc(t_start);

fprintf('Broadband simulation completed in %.2f s.\n', ...
    results.elapsed_seconds);

end

%% ========================================================================
% Local helper functions
% ========================================================================

function y = planck_spectral_shape(lambda_m, temperature_K)
% Relative Planck spectral radiance B_lambda. Absolute scaling is not
% required because the broadband weights are normalized.
h = 6.62607015e-34;
c = 299792458;
kB = 1.380649e-23;

y = (2*h*c^2) ./ (lambda_m.^5) ./ ...
    (exp(h*c ./ (lambda_m*kB*temperature_K)) - 1);
y = y / max(y);
end

function normalized_stokes = normalize_stokes_rows(stokes)
% Normalize each Stokes vector by S0. Rows with negligible S0 are set NaN.
stokes = real(stokes);
s0 = stokes(:,1);
invalid = abs(s0) < 1e-12;
s0(invalid) = NaN;
normalized_stokes = stokes ./ s0;
end

function data_out = select_camera_region( ...
    data_in, azimuth_limits, zenith_limits)
% Apply the same angular-selection helpers used by the repository model.
[data_out,~] = select_data_within_limit_greater_az( ...
    data_in, [], azimuth_limits(1));
[data_out,~] = select_data_within_limit_less_az( ...
    data_out, [], azimuth_limits(2));
[data_out,~] = select_data_within_limit_greater_zen( ...
    data_out, [], zenith_limits(1));
[data_out,~] = select_data_within_limit_less_zen( ...
    data_out, [], zenith_limits(2));
end
