%% PIMS DoLP-to-Reflectance Simulation
% Simulates the relationship between mirror reflectance and Degree of
% Linear Polarization (DoLP) using Rayleigh-polarized skylight and
% Mie-scattered sunlight from particles on the mirror surface.
%
% Run this file from the repository root. Required helper functions are
% added to the MATLAB path automatically.

clearvars;
close all;
clc;
tic;

project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'functions'));

%% Model parameters
wavelength_m = 530e-9;
medium_refractive_index = 1.0;
particle_refractive_index = 1.57 + 0i;
wave_number = 2*pi/wavelength_m;

% Measured maximum DoLP used to scale the simulated curve.
max_measured_dolp = 0.6428;

% Fixed source direction used for this data set.
sun_azimuth_deg = 205.3034;
sun_zenith_deg = 69.6823;

% Relative direct-sunlight and skylight contributions.
sun_fraction = 0.82;
sky_fraction = 1 - sun_fraction;

% Camera angular coverage. The default values retain the full simulated map.
camera_azimuth_center_deg = 180;
camera_zenith_center_deg = 45;
camera_azimuth_span_deg = 360;
camera_zenith_span_deg = 90;

camera_azimuth_limits = camera_azimuth_center_deg + ...
    [-0.5, 0.5]*camera_azimuth_span_deg;
camera_zenith_limits = camera_zenith_center_deg + ...
    [-0.5, 0.5]*camera_zenith_span_deg;

% Angular simulation grid.
angle_resolution_deg = 1;
mirror_azimuth_deg = 0:angle_resolution_deg:360;
mirror_zenith_deg = 0:angle_resolution_deg:90;

% Fraction of scattered/reflected contribution accepted by the camera.
% This model parameter is retained as 0.06 for reproducibility.
camera_scatter_acceptance = 0.06;

% Use attenuation-weighted particle-size probability.
use_cross_section_weighting = true;

%% Load particle-size distribution
particle_table = readmatrix( ...
    fullfile(project_root, 'data', 'particle_size_distribution.txt'), ...
    'NumHeaderLines', 1);

% Column 1 is used as the particle-size quantity expected by the model.
% The calculation below follows the diameter convention used by the Mie
% size-parameter calculation. Verify the radius/diameter convention before
% replacing the supplied distribution with a different data set.
particle_distribution = [particle_table(:,1), particle_table(:,3)*250];
particle_distribution = sortrows(particle_distribution, 1);

%% Calculate Rayleigh skylight Stokes parameters
[~,~,~,~,~,~,skylight_zenith_deg,skylight_azimuth_deg,~,~,~,~, ...
    skylight_I,skylight_Q,skylight_U,skylight_V] = ...
    get_skylight_Rayleigh_Scat( ...
        sun_zenith_deg, sun_azimuth_deg, mirror_azimuth_deg, mirror_zenith_deg);

skylight_angles = [skylight_azimuth_deg, skylight_zenith_deg];
skylight_stokes = [skylight_I, skylight_Q, skylight_U, skylight_V];

% The active model treats mirror reflection as preserving the skylight
% Stokes state. Reflection/transmission Mueller-matrix functions remain in
% functions/ for studies that require interface effects.
reflected_skylight_stokes = skylight_stokes;

% Sort by azimuth for consistent angular ordering.
skylight_data = sortrows([skylight_angles, reflected_skylight_stokes], 1);
skylight_angles = skylight_data(:,1:2);
reflected_skylight_stokes = skylight_data(:,3:6);

%% Apply skylight intensity weighting
intensity_weight = reflected_skylight_stokes(:,1) ./ ...
    max(reflected_skylight_stokes(:,1));
reflected_skylight_stokes = reflected_skylight_stokes .* intensity_weight;

%% Select skylight samples within the camera field of view
skylight_camera_data = [skylight_angles, reflected_skylight_stokes];
[skylight_camera_data,~] = select_data_within_limit_greater_az( ...
    skylight_camera_data, [], camera_azimuth_limits(1));
[skylight_camera_data,~] = select_data_within_limit_less_az( ...
    skylight_camera_data, [], camera_azimuth_limits(2));
[skylight_camera_data,~] = select_data_within_limit_greater_zen( ...
    skylight_camera_data, [], camera_zenith_limits(1));
[skylight_camera_data,~] = select_data_within_limit_less_zen( ...
    skylight_camera_data, [], camera_zenith_limits(2));

camera_angles = skylight_camera_data(:,1:2);
camera_skylight_stokes = skylight_camera_data(:,3:6);

%% Precompute Mie-scattered sunlight maps for each particle size
particle_size_nm = particle_distribution(:,1)';
size_parameter = wave_number * medium_refractive_index * ...
    (particle_size_nm*1e-9/2);
incident_sun_stokes = [1; 0; 0; 0];

[I_forward,Q_forward,U_forward,V_forward] = ...
    get_mirror_stokes_map_size_distribution( ...
        mirror_azimuth_deg, mirror_zenith_deg, ...
        sun_azimuth_deg, sun_zenith_deg, ...
        particle_refractive_index, medium_refractive_index, ...
        size_parameter, wave_number, incident_sun_stokes);

% Second geometric contribution used by the model.
[I_mirrored,Q_mirrored,U_mirrored,V_mirrored] = ...
    get_mirror_stokes_map_size_distribution( ...
        mirror_azimuth_deg, mirror_zenith_deg, ...
        sun_azimuth_deg, 180-sun_zenith_deg, ...
        particle_refractive_index, medium_refractive_index, ...
        size_parameter, wave_number, incident_sun_stokes);

%% Sweep particle concentration
concentration_modifier = ...
    [0.0225 0.025 0.0272 0.03 0.04 0.045 0.05 0.06 0.07 0.08 ...
     0.09 0.10 0.12 0.15 0.18 0.20 0.22 0.25 0.28 0.30 ...
     0.40 0.50 0.65 1.00 1.50 2.00 3.00 8.00] / 10;

simulation_area_m2 = (280e-6)^2;
reflectance = zeros(size(concentration_modifier));
average_dolp = zeros(size(concentration_modifier));

for concentration_index = 1:numel(concentration_modifier)
    number_concentration = particle_distribution(:,2)' ./ ...
        concentration_modifier(concentration_index);

    particle_area_m2 = pi*(particle_size_nm*1e-9/2).^2;
    covered_area_m2 = sum(particle_area_m2 .* number_concentration);
    covered_fraction = covered_area_m2 / simulation_area_m2;
    uncovered_fraction = 1 - covered_fraction;

    % Relative reflectance model: uncovered mirror fraction.
    reflectance(concentration_index) = uncovered_fraction;

    %% Particle-size weighting from Mie efficiencies
    attenuation_coefficient = zeros(size(particle_size_nm));
    relative_index = particle_refractive_index / medium_refractive_index;

    for particle_index = 1:numel(particle_size_nm)
        radius_m = particle_size_nm(particle_index)*1e-9/2;
        x = wave_number * radius_m * medium_refractive_index;
        geometric_area_m2 = pi*radius_m^2;
        surface_number_density = number_concentration(particle_index) / ...
            simulation_area_m2;

        mie_result = mie(relative_index, x);
        scattering_efficiency = mie_result(5);
        absorption_efficiency = mie_result(6);

        attenuation_coefficient(particle_index) = ...
            surface_number_density * geometric_area_m2 * ...
            (scattering_efficiency + absorption_efficiency);
    end

    if use_cross_section_weighting
        particle_probability = attenuation_coefficient .* number_concentration;
    else
        particle_probability = number_concentration;
    end
    particle_probability = particle_probability / sum(particle_probability);

    %% Sum particle-size-resolved scattered Stokes maps
    scattered_I = 0;
    scattered_Q = 0;
    scattered_U = 0;
    scattered_V = 0;

    for particle_index = 1:numel(particle_probability)
        weight = particle_probability(particle_index);
        scattered_I = scattered_I + ...
            (I_forward{particle_index,1} + I_mirrored{particle_index,1}) * weight;
        scattered_Q = scattered_Q + ...
            (Q_forward{particle_index,1} + Q_mirrored{particle_index,1}) * weight;
        scattered_U = scattered_U + ...
            (U_forward{particle_index,1} + U_mirrored{particle_index,1}) * weight;
        scattered_V = scattered_V + ...
            (V_forward{particle_index,1} + V_mirrored{particle_index,1}) * weight;
    end

    scattered_stokes = [scattered_I(:), scattered_Q(:), ...
        scattered_U(:), scattered_V(:)];
    mirror_camera_data = [skylight_angles, scattered_stokes];

    [mirror_camera_data,~] = select_data_within_limit_greater_az( ...
        mirror_camera_data, [], camera_azimuth_limits(1));
    [mirror_camera_data,~] = select_data_within_limit_less_az( ...
        mirror_camera_data, [], camera_azimuth_limits(2));
    [mirror_camera_data,~] = select_data_within_limit_greater_zen( ...
        mirror_camera_data, [], camera_zenith_limits(1));
    [mirror_camera_data,~] = select_data_within_limit_less_zen( ...
        mirror_camera_data, [], camera_zenith_limits(2));

    camera_scattered_stokes = mirror_camera_data(:,3:6);

    %% Combine scattered sunlight and reflected skylight in Stokes space
    normalized_scattered = camera_scattered_stokes ./ ...
        camera_scattered_stokes(:,1);
    normalized_skylight = camera_skylight_stokes ./ ...
        camera_skylight_stokes(:,1);

    total_stokes = ...
        normalized_scattered * sun_fraction * covered_fraction + ...
        normalized_skylight * sky_fraction * ...
        (uncovered_fraction + camera_scatter_acceptance*covered_fraction);

    dolp = sqrt(total_stokes(:,2).^2 + total_stokes(:,3).^2) ./ ...
        total_stokes(:,1);
    average_dolp(concentration_index) = mean(dolp);
end

%% Scale simulated curve to measured maximum DoLP
model_dolp = real(average_dolp / max(average_dolp) * max_measured_dolp);

%% Compare model with measured values used for fitting
measured_reflectance = ...
    [0.9852 0.9634 0.9184 0.9890 0.9326 0.9226 0.9523 0.9547 0.8836];
measured_dolp = ...
    [0.5550 0.4813 0.3817 0.6036 0.5247 0.5639 0.5360 0.4670 0.3363];

interpolated_dolp = interp1(reflectance, model_dolp, ...
    measured_reflectance, 'spline');
dolp_error = interpolated_dolp - measured_dolp;
rms_dolp_error = sqrt(mean(dolp_error.^2));

%% Plot model and comparison measurements
figure('Name','PIMS DoLP-Reflectance Model');
plot(reflectance*100, model_dolp, '-o', 'LineWidth', 2);
hold on;
plot(measured_reflectance*100, measured_dolp, '*', 'LineWidth', 2);
xlabel('Relative Reflectance (%)');
ylabel('DoLP');
legend('Model','Measurement','Location','northwest');
grid on;
box on;

fprintf('RMS DoLP error: %.6f\n', rms_dolp_error);
fprintf('Simulation completed in %.2f s.\n', toc);
