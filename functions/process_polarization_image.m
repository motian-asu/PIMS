function result = process_polarization_image(image_file, roi)
%PROCESS_POLARIZATION_IMAGE Calculate linear Stokes parameters, DoLP, and AoP.
%
% result = process_polarization_image(image_file)
% result = process_polarization_image(image_file, roi)
%
% INPUTS
%   image_file : path to a raw image from a 2x2 polarization-mosaic camera.
%   roi        : optional [x y width height] rectangle for DoLP statistics.
%
% OUTPUTS
%   result.I0, I45, I90, I135 : analyzer sub-images.
%   result.S0, S1, S2          : normalized linear Stokes quantities.
%   result.DoLP                : degree of linear polarization.
%   result.AoP                 : angle of polarization in degrees, [0,180).
%   result.roi_stats           : summary statistics when roi is supplied.
%
% Pixel mapping:
%   I0   = odd rows, odd columns
%   I45  = odd rows, even columns
%   I90  = even rows, even columns
%   I135 = even rows, odd columns

arguments
    image_file (1,1) string
    roi (1,4) double = [NaN NaN NaN NaN]
end

image_data = imread(image_file);
image_data = double(image_data) / 255;

if ndims(image_data) == 3
    image_data = mean(image_data, 3);
end

[n_rows, n_cols] = size(image_data);
n_rows = n_rows - mod(n_rows, 2);
n_cols = n_cols - mod(n_cols, 2);
image_data = image_data(1:n_rows, 1:n_cols);

I0   = image_data(1:2:end, 1:2:end);
I45  = image_data(1:2:end, 2:2:end);
I90  = image_data(2:2:end, 2:2:end);
I135 = image_data(2:2:end, 1:2:end);

S0_raw = (I45 + I135 + I0 + I90) / 2;
S1_raw = I0 - I90;
S2_raw = I45 - I135;

S0 = S0_raw ./ S0_raw;
S1 = S1_raw ./ S0_raw;
S2 = S2_raw ./ S0_raw;

DoLP = sqrt(S1.^2 + S2.^2) ./ S0;

% Convert the Stokes angle to the physical linear-polarization angle.
% The quadrant handling is kept explicit to preserve the analysis logic.
AoP2 = atan(S2 ./ S1) / pi * 180;
AoP2(S1 < 0) = AoP2(S1 < 0) + 180;
AoP2(AoP2 < 0) = 360 + AoP2(AoP2 < 0);
AoP2(S1 == 0 | S0_raw == 0) = 0;
AoP = 0.5 * AoP2;

result = struct( ...
    'I0', I0, ...
    'I45', I45, ...
    'I90', I90, ...
    'I135', I135, ...
    'S0', S0, ...
    'S1', S1, ...
    'S2', S2, ...
    'DoLP', DoLP, ...
    'AoP', AoP, ...
    'roi_stats', []);

if all(isfinite(roi))
    x1 = max(1, round(roi(1)));
    y1 = max(1, round(roi(2)));
    x2 = min(size(DoLP,2), round(roi(1) + roi(3)));
    y2 = min(size(DoLP,1), round(roi(2) + roi(4)));

    values = DoLP(y1:y2, x1:x2);
    result.roi_stats = struct( ...
        'median', median(values, 'all', 'omitnan'), ...
        'min', min(values, [], 'all', 'omitnan'), ...
        'max', max(values, [], 'all', 'omitnan'), ...
        'mean', mean(values, 'all', 'omitnan'), ...
        'std', std(values, 0, 'all', 'omitnan'));
end
end
