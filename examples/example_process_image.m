%% Process a polarization-mosaic image
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'functions'));

image_file = "example.png";  % Replace with an image path.
result = process_polarization_image(image_file);

figure;
imagesc(result.DoLP);
axis image;
colorbar;
title('DoLP');

figure;
imagesc(result.AoP);
axis image;
colorbar;
title('AoP (deg)');
