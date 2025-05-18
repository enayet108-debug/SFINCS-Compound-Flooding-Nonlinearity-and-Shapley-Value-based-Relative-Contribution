%% Code developed by Md Enayet Chowdhury, Graduate Research Assistant, UT Austin
clc;
clear;
close all;

%% Loading data
% Flood depth threshold
hmin = 0.3; % Minimum flood depth (meters)

% Load coordinates (x, y)
x = ncread('sfincs_map_WL+Storm.nc', 'x');
y = ncread('sfincs_map_WL+Storm.nc', 'y');

% Load bed level data (zb) to differentiate ocean and overland areas
zb = ncread('sfincs_map_WL+Storm.nc', 'zb');

% Load depths with the correct time index
time_data = ncread('sfincs_map_WL+Storm.nc', 'time'); % Use any file to get time data
time_ref = datetime(2024, 7, 1, 0, 0, 0); % Reference time
time_units = seconds(time_data);
actual_time = time_ref + time_units; % Compute actual time steps

% Function to read water surface data for ocean and overland areas
read_depth = @(ncfile) ...
    combine_ocean_overland(...
        max(ncread(ncfile, 'zsmax'), [], 3), ... % Ocean: zsmax (collapsed to 2D)
        max(ncread(ncfile, 'hmax'), [], 3), ... % Overland: hmax (collapsed to 2D)
        zb, hmin ...
    );

% Custom function to combine ocean and overland data
function combined_data = combine_ocean_overland(zsmax, hmax, zb, hmin)
    combined_data = hmax; % Initialize with hmax (overland)
    combined_data(combined_data < hmin) = NaN; % Apply flood depth threshold
    combined_data(zb < 0) = zsmax(zb < 0); % Replace with zsmax for ocean areas
end

% Load depths from individual files
% Single-driver depths
h_WL = read_depth('sfincs_map_WL_Only.nc');
h_WindPressure = read_depth('sfincs_map_Storm_Only.nc');
h_Dis = read_depth('sfincs_map_Dis_Only.nc');
h_Precipitation = read_depth('sfincs_map_Precipitation_Only.nc');

% Two-driver combinations
h_WL_WindPressure = read_depth('sfincs_map_WL+Storm.nc');
h_WL_Dis = read_depth('sfincs_map_WL+Dis.nc');
h_WL_Precipitation = read_depth('sfincs_map_WL+Precipitation.nc');
h_WindPressure_Dis = read_depth('sfincs_map_Storm+Dis.nc');
h_WindPressure_Precipitation = read_depth('sfincs_map_Storm+Precipitation.nc');
h_Dis_Precipitation = read_depth('sfincs_map_Dis+Precipitation.nc');

% Three-driver combinations
h_WL_WindPressure_Dis = read_depth('sfincs_map_WL+Storm+Dis.nc');
h_WL_Dis_Precipitation = read_depth('sfincs_map_WL+Dis+Precipitation.nc');
h_WL_WindPressure_Precipitation = read_depth('sfincs_map_WL+Storm+Precipitation.nc');
h_WindPressure_Dis_Precipitation = read_depth('sfincs_map_Storm+Dis+Precipitation.nc');

% All four drivers
h_All = read_depth('sfincs_map_WL+Storm+Dis+Precipitation.nc');

%% Calculate Shapley values for each driver

% Define weights for each subset size
weights = [1/4, 1/12, 1/12, 1/4]; % |S| = 0, 1, 2, 3

% Define datasets
datasets = {'h_WL', 'h_WindPressure', 'h_Dis', 'h_Precipitation', ...
            'h_WL_WindPressure', 'h_WL_Dis', 'h_WL_Precipitation', ...
            'h_WindPressure_Dis', 'h_WindPressure_Precipitation', ...
            'h_Dis_Precipitation', 'h_WL_WindPressure_Dis', ...
            'h_WL_Dis_Precipitation', 'h_WL_WindPressure_Precipitation', ...
            'h_WindPressure_Dis_Precipitation', 'h_All'};

% Replace NaN values with 0 in all datasets
for i = 1:length(datasets)
    eval([datasets{i} '(isnan(' datasets{i} ')) = 0;']);
end

% Function to safely subtract two arrays, handling NaNs as specified
function result = safe_subtract(a, b)
    % Replace NaN in only one of the inputs
    a(isnan(a) & ~isnan(b)) = 0;
    b(isnan(b) & ~isnan(a)) = 0;
    % Perform subtraction
    result = a - b;
    % Ensure result is NaN if both a and b are NaN
    result(isnan(a) & isnan(b)) = NaN;
end

% Calculation using weights and safe subtraction
phi_WL = weights(1) * safe_subtract(h_WL, 0) + ...
         weights(2) * (safe_subtract(h_WL_WindPressure, h_WindPressure) + safe_subtract(h_WL_Dis, h_Dis) + safe_subtract(h_WL_Precipitation, h_Precipitation)) + ...
         weights(3) * (safe_subtract(h_WL_WindPressure_Dis, h_WindPressure_Dis) + safe_subtract(h_WL_WindPressure_Precipitation, h_WindPressure_Precipitation) + safe_subtract(h_WL_Dis_Precipitation, h_Dis_Precipitation)) + ...
         weights(4) * safe_subtract(h_All, h_WindPressure_Dis_Precipitation);

% Repeat similar operations for phi_WindPressure, phi_Dis, and phi_Precipitation
phi_WindPressure = weights(1) * safe_subtract(h_WindPressure, 0) + ...
                   weights(2) * (safe_subtract(h_WL_WindPressure, h_WL) + safe_subtract(h_WindPressure_Dis, h_Dis) + safe_subtract(h_WindPressure_Precipitation, h_Precipitation)) + ...
                   weights(3) * (safe_subtract(h_WL_WindPressure_Dis, h_WL_Dis) + safe_subtract(h_WL_WindPressure_Precipitation, h_WL_Precipitation) + safe_subtract(h_WindPressure_Dis_Precipitation, h_Dis_Precipitation)) + ...
                   weights(4) * safe_subtract(h_All, h_WL_Dis_Precipitation);

phi_Dis = weights(1) * safe_subtract(h_Dis, 0) + ...
          weights(2) * (safe_subtract(h_WL_Dis, h_WL) + safe_subtract(h_WindPressure_Dis, h_WindPressure) + safe_subtract(h_Dis_Precipitation, h_Precipitation)) + ...
          weights(3) * (safe_subtract(h_WL_WindPressure_Dis, h_WL_WindPressure) + safe_subtract(h_WL_Dis_Precipitation, h_WL_Precipitation) + safe_subtract(h_WindPressure_Dis_Precipitation, h_WindPressure_Precipitation)) + ...
          weights(4) * safe_subtract(h_All, h_WL_WindPressure_Precipitation);

phi_Precipitation = weights(1) * safe_subtract(h_Precipitation, 0) + ...
                    weights(2) * (safe_subtract(h_WL_Precipitation, h_WL) + safe_subtract(h_WindPressure_Precipitation, h_WindPressure) + safe_subtract(h_Dis_Precipitation, h_Dis)) + ...
                    weights(3) * (safe_subtract(h_WL_WindPressure_Precipitation, h_WL_WindPressure) + safe_subtract(h_WL_Dis_Precipitation, h_WL_Dis) + safe_subtract(h_WindPressure_Dis_Precipitation, h_WindPressure_Dis)) + ...
                    weights(4) * safe_subtract(h_All, h_WL_WindPressure_Dis);

%% Normalized Shapley Values

% Define a threshold
threshold = 0.0; % Values below this threshold will not be shown

% Apply the threshold to exclude low values
phi_WL(phi_WL == threshold) = NaN;
phi_WindPressure(phi_WindPressure == threshold) = NaN;
phi_Dis(phi_Dis == threshold) = NaN;
phi_Precipitation(phi_Precipitation == threshold) = NaN;

denominator = phi_WL + phi_WindPressure + phi_Dis + phi_Precipitation;

Norm_phi_WL = phi_WL ./ denominator;
Norm_phi_WindPressure = phi_WindPressure ./ denominator;
Norm_phi_Dis = phi_Dis ./ denominator;
Norm_phi_Precipitation = phi_Precipitation ./ denominator;

% Define a custom diverging colormap (light colors near zero, dark at extremes)
n_colors = 256; % Total number of colors in the colormap
half_colors = n_colors / 2;

% Red shades for negative values (light red at -1, dark red near 0)
red_colormap = [linspace(1, 0.8, half_colors)', linspace(0.8, 0, half_colors)', linspace(0.8, 0, half_colors)'];

% Blue shades for positive values (light blue near 0, dark blue at +1)
blue_colormap = [linspace(0.8, 0, half_colors)', linspace(0.8, 0, half_colors)', linspace(1, 1, half_colors)'];

% Combine into a diverging colormap
custom_colormap = [flipud(blue_colormap); red_colormap];

% Create a single figure with a 2x2 layout
figure;
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Set colormap to 'jet'
colormap(custom_colormap);

% Water Level (WL)
nexttile;
pcolor(x, y, Norm_phi_WL);
shading interp; % Smooth shading
caxis([-1, 1]); % Global color scale
axis equal; % Ensure consistent aspect ratio
title('Water Level (Normalized)');
set(gca, 'XTick', [], 'YTick', []);

% WindPressure Surge (WindPressure)
nexttile;
pcolor(x, y, Norm_phi_WindPressure);
shading interp; % Smooth shading
caxis([-1, 1]); % Global color scale
axis equal; % Ensure consistent aspect ratio
title('Wind and Pressure (Normalized)');
set(gca, 'XTick', [], 'YTick', []);

% Discharge (Dis)
nexttile;
pcolor(x, y, Norm_phi_Dis);
shading interp; % Smooth shading
caxis([-1, 1]); % Global color scale
axis equal; % Ensure consistent aspect ratio
title('Discharge (Normalized)');
set(gca, 'XTick', [], 'YTick', []);

% Precipitation
nexttile;
pcolor(x, y, Norm_phi_Precipitation);
shading interp; % Smooth shading
caxis([-1, 1]); % Global color scale
axis equal; % Ensure consistent aspect ratio
title('Precipitation (Normalized)');
set(gca, 'XTick', [], 'YTick', []);

% Add a common colorbar for all subplots
cb = colorbar('Location', 'eastoutside'); % Place the colorbar outside the tiles
cb.Layout.Tile = 'east'; % Attach the colorbar to the tiled layout
cb.Label.String = 'Normalized Shapley Value';
cb.Label.FontSize = 14; % Set the font size of the color bar label
cb.Label.FontSize = 12;