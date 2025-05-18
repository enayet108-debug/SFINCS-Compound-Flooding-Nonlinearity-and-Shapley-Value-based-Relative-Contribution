%% Code developed by Md Enayet Chowdhury, Graduate Research Assistant, UT Austin
clc;
clear;
close all;

%% Load coordinates (x, y)
x = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'x');
y = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'y');

% Flood depth threshold
hmin = 0.3; % Minimum flood depth (meters)

% Load bed level data (zb) to differentiate ocean and overland areas
zb = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'zb');

% Load time data
time_data = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'time');
time_ref = datetime(2024, 7, 1, 0, 0, 0); % Reference time
time_units = seconds(time_data); % Ensure this converts the numeric data in 'time_data' to hours
actual_time = time_ref + time_units; % Compute actual time steps

% Define zoomed-in time range
zoom_start = datetime(2024, 7, 6, 0, 0, 0);
zoom_indices = actual_time >= zoom_start;

%% Define UTM coordinates for analysis
% NOAA Stations
noaa_coords = [
    353331.6, 3266110.0;
    307924.2, 3285159.7;
    280921.3, 3290645.8;
    315789.7, 3243084.8;
    325839.2, 3243665.7;
    332552.6, 3248737.7;
    314006.9, 3262699.5
];
noaa_names = {"N1", "N2", "N3", "N4", "N5", "N6", "N7"};

% USGS Stations
usgs_coords = [
    277194.4, 3276214.9;
    262790.0, 3278874.2;
    266445.7, 3287720.1;
    268170.6, 3296280.5;
    276997.2, 3311990.1;
    274837.2, 3336670.1;
    308275.9, 3317571.9;
    325009.1, 3317361.9
];
usgs_names = {"D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8"};

% Selected Locations
selected_coords = [
    311449.18, 3251615.78;
    280082.43, 3291275.66;
    325147.58, 3243380.27;
    268006.13, 3292121.87;
    268703.26, 3294665.56;
    267640.9901, 3288800.206
];
selected_names = {"L1", "L2", "L3", "L4", "L5", "L6"};

% Combine all station data
station_groups = {noaa_coords, usgs_coords, selected_coords};
station_names = {noaa_names, usgs_names, selected_names};
group_titles = {"NOAA Stations", "USGS Stations", "Selected Locations"};

%% Helper Functions
function combined_data = read_depth(ncfile, time_idx, zb, hmin)
    zs = ncread(ncfile, 'zs', [1, 1, time_idx], [Inf, Inf, 1]); % Ocean: zs for the specified time step
    h = ncread(ncfile, 'h', [1, 1, time_idx], [Inf, Inf, 1]); % Overland: h for the specified time step
    combined_data = combine_ocean_overland(zs, h, zb, hmin);
end

function combined_data = combine_ocean_overland(zs, h, zb, hmin)
    combined_data = h; % Initialize with h (overland)
    combined_data(zb < 0) = zs(zb < 0); % Replace with zs for ocean areas
end

function result = safe_subtract(a, b)
    % Replace NaN in only one of the inputs
    a(isnan(a) & ~isnan(b)) = 0;
    b(isnan(b) & ~isnan(a)) = 0;
    % Perform subtraction
    result = a - b;
    % Ensure result is NaN if both a and b are NaN
    result(isnan(a) & isnan(b)) = NaN;
end

% Set all fonts to Arial
set(groot, 'DefaultAxesFontName', 'Arial'); % For axes
set(groot, 'DefaultTextFontName', 'Arial'); % For text elements

%% Calculate Nonlinearity Time Series for Each Group and Plot
figure;
tiledlayout(2, 3); % Layout for two rows and three columns

for group_idx = 1:length(station_groups)
    coords = station_groups{group_idx};
    names = station_names{group_idx};
    num_locations = size(coords, 1);
    num_time_steps = length(actual_time);

    % Recalculate nonlinearity_values for the current group
    nonlinearity_values = zeros(num_locations, num_time_steps);

    for t = 1:num_time_steps
        % Load linear components for each driver
        h_WL = read_depth('sfincs_map_WL_Only.nc', t, zb, hmin);
        h_Dis = read_depth('sfincs_map_Dis_Only.nc', t, zb, hmin);
        h_Storm = read_depth('sfincs_map_Storm_Only.nc', t, zb, hmin);
        h_Precipitation = read_depth('sfincs_map_Precipitation_Only.nc', t, zb, hmin);

        % Calculate linearly added water levels
        h_linear_combined = h_WL + h_Dis + h_Storm + h_Precipitation;

        % Load nonlinear water levels
        h_nonlinear = read_depth('sfincs_map_WL+Storm+Dis+Precipitation.nc', t, zb, hmin);

        % Calculate nonlinearity
        h_nonlinearity = safe_subtract(h_nonlinear, h_linear_combined);

        % Store nonlinearity values at specified locations
        for loc_idx = 1:num_locations
            distances = hypot(x - coords(loc_idx, 1), y - coords(loc_idx, 2));
            [~, idx] = min(distances(:));
            [i, j] = ind2sub(size(x), idx);
            nonlinearity_values(loc_idx, t) = h_nonlinearity(i, j);
        end
    end

    % Plot all stations in one subplot for the zoomed-in time range
    nexttile(group_idx + 3); % Shift to the second row of the layout
    hold on;

    for loc_idx = 1:num_locations
        plot(actual_time(zoom_indices), nonlinearity_values(loc_idx, zoom_indices), 'LineWidth', 1.0, ...
             'DisplayName', names{loc_idx});
    end

    % Customize the zoomed-in plot
    xlabel('Date', 'FontSize', 10);
    if group_idx == 1
        ylabel('Marginal Water Contribution (m)', 'FontSize', 10); % Only for the leftmost plot
    end
    % title([group_titles{group_idx}, ' (Zoomed)'], 'FontSize', 14);
    legend('show', 'Location', 'best', 'FontSize', 8); % Add legend for each group
    grid on;
    hold off;
end
