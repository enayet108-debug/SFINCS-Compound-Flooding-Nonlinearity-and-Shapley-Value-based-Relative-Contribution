%% Code developed by Md Enayet Chowdhury, Graduate Research Assistant, UT Austin
clc;
clear;
close all;

%% Load coordinates (x, y)
x = ncread('sfincs_map_WL+Storm.nc', 'x');
y = ncread('sfincs_map_WL+Storm.nc', 'y');

% Flood depth threshold
hmin = 0.3; % Minimum flood depth (meters)

% Load bed level data (zb) to differentiate ocean and overland areas
zb = ncread('sfincs_map_WL+Storm.nc', 'zb');

% Load time data
time_data = ncread('sfincs_map_WL+Storm.nc', 'time');
time_ref = datetime(2024, 7, 1, 0, 0, 0); % Reference time
time_units = seconds(time_data); % Ensure this converts the numeric data in 'time_data' to hours
actual_time = time_ref + time_units; % Compute actual time steps

utm_coords = [
    353331.6, 3266110.0;
    307924.2, 3285159.7;
    280921.3, 3290645.8;
    315789.7, 3243084.8;
    325839.2, 3243665.7;
    332552.6, 3248737.7;
    314006.9, 3262699.5;
    277194.4, 3276214.9;
    262790.0, 3278874.2;
    266445.7, 3287720.1;
    268170.6, 3296280.5;
    276997.2, 3311990.1;
    274837.2, 3336670.1;
    308275.9, 3317571.9;
    325009.1, 3317361.9;
    311449.18, 3251615.78;
    280082.43, 3291275.66;
    325147.58, 3243380.27;
    268006.13, 3292121.87;
    268703.26, 3294665.56;
    267640.9901, 3288800.206;
];

%% Local Functions
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

num_locations = size(utm_coords, 1);
num_time_steps = length(actual_time);
phi_values = zeros(num_locations, num_time_steps, 5); % 4 drivers + h_All

for t = 1:num_time_steps
    % Load depths from individual and combination files for each driver at time t
    h_WL = read_depth('sfincs_map_WL_Only.nc', t, zb, hmin);
    h_WindPressure = read_depth('sfincs_map_Storm_Only.nc', t, zb, hmin);
    h_Dis = read_depth('sfincs_map_Dis_Only.nc', t, zb, hmin);
    h_Precipitation = read_depth('sfincs_map_Precipitation_Only.nc', t, zb, hmin);

    % Two-driver combinations
    h_WL_WindPressure = read_depth('sfincs_map_WL+Storm.nc', t, zb, hmin);
    h_WL_Dis = read_depth('sfincs_map_WL+Dis.nc', t, zb, hmin);
    h_WL_Precipitation = read_depth('sfincs_map_WL+Precipitation.nc', t, zb, hmin);
    h_WindPressure_Dis = read_depth('sfincs_map_Storm+Dis.nc', t, zb, hmin);
    h_WindPressure_Precipitation = read_depth('sfincs_map_Storm+Precipitation.nc', t, zb, hmin);
    h_Dis_Precipitation = read_depth('sfincs_map_Dis+Precipitation.nc', t, zb, hmin);

    % Three-driver combinations
    h_WL_WindPressure_Dis = read_depth('sfincs_map_WL+Storm+Dis.nc', t, zb, hmin);
    h_WL_Dis_Precipitation = read_depth('sfincs_map_WL+Dis+Precipitation.nc', t, zb, hmin);
    h_WL_WindPressure_Precipitation = read_depth('sfincs_map_WL+Storm+Precipitation.nc', t, zb, hmin);
    h_WindPressure_Dis_Precipitation = read_depth('sfincs_map_Storm+Dis+Precipitation.nc', t, zb, hmin);

    % All four drivers
    h_All = read_depth('sfincs_map_WL+Storm+Dis+Precipitation.nc', t, zb, hmin);

    % Define weights for each subset size
    weights = [1/4, 1/12, 1/12, 1/4]; % |S| = 0, 1, 2, 3

    % Calculation using weights and safe subtraction
    phi_WL = weights(1) * safe_subtract(h_WL, 0) + ...
             weights(2) * (safe_subtract(h_WL_WindPressure, h_WindPressure) + safe_subtract(h_WL_Dis, h_Dis) + safe_subtract(h_WL_Precipitation, h_Precipitation)) + ...
             weights(3) * (safe_subtract(h_WL_WindPressure_Dis, h_WindPressure_Dis) + safe_subtract(h_WL_WindPressure_Precipitation, h_WindPressure_Precipitation) + safe_subtract(h_WL_Dis_Precipitation, h_Dis_Precipitation)) + ...
             weights(4) * safe_subtract(h_All, h_WindPressure_Dis_Precipitation);

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

    % Store Shapley values and h_All in the phi_values matrix for each location
    for loc_idx = 1:num_locations
        % Calculate the distances from the current UTM coordinates to each grid point
        distances = hypot(x - utm_coords(loc_idx, 1), y - utm_coords(loc_idx, 2));

        % Find the index of the minimum distance
        [~, idx] = min(distances(:)); % Flatten the array to get a linear index

        % Convert linear index to subscript indices
        [i, j] = ind2sub(size(x), idx); % Ensure x and y are appropriately sized

        % Assign values using the corrected indices
        phi_values(loc_idx, t, 1) = phi_WL(i, j);
        phi_values(loc_idx, t, 2) = phi_WindPressure(i, j);
        phi_values(loc_idx, t, 3) = phi_Dis(i, j);
        phi_values(loc_idx, t, 4) = phi_Precipitation(i, j);
        phi_values(loc_idx, t, 5) = h_All(i, j);
    end
end

%% Plot with 21 Locations and h_All (from July 6 to end)
% Define distinct colors for drivers and h_All
driver_titles = {'Water Level', 'Wind and Pressure', 'Discharge', 'Precipitation', 'h\_All'};
location_names = {'N1', 'N2', 'N3', 'N4', 'N5', 'N6','N7', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6','D7', 'D8', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6'};
driver_colors = {'b', 'g', 'c', 'k'}; % Blue, Green, Cyan, Black for drivers
h_All_color = 'r'; % Red exclusively for h_All

% Ensure Proper Date Indexing
date_indices = find(actual_time >= datetime(2024, 7, 6)); % Indices from July 6 onward
valid_date_indices = date_indices(date_indices <= size(phi_values, 2)); % Ensure within bounds
filtered_time = actual_time(valid_date_indices); % Extract corresponding times

% Create a figure and set size
figure;
set(gcf, 'Position', [100, 100, 1000, 2000]); % Increase figure size

% Create a tiled layout with adjustable spacing
t = tiledlayout(7, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Adjust gaps between subplots

% Loop Through Each Location
for loc_idx = 1:21
    nexttile;
    hold on;

    % Loop through each flood driver (Water Level, Wind & Pressure, Discharge, Precipitation)
    for k = 1:4
        % Check if valid indices exist before indexing
        if ~isempty(valid_date_indices)
            filtered_phi = squeeze(phi_values(loc_idx, valid_date_indices, k));

            % Remove zero values for better visualization
            non_zero_indices = find(filtered_phi ~= 0);
            if ~isempty(non_zero_indices)
                plot_start = min(non_zero_indices);
                plot_end = max(non_zero_indices);
                plot(filtered_time(plot_start:plot_end), filtered_phi(plot_start:plot_end), ...
                     'Color', driver_colors{k}, 'LineWidth', 0.8, 'DisplayName', driver_titles{k});
            end
        end
    end

    % Plot h_All (Combination of all drivers)
    if ~isempty(valid_date_indices)
        filtered_h_All = squeeze(phi_values(loc_idx, valid_date_indices, 5)); % h_All is at index 5
        if any(filtered_h_All ~= 0)
            plot(filtered_time, filtered_h_All, 'Color', h_All_color, 'LineWidth', 0.8, ...
                 'DisplayName', 'h (WL, WP, Dis, Precip)');
        end
    end

    hold off;

    % Add title for each subplot
    title(location_names{loc_idx}, 'FontSize', 12); 
    grid on;
    axis tight;

    % Set font size for axes numbers
    ax = gca;
    ax.FontSize = 10; % Increase font size for x and y axes numbers

    % Add X-label only in the last row
    if loc_idx > 18 % Last row (since there are 21 plots arranged in a 7x3 layout)
        xlabel('Date');
    else
        ax.XTickLabel = []; % Remove x-axis labels for all other rows
    end
end

% Add a single legend
lgd = legend({'Water Level', 'Wind and Pressure', 'Discharge', 'Precipitation', 'h (WL, WP, Dis, Precip)'}, ...
    'Location', 'northeastoutside', 'Orientation', 'horizontal', 'FontSize', 10);
lgd.Layout.Tile = 'south';

% Save the figure as a PNG file with 300 DPI
saveas(gcf, 'Shapley_Values_Plot.png');
print(gcf, 'Shapley_Values_Plot', '-dpng', '-r300'); % 300 DPI resolution
