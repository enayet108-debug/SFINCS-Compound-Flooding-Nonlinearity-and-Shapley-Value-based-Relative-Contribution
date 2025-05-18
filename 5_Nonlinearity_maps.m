%% Code developed by Md Enayet Chowdhury, Graduate Research Assistant, UT Austin
clc;
clear;
close all;

%% File names for combined and individual drivers after running SFINCS
file_combinations = {
    'sfincs_map_WL+Storm.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Storm_Only.nc'};
    'sfincs_map_WL+Dis.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Dis_Only.nc'};
    'sfincs_map_WL+Precipitation.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Precipitation_Only.nc'};
    'sfincs_map_Storm+Precipitation.nc', {'sfincs_map_Storm_Only.nc', 'sfincs_map_Precipitation_Only.nc'};
    'sfincs_map_Storm+Dis.nc', {'sfincs_map_Storm_Only.nc', 'sfincs_map_Dis_Only.nc'};
    'sfincs_map_Dis+Precipitation.nc', {'sfincs_map_Dis_Only.nc', 'sfincs_map_Precipitation_Only.nc'};
    'sfincs_map_WL+Storm+Dis.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Storm_Only.nc', 'sfincs_map_Dis_Only.nc'};
    'sfincs_map_WL+Dis+Precipitation.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Dis_Only.nc', 'sfincs_map_Precipitation_Only.nc'};
    'sfincs_map_WL+Storm+Precipitation.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Storm_Only.nc', 'sfincs_map_Precipitation_Only.nc'};
    'sfincs_map_Storm+Dis+Precipitation.nc', {'sfincs_map_Storm_Only.nc', 'sfincs_map_Dis_Only.nc', 'sfincs_map_Precipitation_Only.nc'};
    'sfincs_map_WL+Storm+Dis+Precipitation.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Storm_Only.nc', 'sfincs_map_Dis_Only.nc', 'sfincs_map_Precipitation_Only.nc'};   
};

% Define readable titles for reordered file combinations
titles = {
    'sfincs_map_WL+Storm.nc', '(a) WL, WP';
    'sfincs_map_WL+Dis.nc', '(b) WL, Dis';
    'sfincs_map_WL+Precipitation.nc', '(c) WL, Precip';
    'sfincs_map_Storm+Precipitation.nc', '(d) WP, Precip';
    'sfincs_map_Storm+Dis.nc', '(e) WP, Dis';
    'sfincs_map_Dis+Precipitation.nc', '(f) Dis, Precip';
    'sfincs_map_WL+Storm+Dis.nc', '(g) WL, WP, Dis';
    'sfincs_map_WL+Dis+Precipitation.nc', '(h) WL, Dis, Precip';
    'sfincs_map_WL+Storm+Precipitation.nc', '(i) WL, WP, Precip';
    'sfincs_map_Storm+Dis+Precipitation.nc', '(j) WP, Dis, Precip';
    'sfincs_map_WL+Storm+Dis+Precipitation.nc', '(k) WL, WP, Dis, Precip';
};

% Flood depth threshold
hmin = 0.3; % Minimum flood depth (meters)

% Reference time
time_ref = datetime(2024, 7, 1, 0, 0, 0);

% Precompute global min and max for clim
global_min = inf;
global_max = -inf;

% Load bed level data (zb) to differentiate ocean and overland areas
zb = ncread('sfincs_map_WL_Only.nc', 'zb');

% Calculate D_difference for all files to find global min and max
D_differences = cell(size(file_combinations, 1), 1); % Store D_difference for plotting
x = [];
y = [];
for i = 1:size(file_combinations, 1)
    combined_file = file_combinations{i, 1};
    individual_files = file_combinations{i, 2};

    % Read zsmax and hmax from the combined file
    zsmax = ncread(combined_file, 'zsmax'); % Maximum water level
    hmax = ncread(combined_file, 'hmax'); % Maximum water depth

    % Collapse the third dimension for both variables
    zsmax = max(zsmax, [], 3); % Collapse to 2D
    hmax = max(hmax, [], 3); % Collapse to 2D
    
    % Combine zsmax and hmax based on zb values
    h_combined = hmax; % Initialize with hmax
    h_combined(h_combined < hmin) = NaN; % Apply flood depth threshold
    h_combined(zb < 0) = zsmax(zb < 0); % Replace with zsmax where zb < 0 (ocean areas)
    
    % Initialize total disturbance for linear addition
    D_total = [];
    for j = 1:numel(individual_files)
        ncfile = individual_files{j};

        % Read zsmax and hmax from individual files
        zsmax = ncread(ncfile, 'zsmax'); % Maximum water level
        hmax = ncread(ncfile, 'hmax'); % Maximum water depth

        % Collapse the third dimension for both variables
        zsmax = max(zsmax, [], 3); % Collapse to 2D
        hmax = max(hmax, [], 3); % Collapse to 2D
        
        % Combine zsmax and hmax based on zb values
        h = hmax; % Initialize with hmax
        h(h < hmin) = NaN; % Apply flood depth threshold
        h(zb < 0) = zsmax(zb < 0); % Replace with zsmax where zb < 0 (ocean areas)
        
        % Initialize D_total on the first iteration
        if j == 1
            D_total = zeros(size(h)); % Initialize total disturbance
            if isempty(x) || isempty(y) % Only read x and y once
                x = ncread(ncfile, 'x');
                y = ncread(ncfile, 'y');
            end
        end
        
        % Add to total disturbance (ignoring NaN values)
        D_total = nansum(cat(3, D_total, h), 3);
    end
    
    % Set all values of D_total < hmin to NaN
    % D_total(D_total < hmin) = NaN;
    D_total(zb > 0 & D_total < hmin) = NaN;
 
    % Create modified copies where NaNs are replaced with zeros only if the other is not NaN
    D_total_mod = D_total;
    h_combined_mod = h_combined;
    
    % Find indices where D_total is NaN and h_combined is not NaN
    idx1 = isnan(D_total) & ~isnan(h_combined);
    % Replace NaNs in D_total_mod with zeros at these indices
    D_total_mod(idx1) = 0;
    
    % Find indices where D_total is not NaN and h_combined is NaN
    idx2 = ~isnan(D_total) & isnan(h_combined);
    % Replace NaNs in h_combined_mod with zeros at these indices
    h_combined_mod(idx2) = 0;
    
    % Compute D_difference
    D_difference = D_total_mod - h_combined_mod;
    
    % For indices where both D_total and h_combined are NaN, D_difference should remain NaN
    idx_both_nan = isnan(D_total) & isnan(h_combined);
    D_difference(idx_both_nan) = NaN;
    
    % Store D_difference for plotting
    D_differences{i} = D_difference;
    
    % Update global min and max
    global_min = min(global_min, nanmin(D_difference(:)));
    global_max = max(global_max, nanmax(D_difference(:)));
end

%% Nonlinearity map plotting

% Define a custom diverging colormap (light colors near zero, dark at extremes)
n_colors = 256; % Total number of colors in the colormap
half_colors = n_colors / 2;

% Red shades for negative values (light red at -1, dark red near 0)
red_colormap = [linspace(1, 0.8, half_colors)', linspace(0.8, 0, half_colors)', linspace(0.8, 0, half_colors)'];

% Blue shades for positive values (light blue near 0, dark blue at +1)
blue_colormap = [linspace(0.8, 0, half_colors)', linspace(0.8, 0, half_colors)', linspace(1, 1, half_colors)'];

% Combine into a diverging colormap
custom_colormap = [flipud(blue_colormap); red_colormap];

% Manually set color bar limits
manual_min = -0.5; % Replace with your desired minimum value
manual_max = 0.5;  % Replace with your desired maximum value

% Plot all difference maps
figure;

% Set the colormap for the entire figure
colormap(custom_colormap);

t = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(file_combinations, 1)
    nexttile;
    flood_map_diff = pcolor(x, y, D_differences{i}); % Plot the difference map
    set(flood_map_diff, 'EdgeColor', 'none'); % Remove grid lines for clarity
    shading interp;
    caxis([manual_min, manual_max]); % Use manually set color bar limits

    % Get readable title from the mapping
    file_name = file_combinations{i, 1};
    title_index = strcmp(titles(:, 1), file_name); % Find the matching title
    readable_title = titles{title_index, 2}; % Retrieve the simplified title

    % Set the title for the current plot
    title(readable_title, 'Interpreter', 'none', 'FontSize', 12);

    % Remove x and y ticks and labels, but keep the black box (axis outline)
    set(gca, 'XTick', [], 'YTick', [], 'TickLength', [0 0]); % Remove ticks
    axis equal; % Keep the aspect ratio equal
    box on; % Ensure the black box outline is visible
end

% Add a common colorbar for all subplots
cb = colorbar('Location', 'southoutside'); % Place the colorbar outside the tiles
cb.Layout.Tile = 'south'; % Attach the colorbar to the tiled layout
cb.Label.String = 'Nonlinearity (m)'; % Update the label to your needs
cb.Label.FontSize = 14; % Set the font size of the color bar label
cb.FontSize = 12; % Set the font size of the color bar tick labels
