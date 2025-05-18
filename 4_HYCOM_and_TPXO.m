%% Code developed by Md Enayet Chowdhury, Graduate Research Assistant, UT Austin
clc;
close all;
clear all;

% Load the data from 'sfincs_wl_output.txt' derived from the HYCOM data preprocessing code
wl_data = load('sfincs_wl_output.txt');
wl_time = wl_data(:, 1);
wl_values = wl_data(:, 2:end);

% Replace NaNs with 0
wl_values(isnan(wl_values)) = 0;

% Load the data from 'sfincs_bzs.txt' derived from Delft Dashboard
bzs_data = load('sfincs_bzs.txt');
bzs_time = bzs_data(:, 1);
bzs_values = bzs_data(:, 2:end) + 0.3; % bias correction of 0.3m. You can change it according to your need.

% Initialize the output matrix for combined water levels
combined_wl = zeros(size(wl_values));

% Loop through each time step in 'sfincs_wl_output.txt'
for i = 1:length(wl_time)
    % Find the index in 'sfincs_bzs.txt' that matches the current time step
    idx = find(bzs_time == wl_time(i));
    
    if ~isempty(idx)
        % If matching time step found, add the values from both files
        combined_wl(i, :) = wl_values(i, :) + bzs_values(idx, :);
    else
        warning('Time step %f not found in sfincs_bzs.txt', wl_time(i));
    end
end

% Save the combined water levels to a new file
output_file = 'combined_water_levels.txt';
fileID = fopen(output_file, 'w');

for i = 1:length(wl_time)
    fprintf(fileID, '%.1f', wl_time(i));  % Write time step
    fprintf(fileID, ' %.4f', combined_wl(i, :));  % Write combined water levels
    fprintf(fileID, '\n');
end

fclose(fileID);
disp('Combined water level data written to combined_water_levels.txt');
