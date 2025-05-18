%% This code was primarily developed by Dr Wonhyun Lee, Research Assistant Professor, UT Austin
% clear all the previous data
clc;
close all;
clear all;

%% SECTION 1. Loading the HYCOM data

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ========================================================================
% BETWEEN CODES, USER NEED TO INPUT or EDIT FILE NAMES OR PATH
% ========================================================================
disp 'Handle with HYCOM .nc files --> .mat files'
% Loading HYCOM after downloading the HYCOM data (non-tidal water level)
% from the website
% ========================================================================
fpath=('ssh_2024_HYCOM.nc4');
% ========================================================================
% Display the contents of the NetCDF file
ncdisp(fpath);

%% SECTION 2. Setting up the time frame

clearvars -except header fpath
close all;
clc;

% Read the necessary variables from the NetCDF file
time = ncread(fpath, 'time');
lat = ncread(fpath, 'lat');
lon = ncread(fpath, 'lon');
surf_el = ncread(fpath, 'surf_el');

% Reference Time
Refer_year = 2024;  Refer_month = 06; Refer_day = 25;
Refer_hour = 00;    Refer_min = 00;   Refer_sec = 00;

% Start Time
Start_year = 2024;  Start_month = 06; Start_day = 25;
Start_hour = 00;    Start_min = 00;   Start_sec = 00;

End_year = 2024;    End_month = 07;   End_day = 08;
End_hour = 00;      End_min = 00;     End_sec = 00;

% Reference Time
ttr = datenum([Refer_year, Refer_month, Refer_day, Refer_hour, Refer_min, Refer_sec]);
% Start Time
tt0 = datenum([Start_year, Start_month, Start_day, Start_hour, Start_min, Start_sec]);
% End Time
tt1 = datenum([End_year, End_month, End_day, End_hour, End_min, End_sec]);

% Time Step (seconds)
ddt_hr = 3;  % Time step in hours
ddt = ddt_hr * 60 * 60;  % Time step in seconds
% ========================================================================

%% SECTION 3. Extracting the Latitude, Longitude, and Sea Surface Elevation data

Len_lat = length(lat);  % 'lat' read from the NetCDF file using ncread
Len_lon = length(lon);  % 'lon' read from the NetCDF file using ncread
Lat_dy = abs(lat(1) - lat(2));  % Calculate latitude spacing
Lon_dx = abs(lon(1) - lon(2));  % Calculate longitude spacing

% Meshgrid for lat/lon coordinates
[X_Lon, Y_Lat] = meshgrid(lon, lat);

% Extract the surface elevation (water level) data
levels = surf_el(:,:,1);  % Assuming the first time step for illustration

% Set the long_name (as indicated in the comment)
long_name = 'waterlevel';

% TNNC is the length of the NetCDF content in terms of time steps
TNNC = length(time);  % 'time' is read from the NetCDF file using ncread
disp(TNNC);  % This will print the total number of time steps

% Get the size of the Sur_El matrix
[TSz_lon, TSz_lat, TSz_Tim] = size(surf_el);
% surf0_el = double(surf_el);
surf3_el =permute(surf_el,[2,1,3]);

%% SECTION 4. Loading the SFINCS boundary data

% Load boundary coordinates from sfincs_bnd.txt. You will get this file
% after downloading the TPXO data from the Delft Dashboard. Then it needs
% to be converted into a normal text file.
boundary_coords = load('sfincs_bnd.txt');

% Prepare output matrix for interpolated water levels
output_wl = zeros(TNNC, size(boundary_coords, 1));

%% Function to convert the coordinate system

function  [Lat,Lon] = utm2deg(xx,yy,utmzone)
% Argument checking
%
error(nargchk(3, 3, nargin)); %3 arguments required
n1=length(xx);
n2=length(yy);
n3=size(utmzone,1);
if (n1~=n2 || n1~=n3)
   error('x,y and utmzone vectors should have the same number or rows');
end
c=size(utmzone,2);
if (c~=4)
   error('utmzone should be a vector of strings like "30 T"');
end

% Memory pre-allocation
%
Lat=zeros(n1,1);
Lon=zeros(n1,1);

% Main Loop
%
for i=1:n1
   if (utmzone(i,4)>'X' || utmzone(i,4)<'C')
      fprintf('utm2deg: Warning utmzone should be a vector of strings like "30 T", not "30 t"\n');
   end
   if (utmzone(i,4)>'M')
      hemis='N';   % Northern hemisphere
   else
      hemis='S';
   end

   x=xx(i);
   y=yy(i);
   zone=str2double(utmzone(i,1:2));

   sa = 6378137.000000 ; sb = 6356752.314245;
  
%   e = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sa;
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
%   alpha = ( sa - sb ) / sa;             %f
%   ablandamiento = 1 / alpha;   % 1/f

   X = x - 500000;
   
   if hemis == 'S' || hemis == 's'
       Y = y - 10000000;
   else
       Y = y;
   end
    
   S = ( ( zone * 6 ) - 183 ); 
   lat =  Y / ( 6366197.724 * 0.9996 );                                    
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   a = X / v;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   b = ( Y - Bm ) / v;
   Epsi = ( ( e2cuadrada * a^ 2 ) / 2 ) * ( cos(lat) )^ 2;
   Eps = a * ( 1 - ( Epsi / 3 ) );
   nab = ( b * ( 1 - Epsi ) ) + lat;
   senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
   Delt = atan(senoheps / (cos(nab) ) );
   TaO = atan(cos(Delt) * tan(nab));
   longitude = (Delt *(180 / pi ) ) + S;
   latitude = ( lat + ( 1 + e2cuadrada* (cos(lat)^ 2) - ( 3 / 2 ) * e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) ) * ( TaO - lat ) ) * ...
                    (180 / pi);
   
   Lat(i)=latitude;
   Lon(i)=longitude;
end
   
end
%% SECTION 5. Coordinate conversion

% Load boundary coordinates from sfincs_bnd.txt
boundary_coords = load('sfincs_bnd.txt');

% Generate a UTM zone array for all points. Assign your area's correct UTM zone here
utmZone = repmat('16 N', size(boundary_coords, 1), 1);  % Create a column of '15 N'

% Convert UTM coordinates to latitude/longitude (WGS 84)
[bnd_lat, bnd_lon] = utm2deg(boundary_coords(:, 1), boundary_coords(:, 2), utmZone);  % Convert UTM to Lat/Lon

% Convert boundary longitudes from degrees west to degrees east
bnd_lon_east = 360 + bnd_lon;  % Convert degrees west to degrees east

%% SECTION 6. Generating the output file

% Prepare output matrix for interpolated water levels
output_wl = zeros(TNNC, size(boundary_coords, 1));

% Interpolate Water Level for Each Boundary Point at Each Time Step
for t_idx = 1:TNNC
    % Get water level data for the current time step
    current_wl_data = surf3_el(:, :, t_idx);  % Assuming 'surf_el' contains water level data

    % Interpolate water level at each boundary point
    for b_idx = 1:size(boundary_coords, 1)
        % Perform spatial interpolation using griddata
        output_wl(t_idx, b_idx) = griddata(X_Lon, Y_Lat, current_wl_data, bnd_lon_east(b_idx), bnd_lat(b_idx), 'linear');
    end
end

% Write the interpolated water levels to a new file in a format like 'sfincs_bzs.txt'
fileID = fopen('sfincs_wl_output.txt', 'w');
for t_idx = 1:TNNC
    fprintf(fileID, '%.1f', (t_idx-1) * ddt);  % Time in seconds
    fprintf(fileID, ' %.4f', output_wl(t_idx, :));  % Water levels
    fprintf(fileID, '\n');
end
fclose(fileID);

disp('Water level data for each boundary point and time step has been written to sfincs_wl_output.txt. Now you can combine it with the TPXO data.');
