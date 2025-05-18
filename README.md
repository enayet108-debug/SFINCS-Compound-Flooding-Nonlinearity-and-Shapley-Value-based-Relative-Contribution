# Compound Flooding using SFINCS: Nonlinearity and Shapley Value-based Relative Contribution
This repository provides the complete workflow and analysis codes for quantifying *spatiotemporal nonlinearity* and *Shapley-value based relative contribution* of four compound flood drivers: 
- Offshore Water Level (WL)
- Wind-Pressure (WP)
- River Discharge (Dis)
- Precipitation (Precip). 

The case study is based on *Hurricane Beryl (2024)* in the *Houston–Galveston region, Texas, USA*, using numerical simulations with the **Super-Fast INundation of CoastS (SFINCS)** model.

---
## 📁 Project Structure

```text
SFINCS-Compound-Flooding-Nonlinearity-and-Shapley-Value-based-Relative-Contribution/
├── notebooks/                     # Python notebooks for preprocessing
│   ├── HydroMT_Preprocessing_for_SFINCS.ipynb
│   └── Precipitation_Code.ipynb
│
├── matlab/                        # MATLAB scripts for analysis
│   ├── 1_Preprocessing_HYCOM_Data.m
│   ├── 2_HYCOM_and_TPXO.m
│   ├── 3_Nonlinearity_maps.m
│   ├── 4_Nonlinearity_timeseries.m
│   ├── 5_Shapley_Values_maps.m
│   └── 6_Shapley_values_timeseries.m
│
├── LICENSE                       # Open-source license
└── README.md                     # Project overview

```
---
## 📥 Input Data

| Input Type         | Source / Model                                        | Resolution / Format             | Notes                                                                 |
|--------------------|--------------------------------------------------------|----------------------------------|-----------------------------------------------------------------------|
| **Digital Elevation Model (DEM)** | NOAA CUDEM + CRM                                | 3 m (coast), 10 m (offshore), 30 m (beyond CUDEM) | Integrated topobathy DEMs; 5-pixel buffer used           |
| **Manning’s Roughness** | NLCD (2016)                                          | 30 m                             | Mapped from 16 land cover classes                                     |
| **Infiltration (Curve Number)** | GCN250                     | 250 m                            | Based on ESA CCI land cover and HYSOGs250m soil hydrologic groups     |
| **River Network**     | MERIT Hydro                           | ~90 m (3 arc-sec)               | Flow direction, drainage area, and river width; threshold ≥1 km/50 km² |
| **Offshore Water Level (WL)** | HYCOM + TPXO 8                                  | HYCOM: ~8 km (3-hourly)<br>TPXO: ~24–28 km (10 min) | Linearly added tidal + non-tidal components (at 3-hourly intervals in the final output) for 29 ocean boundary points |
| **Wind & Pressure (WP)** | NHC-JTWC best track + WES (Deltares)               | Spiderweb grid | Holland model; 6-hourly fields extended 1000 km from storm center     |
| **Discharge (Dis)**   | USGS Stream Gauges                                   | 15-minute intervals              | 8 inflow points matched to rivers      |
| **Precipitation (Precip)** | ERA5 (ECMWF Reanalysis)                            | ~31 km, hourly                   | Downscaled hourly rainfall globally                                   |

---

## 🌐 Data Sources

All datasets used in this study are publicly available from the following sources:

- **CUDEM 1/9″ and 1/3″ Digital Elevation Models (NOAA)**  
  https://coast.noaa.gov/htdata/raster2/elevation/

- **HYCOM Sea Surface Height (SSH) Data**  
  https://www.hycom.org/hycom/overview

- **TPXO Global Tidal Model**  
  https://www.tpxo.net/global

- **National Land Cover Dataset (NLCD – 2016)**  
  https://www.mrlc.gov/data/nlcd-2016-land-cover-conus

- **Global Curve Number Dataset (GCN250)**  
  https://www.ndbc.noaa.gov/

- **ERA5 Reanalysis Data (Precipitation)**  
  https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download

- **Best Storm Track – NHC & JTWC**  
  - NHC: https://www.nhc.noaa.gov/data/  
  - JTWC: https://www.metoc.navy.mil/jtwc/jtwc.html?western-pacific

- **NOAA Water Level Observations (Tide Gauges)**  
  https://tidesandcurrents.noaa.gov/stations.html?type=Water+Levels

- **USGS Streamflow Gauges (for discharge)**  
  https://waterdata.usgs.gov/nwis

---

## 🧾 Code Overview

| File Name                               | Description |
|----------------------------------------|-------------|
| **HydroMT_Preprocessing_for_SFINCS.ipynb** | Python notebook using HydroMT plugin to prepare gridded input files (e.g., DEM, land cover, roughness, subgrid, etc.) required for the SFINCS model setup. |
| **Precipitation_Code.ipynb**  | Python notebook for processing ERA5 precipitation data. Extracts and formats rainfall inputs to be used in the SFINCS rainfall forcing file. |
| **1_Preprocessing_HYCOM_Data.m** | Extracts and processes HYCOM non-tidal sea surface height data for selected offshore points. Outputs time series for use as offshore boundary conditions in SFINCS. |
| **2_HYCOM_and_TPXO.m**           | Combines non-tidal water levels from HYCOM with tidal constituents from TPXO. Produces total offshore water level boundary time series for 29 locations. |
| **3_Nonlinearity_maps.m**        | Computes and visualizes spatial nonlinearity maps using SFINCS output. Calculates nonlinear interaction effects based on marginal water contribution. |
| **4_Nonlinearity_timeseries.m**  | Extracts and plots time series of nonlinearity at NOAA, USGS, and selected locations for all-driver simulations. |
| **5_Shapley_Values_maps.m**      | Calculates Shapley values for all 15 driver combinations, then normalizes and maps the spatial distribution of each driver’s contribution. |
| **6_Shapley_values_timeseries.m**| Computes and plots temporal evolution of Shapley values for all four flood drivers at 21 observation locations. |

---

## 🔁 Workflow Summary

### 1. Preprocessing

- Use `HydroMT_Preprocessing_for_SFINCS.ipynb` to generate SFINCS input rasters (DEM, land cover, Manning's n, infiltration, river mask, subgrid) using HydroMT plugins.
- Process ERA5 precipitation data using `Precipitation_Code.ipynb` to create rainfall input time series for the SFINCS model.
- Extract non-tidal sea surface height data from HYCOM using `1_Preprocessing_HYCOM_Data.m`.
- Combine HYCOM (non-tidal) and TPXO (tidal) datasets using `2_HYCOM_and_TPXO.m` to create complete offshore water level boundary inputs for 29 locations.

### 2. Compound Flood Simulation

- Run SFINCS v2.1.1 for all 15 combinations of flood drivers: Offshore Water Level (WL), Wind-Pressure (WP), Discharge (Dis), and Precipitation (Precip).
- Use marginal water contribution output (from water level and water depth) to capture nonlinear flood responses.

### 3. Nonlinearity Analysis

- Compute spatial maps of nonlinearity using `3_Nonlinearity_maps.m` by comparing nonlinear vs. linear water level superpositions.
- Extract and visualize time series of nonlinearity at NOAA, USGS, and selected locations using `4_Nonlinearity_timeseries.m`.

### 4. Shapley Value-Based Contribution

- Calculate Shapley values from all possible driver combinations using `5_Shapley_Values_maps.m`.
- Normalize and map spatial distribution of each flood driver's relative contribution.
- Analyze temporal variation in driver dominance using `6_Shapley_values_timeseries.m`.

---

## ▶️ Example Usage

### HydroMT Preprocessing (Python)

```python
#HydroMT_Preprocessing_for_SFINCS.ipynb

from hydromt_sfincs import SfincsModel

model = SfincsModel(root="sfincs_model")

# Specify an input dictionary with the grid settings x0,y0,dx,dy,nmax,mmax,rotation and epsg code.
# create SFINCS model with regular grid and characteristics of the input dictionary:
sf.setup_grid(
    x0=294515.3984,
    y0=3205240.0993,
    dx=200.0,
    dy=200.0,
    nmax=600,
    mmax=600,
    rotation=35,
    epsg=26915, #HydroMT is very sensitive to epsg code. Please use the right epsg code for your study area
)

# sf.region.boundary.plot(figsize=(6,6))
_ = sf.plot_basemap(plot_region=True, bmap="sat", figsize=(5,20), zoomlevel=10)

```

### Preprocessing Precipitation (Python)

```python
#Precipitation_Code.ipynb

root = "era5_hourly_data_beryl.nc"

# Open the dataset
dataset = xr.open_dataset(root)
dataset

lat = dataset['latitude'].values
lon = dataset['longitude'].values
precip = dataset['tp'].values
time = dataset['valid_time'].values

# Create a DataArray with the precipitation data
n1_zsmax = xr.DataArray(
    precip, 
    name="precip", 
    dims=("time", "latitude", "longitude"), 
    coords={"latitude": lat, "longitude": lon, "time": time}, 
    attrs={"units": "m", "long_name": "Total Precipitation", "_FillValue": "NaN", "coordinates": "spatial_ref"}
)

```

### Preprocessing HYCOM Data (MATLAB)

```matlab
% 1_Preprocessing_HYCOM_Data.m

fpath = ('ssh_2024_HYCOM.nc4');
ncdisp(fpath);

time = ncread(fpath, 'time');
lat = ncread(fpath, 'lat');
lon = ncread(fpath, 'lon');
surf_el = ncread(fpath, 'surf_el');

boundary_coords = load('sfincs_bnd.txt');
utmZone = repmat('16 N', size(boundary_coords, 1), 1);
[bnd_lat, bnd_lon] = utm2deg(boundary_coords(:, 1), boundary_coords(:, 2), utmZone);
bnd_lon_east = 360 + bnd_lon;

output_wl = zeros(length(time), size(boundary_coords, 1));
for t_idx = 1:length(time)
    current_wl_data = surf_el(:, :, t_idx);
    for b_idx = 1:size(boundary_coords, 1)
        output_wl(t_idx, b_idx) = griddata(lon, lat, current_wl_data', bnd_lon_east(b_idx), bnd_lat(b_idx), 'linear');
    end
end

fileID = fopen('sfincs_wl_output.txt', 'w');
for t_idx = 1:length(time)
    fprintf(fileID, '%.1f', (t_idx-1)*3*3600); % assuming 3-hour steps
    fprintf(fileID, ' %.4f', output_wl(t_idx, :));
    fprintf(fileID, '\n');
end
fclose(fileID);
```

### Combine HYCOM and TPXO Water Levels (MATLAB)

```matlab
% 2_HYCOM_and_TPXO.m

wl_data = load('sfincs_wl_output.txt');
wl_time = wl_data(:, 1);
wl_values = wl_data(:, 2:end);
wl_values(isnan(wl_values)) = 0;

bzs_data = load('sfincs_bzs.txt');
bzs_time = bzs_data(:, 1);
bzs_values = bzs_data(:, 2:end) + 0.3;

combined_wl = zeros(size(wl_values));
for i = 1:length(wl_time)
    idx = find(bzs_time == wl_time(i));
    if ~isempty(idx)
        combined_wl(i, :) = wl_values(i, :) + bzs_values(idx, :);
    end
end

fileID = fopen('combined_water_levels.txt', 'w');
for i = 1:length(wl_time)
    fprintf(fileID, '%.1f', wl_time(i));
    fprintf(fileID, ' %.4f', combined_wl(i, :));
    fprintf(fileID, '\n');
end
fclose(fileID);
```

### Nonlinearity Map Generation (MATLAB)

```matlab
% 3_Nonlinearity_maps.m

file_combinations = {
    'sfincs_map_WL+Storm.nc', {'sfincs_map_WL_Only.nc', 'sfincs_map_Storm_Only.nc'};
    ...
};

zb = ncread('sfincs_map_WL_Only.nc', 'zb');
zsmax = ncread(file_combinations{1, 1}, 'zsmax');
hmax = ncread(file_combinations{1, 1}, 'hmax');
h_combined = max(hmax, [], 3);
h_combined(zb < 0) = max(zsmax, [], 3)(zb < 0);
```

### Nonlinearity Time Series (MATLAB)

```matlab
% 4_Nonlinearity_timeseries.m

x = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'x');
y = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'y');
zb = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'zb');
time_data = ncread('sfincs_map_WL+Storm+Dis+Precipitation.nc', 'time');
actual_time = datetime(2024, 7, 1) + seconds(time_data);

noaa_coords = [...]; % Coordinates provided in code
nonlinearity_values = zeros(numel(noaa_coords)/2, length(actual_time));

% Each time step
for t = 1:length(actual_time)
    h_WL = read_depth('sfincs_map_WL_Only.nc', t, zb, 0.3);
    ...
    h_nonlinearity = safe_subtract(h_nonlinear, h_WL + h_Dis + h_Storm + h_Precipitation);
end

```

### Shapley Value Map (MATLAB)

```matlab
% 5_Shapley_Values_maps.m

h_WL = read_depth('sfincs_map_WL_Only.nc');
h_WindPressure = read_depth('sfincs_map_Storm_Only.nc');
...
phi_WL = weights(1)*safe_subtract(h_WL, 0) + weights(2)*... + weights(3)*... + weights(4)*...;

Norm_phi_WL = phi_WL ./ (phi_WL + phi_WindPressure + phi_Dis + phi_Precipitation);
pcolor(x, y, Norm_phi_WL); shading interp; caxis([-1, 1]);

```
### Shapley Value Time Series (MATLAB)

```matlab
% 6_Shapley_values_timeseries.m

utm_coords = [...]; % 21 locations
phi_values = zeros(21, num_time_steps, 5);
for t = 1:num_time_steps
    h_WL = read_depth('sfincs_map_WL_Only.nc', t, zb, hmin);
    ...
    phi_WL = weights(1)*... + weights(2)*... + weights(3)*... + weights(4)*...;
end

plot(actual_time, squeeze(phi_values(1, :, 1))); % WL at location 1

```
---

## 💻 Requirements

| Package       | Version |
|:--------------|:--------|
| hydromt       | >=0.7   |
| hydromt-sfincs| >=0.4   |
| xarray        | >=2023.1 |
| numpy         | >=1.21  |
| pandas        | >=1.3   |
| netCDF4       | >=1.5   |
| matplotlib    | >=3.5   |
| rasterio      | >=1.2   |
| geopandas     | >=0.10  |
| pyproj        | >=3.2   |

---

### 🔧 Installation

Install all dependencies using:

```bash
pip install hydromt hydromt-sfincs xarray numpy pandas netCDF4 matplotlib rasterio geopandas pyproj
```

---

## 🚀 Future Enhancements

- Integrate wave interactions to better capture their influence on compound flooding dynamics.
- Incorporate groundwater processes to account for subsurface contributions to flooding.
- Include urban drainage networks to improve representation of pluvial flooding in developed areas.
- Expand the framework to analyze multiple hurricane events under varying hydrodynamic conditions.
- Generalize findings by applying the methodology to different geographic regions and flood scenarios.

---

## Author

**MD Enayet Chowdhury**  
Email: enayet108@utexas.edu

Distributed under the Creative Commons CC0 1.0 Universal License. See LICENSE for more information.
