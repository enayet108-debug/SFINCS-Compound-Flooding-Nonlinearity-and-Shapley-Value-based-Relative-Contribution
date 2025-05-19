# Compound Flooding using SFINCS: Spatiotemporal Analysis of Nonlinearity and Shapley Value-based Relative Contribution
This repository provides the complete workflow and analysis codes for quantifying *spatiotemporal nonlinearity* and *Shapley-value based relative contribution* across space and time of four compound flood drivers: 
- Offshore Water Level (WL)
- Wind-Pressure (WP)
- River Discharge (Dis)
- Precipitation (Precip). 

The case study is based on *Hurricane Beryl (2024)* in the *Houston–Galveston region, Texas, USA*, using numerical simulations with the **Super-Fast INundation of CoastS (SFINCS)** model.

---
## 📁 Project Structure

```text
SFINCS-Compound-Flooding-Nonlinearity-and-Shapley-Value-based-Relative-Contribution/
├── notebooks/                            # Python notebooks for preprocessing
│   ├── 1_HydroMT_Preprocessing_for_SFINCS.ipynb
│   └── 2_Precipitation_Code.ipynb
│
├── matlab/                               # MATLAB scripts for analysis
│   ├── 3_Preprocessing_HYCOM_Data.m
│   ├── 4_HYCOM_and_TPXO.m
│   ├── 5_Nonlinearity_maps.m
│   ├── 6_Nonlinearity_timeseries.m
│   ├── 7_Shapley_Values_maps.m
│   └── 8_Shapley_Values_timeseries.m
│
├── input/                                # Input files for SFINCS model
│   ├── dem3.tif
│   ├── crm.tif
│   ├── tx_man400_utm.tif
│   ├── GCN250_ARCI.tif
│   ├── GCN250_ARCII.tif
│   ├── GCN250_ARCIII.tif
│   ├── flwdir.tif
│   ├── rivwth.tif
│   ├── uparea.tif
│   ├── sfincs_grid_region.shp
│   ├── ssh_2024_2.nc4
│   ├── Beryl.spw
│   ├── Beryl.pol
│   ├── sfincs.bca
│   ├── sfincs.bnd
│   ├── sfincs.bzs
│   ├── sfincs.dis
│   ├── sfincs.ind
│   ├── sfincs.inp
│   ├── sfincs.msk
│   ├── sfincs.obs
│   ├── sfincs.scs
│   ├── sfincs.src
│   ├── sfincs_subgrid.nc
│   ├── sfincs_map.nc
│   ├── sfincs_his.nc
│   ├── era5_hourly_data_beryl.nc
│   └── Precip_2d.nc
│
├── LICENSE                               # Open-source license
└── README.md                             # Project overview

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

<sub>
Abbreviations:  
DEM = Digital Elevation Model,  
CUDEM = Continuously Updated Digital Elevation Model,  
CRM = Coastal Relief Model,  
NLCD = National Land Cover Database,  
GCN = Global Curve Number,  
CCI = Climate Change Initiative,  
HYSOGs = Global Hydrologic Soil Groups (HYSOGs250m),  
MERIT = Multi-Error-Removed Improved-Terrain Hydro,  
HYCOM = HYbrid Coordinate Ocean Model,  
TPXO = Tidal Prediction using TOPEX/POSEIDON,  
NHC = National Hurricane Center,  
JTWC = Joint Typhoon Warning Center,  
WES = Wind Enhance Scheme,  
USGS = United States Geological Survey,  
ERA5 = ECMWF Reanalysis 5th Generation,  
ECMWF = European Centre for Medium-Range Weather Forecasts.  
</sub>

---
## 📂 Input Files Overview

Example input files can be downloaded from:
[SFINCS Spatiotemporal Analysis Zenodo](https://zenodo.org/doi/10.5281/zenodo.15460823)

| Related Item              | File Name                        | File Type | Resolution                 | Description |
|---------------------------|----------------------------------|-----------|----------------------------|-------------|
| Digital Elevation Model   | dem3.tif                         | TIFF      | 3 m (coast), 10 m (offshore) | A 32-bit floating point TIFF raster based on CUDEM. The uncompressed file size is 13.45 GB, contains no color map with basic mensuration capabilities, and 8 pyramid levels for faster rendering. |
| Digital Elevation Model   | crm.tif                          | TIFF      | 30 m                       | A 32-bit floating point TIFF raster from the CRM bathymetry datasets. The uncompressed size is 2.25 GB, with no colormap, and includes 5 pyramid levels for efficient rendering and basic mensuration capabilities. |
| Manning’s Roughness       | tx_man400_utm.tif                | TIFF      | 30 m                       | A TIFF raster file composed of 64-bit double-precision values. The file has an uncompressed size of 5.67 MB, contains no colormap or compression, uses -9999 as the NoData value, supports basic mensuration capabilities, and includes arbitrary pyramid levels for improved rendering performance. |
| Infiltration              | GCN250_ARCI/ARCII/ARCIII.tif     | TIFF      | 250 m                      | A TIFF raster file composed of 8-bit unsigned integer values. The file has an uncompressed size of 10.81 GB, has no colormap or pyramids, designates 255 as the NoData value, and supports basic mensuration capabilities. |
| Flow Direction            | flwdir                           | TIFF      | ~90 m (3 arc-sec)          | A MERIT Hydro derived flow direction raster from the latest elevation data (MERIT DEM) and water body datasets (G1WBM, GSWO, and OpenStreetMap). Contains a single 8-bit unsigned band, and no colormap, no compression, a NoData value of 255, and includes 4 pyramid levels with nearest neighbor resampling for enhanced display performance, along with basic mensuration capabilities. |
| River Width               | rivwth                           | TIFF      | ~90 m (3 arc-sec)          | Contains a single 32-bit floating point band. The file is uncompressed with a size of 549.39 MB, has no colormap, includes 4 pyramid levels with nearest neighbor resampling, and supports basic mensuration capabilities. |
| Upstream Drainage Area    | uparea                           | TIFF      | ~90 m (3 arc-sec)          | A single 32-bit floating point band with an uncompressed size of 549.39 MB, no compression or colormap, and includes 4 pyramid levels with nearest neighbor resampling for improved visualization, along with basic mensuration capabilities. |
| Study Area                | sfincs_grid_region.shp           | SHP       | N/A                        | A shapefile representing a rectangular study area encompassing the Houston–Galveston region of Texas, USA. |
| Offshore Water Level      | ssh_2024_2.nc4                   | netCDF4   | ~8 km (3-hourly)           | Downloaded HYCOM sea surface elevation data extracted for the coordinates specified in the ‘sfincs.bnd’ file at a 3-hourly interval. |
| Wind-Pressure (Hurricane) | Beryl.spw                        | SPW       | 6-hourly                   | An NHC-JTWC best track data for Hurricane Beryl as a spiderweb file, including wind speed, direction, and pressure. Coordinates: meters in projected UTM zone, data: m/s, wind_from_direction in degrees, p_drop in Pa. |
| Wind-Pressure (Hurricane) | Beryl.pol                        | POL       | 6-hourly                   | A polygon-based file representing the Hurricane’s best track. |
| SFINCS Input              | sfincs.bca                       | BCA       | N/A                        | A file carrying the astronomical tide constituents’ amplitude and phase information. |
| SFINCS Input              | sfincs.bnd                       | BND       | N/A                        | To specify water-level time-series to the boundary cells (msk=2). The input locations are specified here. Units: meters in projected UTM zone. |
| SFINCS Input              | sfincs.bzs                       | BZS       | 3-hourly                   | The (slowly varying) water level time series specified per input location. Units: m above reference level. This data is achieved after combining the HYCOM and TPXO data. |
| SFINCS Input              | sfincs.dis                       | DIS       | 15-minutes                 | The discharge time series specified per input location. Unit: cubic meters per second. |
| SFINCS Input              | sfincs.ind                       | IND       | N/A                        | Describes the indices of active grid cells within the overall grid. Not used by SFINCS with ASCII input. |
| SFINCS Input              | sfincs.inp                       | INP       | N/A                        | General input file of SFINCS describing all model settings, the domain, forcing, and structures. |
| SFINCS Input              | sfincs.msk                       | MSK       | N/A                        | Indicates for every cell whether it is an inactive cell (msk=0), active cell (msk=1), WL boundary cell (msk=2), or outflow boundary cell (msk=3). |
| SFINCS Input              | sfincs.obs                       | OBS       | N/A                        | Observation points for output time series. Units: meters in projected UTM zone. |
| SFINCS Input              | sfincs.scs                       | SCS       | N/A                        | Grid-based input using Curve Number method A (without recovery), same structure as `depfile`, binary input. |
| SFINCS Input              | sfincs.src                       | SRC       | N/A                        | Discharge point coordinates in UTM meters. |
| SFINCS Input              | sfincs_subgrid.nc                | netCDF    | N/A                        | A netCDF file including all the subgrid information, to run SFINCS. |
| SFINCS Output              | sfincs_map.nc                    | netCDF    | 10 m, 1-hourly             | SFINCS output file carrying the study area’s spatial visualization data generated by SFINCS. It can be opened by platforms like Quickplot, Panoply, Matlab, Python, etc. Included variables: spatial coordinates, bed level elevation, instantaneous water level, instantaneous water depth, maximum water level and water depth across all timesteps, etc. |
| SFINCS Output              | sfincs_his.nc                    | netCDF    | 6-minutes                   | This contains point-based time series data for specified coordinates in the `sfincs.obs` file. Included variables: spatial coordinates, bed level elevation of observation points, time step of output (it can be different from the `sfincs_map.nc` file), instantaneous water level, instantaneous water depth, etc. |
| Precipitation             | era5_hourly_data_beryl.nc        | netCDF    | ~31 km, hourly             | Directly downloaded hourly single-level total precipitation data from the ERA5 database. |
| Precipitation (Processed)             | Precip_2d.nc                     | netCDF    | ~31 km, hourly             | Processed netCDF file originally downloaded from ERA5 data, to be used inside the SFINCS environment. |

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
| **1_HydroMT_Preprocessing_for_SFINCS.ipynb** | Python notebook using HydroMT plugin to prepare gridded input files (e.g., DEM, land cover, roughness, subgrid, etc.) required for the SFINCS model setup. |
| **2_Precipitation_Code.ipynb**  | Python notebook for processing ERA5 precipitation data. Extracts and formats rainfall inputs to be used in the SFINCS rainfall forcing file. |
| **3_Preprocessing_HYCOM_Data.m** | Extracts and processes HYCOM non-tidal sea surface height data for selected offshore points. Outputs time series for use as offshore boundary conditions in SFINCS. |
| **4_HYCOM_and_TPXO.m**           | Combines non-tidal water levels from HYCOM with tidal constituents from TPXO. Produces total offshore water level boundary time series for 29 locations. |
| **5_Nonlinearity_maps.m**        | Computes and visualizes spatial nonlinearity maps using SFINCS output. Calculates nonlinear interaction effects based on marginal water contribution. |
| **6_Nonlinearity_timeseries.m**  | Extracts and plots time series of nonlinearity at NOAA, USGS, and selected locations for all-driver simulations. |
| **7_Shapley_Values_maps.m**      | Calculates Shapley values for all 15 driver combinations, then normalizes and maps the spatial distribution of each driver’s contribution. |
| **8_Shapley_Values_timeseries.m**| Computes and plots temporal evolution of Shapley values for all four flood drivers at 21 observation locations. |

---

## 🔁 Workflow Summary

### 1. Preprocessing

- Use `1_HydroMT_Preprocessing_for_SFINCS.ipynb` to generate SFINCS input rasters (DEM, land cover, Manning's n, infiltration, river mask, subgrid) using HydroMT plugins.
- Process ERA5 precipitation data using `2_Precipitation_Code.ipynb` to create rainfall input time series for the SFINCS model.
- Extract non-tidal sea surface height data from HYCOM using `3_Preprocessing_HYCOM_Data.m`.
- Combine HYCOM (non-tidal) and TPXO (tidal) datasets using `4_HYCOM_and_TPXO.m` to create complete offshore water level boundary inputs for 29 locations.

### 2. Compound Flood Simulation

- Run SFINCS v2.1.1 for all 15 combinations of flood drivers: Offshore Water Level (WL), Wind-Pressure (WP), Discharge (Dis), and Precipitation (Precip).
- Use marginal water contribution output (from water level and water depth) to capture nonlinear flood responses.

### 3. Nonlinearity Analysis

- Compute spatial maps of nonlinearity using `5_Nonlinearity_maps.m` by comparing nonlinear vs. linear water level superpositions.
- Extract and visualize time series of nonlinearity at NOAA, USGS, and selected locations using `6_Nonlinearity_timeseries.m`.

### 4. Shapley Value-Based Contribution

- Calculate Shapley values from all possible driver combinations using `7_Shapley_Values_maps.m`.
- Normalize and map spatial distribution of each flood driver's relative contribution.
- Analyze temporal variation in driver dominance using `8_Shapley_Values_timeseries.m`.

---

## ▶️ Example Usage

### HydroMT Preprocessing (Python)

```python
#1_HydroMT_Preprocessing_for_SFINCS.ipynb

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
#2_Precipitation_Code.ipynb

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
% 3_Preprocessing_HYCOM_Data.m

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
% 4_HYCOM_and_TPXO.m

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
% 5_Nonlinearity_maps.m

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
% 6_Nonlinearity_timeseries.m

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
% 7_Shapley_Values_maps.m

h_WL = read_depth('sfincs_map_WL_Only.nc');
h_WindPressure = read_depth('sfincs_map_Storm_Only.nc');
...
phi_WL = weights(1)*safe_subtract(h_WL, 0) + weights(2)*... + weights(3)*... + weights(4)*...;

Norm_phi_WL = phi_WL ./ (phi_WL + phi_WindPressure + phi_Dis + phi_Precipitation);
pcolor(x, y, Norm_phi_WL); shading interp; caxis([-1, 1]);

```
### Shapley Value Time Series (MATLAB)

```matlab
% 8_Shapley_Values_timeseries.m

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
