import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

print('Starting CheckFates_Points')

# File names. Modify manually!
case_name = 'i1850.FATES-NOCOMP.ne30pg3_tn14.alpha08d.20250218_coldstart_adrianas_tuning'
#case_name = 'n1850.FATES-NOCOMP-postAD.ne30pg3_tn14.alpha08d.20250206_fixFincl1'
case_dir = f'/cluster/work/users/kjetisaa/archive/{case_name}/lnd/hist/'

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/{case_name}.clm2.h0.*-*.nc'))
#variables = ["FATES_NCOHORTS", "FATES_NPATCHES", "TLAI","TOTECOSYSC", "FATES_VEGC", "FATES_NONSTRUCTC", "FATES_STRUCTC", "FATES_STOREC", "FATES_SEED_BANK", "FATES_REPROC", "FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"]
#variables = ["FATES_NCOHORTS", "FATES_NPATCHES", "TLAI","TOTECOSYSC", "FATES_NONSTRUCTC", "FATES_STRUCTC", "FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"]
#variables = ["FATES_NPLANT_SZ","FATES_NCOHORTS", "FATES_NPATCHES", "TLAI","TOTECOSYSC", "FATES_NONSTRUCTC", "FATES_STRUCTC", "FATES_MORTALITY_AGESCEN_SZ", "FATES_MORTALITY_BACKGROUND_SZ", "FATES_MORTALITY_CSTARV_SZ", "FATES_MORTALITY_FIRE_SZ", "FATES_MORTALITY_FREEZING_SZ", "FATES_MORTALITY_HYDRAULIC_SZ", "FATES_MORTALITY_IMPACT_SZ", "FATES_MORTALITY_LOGGING_SZ", "FATES_MORTALITY_SENESCENCE_SZ"]
variables = ["FATES_NPLANT_SZ","FATES_NCOHORTS", "FATES_NPATCHES", "TLAI","TOTECOSYSC", 
"FATES_NONSTRUCTC", "FATES_STRUCTC","FATES_MORTALITY_BACKGROUND_SZ", "FATES_MORTALITY_CSTARV_SZ", 
"FATES_MORTALITY_FREEZING_SZ", "FATES_MORTALITY_HYDRAULIC_SZ", "FATES_MORTALITY_IMPACT_SZ"]

# Locations
locations = {
    'Amazon': {'lat': -0.5, 'lon': -65},
    'Temperate': {'lat': 47.5, 'lon': 1.8},
    'Boreal': {'lat': 57, 'lon': -113},
    'Siberian': {'lat': 67.5, 'lon': 100.0}
}

#convert negative longitudes to positive
for loc in locations:
    if locations[loc]['lon'] < 0:
        locations[loc]['lon'] = 360 + locations[loc]['lon']

# Initialize a dictionary to store the accumulated results for each location
results = {loc: {var: [] for var in variables} for loc in locations}
time = {loc: [] for loc in locations}

# Loop through all timeseries files except the last one
first_file = True
if len(timeseries_files) > 1:
    files_to_read = timeseries_files[:-1]
else:
    files_to_read = timeseries_files

for filename in files_to_read: 
    print(filename)   
    try:        
        with xr.open_dataset(filename, engine='netcdf4') as data:
            if first_file:        
                print('Reading latitude and longitude')
                lats = data['lat'].values
                lons = data['lon'].values
                first_file = False

                for loc, coords in locations.items():
                    lat = coords['lat']
                    lon = coords['lon']
                    # Find the absolute difference between the given lat/lon and the lat/lon values in the dataset
                    lat_diff = np.abs(lats - lat)
                    lon_diff = np.abs(lons - lon)
                    # Find the index of the minimum difference for lat and lon, ignoring NaNs
                    closest_idx = np.nanargmin(lat_diff + lon_diff)
                    locations[loc]['closest_idx'] = closest_idx
                    print(f"{loc}: The index of the closest point is: {closest_idx}")
                    print(f"{loc}: Closest point is at: {lats[closest_idx]}, {lons[closest_idx]}")
                    if 'landfrac' in data and 'PCT_LANDUNIT' in data:
                        landfrac = data['landfrac'][closest_idx].values
                        pct_landunit = data['PCT_LANDUNIT'][:,:, closest_idx].values
                        print(f"{loc}: Land fraction (landfrac) is: {landfrac}")
                        print(f"{loc}: Percentage of each landunit (PCT_LANDUNIT) is: {pct_landunit}")

            for loc in locations:
                closest_idx = locations[loc]['closest_idx']
                for var in variables:
                    if var in data:
                        var_data = data[var]
                        # Check if the variable has more than 2 dimensions
                        if var_data.ndim > 2:
                            # Sum over the second dimension                            
                            var_data = var_data.sum(axis=1)
                        var_data = var_data[:, closest_idx]
                        results[loc][var].append(var_data.values)
                    
                time[loc].append(data['time'].values)
    except FileNotFoundError:
        print(f"File not found: {filename}")
    except ValueError as e:
        print(f"Error reading the file: {e}")

# Convert lists to numpy arrays
for loc in locations:
    for var in variables:
        if results[loc][var]:  # Check if the list is not empty
            results[loc][var] = np.concatenate(results[loc][var])
        else:
            results[loc][var] = np.array([])

    if time[loc]:  # Check if the list is not empty
        time[loc] = np.concatenate(time[loc])
    else:
        time[loc] = np.array([])

plot_reduced_list = False
if plot_reduced_list:
    #Collect variables to total carbon pools
    fates_litter_vars = ["FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"]
    #fates_veg_vars = ["FATES_STRUCTC", "FATES_NONSTRUCTC", "FATES_STOREC", "FATES_REPROC"]
    for loc in locations:   
        if all(var in results[loc] for var in fates_litter_vars):
            results[loc]["FATES_LITTERC"] = results[loc]["FATES_LITTER_CWD_ELDC"] + results[loc]["FATES_LITTER_AG_CWD_EL"] + results[loc]["FATES_LITTER_AG_FINE_EL"] + results[loc]["FATES_LITTER_BG_CWD_EL"] + results[loc]["FATES_LITTER_BG_FINE_EL"]            
    #    if all(var in results[loc] for var in fates_veg_vars):
    #        results[loc]["FATES_TOTVEGC"] = results[loc]["FATES_STRUCTC"] + results[loc]["FATES_NONSTRUCTC"] + results[loc]["FATES_STOREC"] + results[loc]["FATES_REPROC"]  
    #plot_variables = [var for var in variables if var not in ["FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL", "FATES_LEAFC", "FATES_FROOTC", "FATES_STOREC", "FATES_SAPWOODC", "FATES_REPROC"]] + ["FATES_LITTERC", "FATES_TOTVEGC"]
    plot_variables = [var for var in variables if var not in ["FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"]] + ["FATES_LITTERC"]
else:
    plot_variables = variables

# Create the suptitle with location names and coordinates
location_titles = ', '.join([f"{loc} (lat: {coords['lat']}, lon: {coords['lon']})" for loc, coords in locations.items()])
fig_title = f"Timeseries Data at {location_titles}"

# Plotting
fig, axes = plt.subplots(len(plot_variables), len(locations), figsize=(20, 15), sharex='col')
fig.suptitle(fig_title)

# Print the specified variables if they exist
for loc in locations:
    print(f"\nLocation: {loc}")
    if all(var in results[loc] for var in ["FATES_VEGC", "FATES_TOTVEGC"]):
        print(f"FATES_VEGC: {results[loc]['FATES_VEGC']}")
        print(f"FATES_TOTVEGC: {results[loc]['FATES_TOTVEGC']}")
        print(f"FATES_VEGC-FATES_TOTVEGC: {results[loc]['FATES_VEGC'] - results[loc]['FATES_TOTVEGC']}")

for i, var in enumerate(plot_variables):
    for j, loc in enumerate(locations):
        axes[i, j].plot(results[loc][var])
        if var in ["FATES_LITTERC", "FATES_TOTVEGC"]:
            axes[i, j].set_title(f"{var} (kg C m-2)")
        else:   
            axes[i, j].set_title(f"{var} ({data[var].units})")
        axes[i, j].grid(True)
        if i == len(plot_variables) - 1:
            axes[i, j].set_xlabel('Time')
        if j == 0:
            axes[i, j].set_ylabel('')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(f"figs/Point_Timeseries_{case_name}.png")

print('Finished CheckFates_Points')