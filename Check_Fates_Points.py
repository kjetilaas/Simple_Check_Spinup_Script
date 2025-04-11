import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

print('Starting CheckFates_Points')

# File names. Modify manually!
#case_name = 'i1850.FATES-NOCOMP-LUH2-evolve.ne30pg3_tn14.alpha09a.20250321'
case_name = 'i1850.FATES-NOCOMP-noLU_cal15_twosteam.ne30pg3_tn14.noresm3_0_alpha01.20250326'
#case_name = 'n1850.ne30_tn14.hybrid_fates-nocomp.3.0a02.20250408'
case_dir = f'/cluster/work/users/kjetisaa/archive/{case_name}/lnd/hist/'
#case_dir = f'/cluster/work/users/tomast/archive/{case_name}/lnd/hist/'
#case_dir = f'/cluster/work/users/kjetisaa/noresm/{case_name}/run/'
#case_dir = f'/cluster/work/users/tomast/noresm/{case_name}/run/'

obs_dir = f'/cluster/home/kjetisaa/OBS_ILAMB/lai/'


# Option to process only the last 10 years, only works for monthly data
process_last_10_years = True

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/{case_name}.clm2.h0.*-*.nc'))

calc_annual = False
if len(timeseries_files) > 120 and not process_last_10_years:
    calc_annual = True
    print(f"Calculating annual means for {len(timeseries_files)} files.")

if process_last_10_years:
    # Filter files to include only the last 10 years
    timeseries_files = timeseries_files[-120:]  # Assuming monthly data, 10 years = 120 months

plot_structure = False
if plot_structure:
    variables = ["FATES_NPLANT_SZ","FATES_NCOHORTS", "FATES_NPATCHES", "TLAI","TOTECOSYSC", 
    "FATES_NONSTRUCTC", "FATES_STRUCTC"]
else:
    variables = ["TLAI", "FATES_GPP", "FATES_NPP", "BTRAN", "SOILWATER_10CM", "TOTSOMC", "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_MORTALITY_BACKGROUND_SZ", 
    "FATES_MORTALITY_CSTARV_SZ", "FATES_MORTALITY_FREEZING_SZ", "FATES_MORTALITY_HYDRAULIC_SZ", 
    "FATES_MORTALITY_IMPACT_SZ", "FATES_MORTALITY_CANOPY_SZ", "FATES_MORTALITY_USTORY_SZ"] #"FATES_MORTALITY_FIRE_SZ"

# Locations
plot_region = 'Global'# 'Norway', 'Nordic', 'Global', 'Biased'
if plot_region == 'Global':
    locations = {
        'Amazon': {'lat': -0.5, 'lon': -65},  # Tropical rainforest
        'C4_Gr': {'lat': 13.0, 'lon': 16.5},  # C4 grassland in Africa
        'BL_CD': {'lat': 39.0, 'lon': -80.5},  # Broadleaf temperate forest in Eastern US
        'Boreal': {'lat': 57.5, 'lon': -121.5},  # Boreal forest
        'Cool_C3': {'lat': 53, 'lon': -7.5},  # Cool C3 grassland
        'Larch': {'lat': 61, 'lon': 122},  # Larch forest
        'Arc_grass': {'lat': 68, 'lon': 120.0},  # Siberian grassland
    }
elif plot_region == 'Biased':    
    locations = {        
        'C4_Gr': {'lat': 13.0, 'lon': 16.5},  # C4 grassland in Africa        
        'Boreal': {'lat': 57.5, 'lon': -121.5},  # Boreal forest
        'Cool_C3': {'lat': 53, 'lon': -7.5},  # Cool C3 grassland
        'Larch': {'lat': 61, 'lon': 122},  # Larch forest        
        'GUICHOU': {'lat': 27.0, 'lon': 106.0}  # Guichou, China
    }
elif plot_region == 'Norway':
    #Plot locations in Norway   
    locations = {
        'Oslo': {'lat': 60.0, 'lon': 10.0},  # Oslo
        'Tromsø': {'lat': 70.0, 'lon': 20.0},  # Tromsø
        'Bergen': {'lat': 60.5, 'lon': 5.5},  # Bergen
        'Stavanger': {'lat': 58.5, 'lon': 5.5},  # Stavanger
        'Kristiansand': {'lat': 58.0, 'lon': 8.0},  # Kristiansand
        'Trondheim': {'lat': 63.5, 'lon': 10.5},  # Trondheim
        'Kirkenes': {'lat': 69.0, 'lon': 29.0}  # Kirkenes
    }
elif plot_region == 'Nordic':
    locations = {
        'Hyytiälä': {'lat': 61.5, 'lon': 24.0},  # Hyytiälä, Finland
        'Hurdal': {'lat': 60.5, 'lon': 10.0},  # Hurdal, Norway
        'Abisko': {'lat': 68.4, 'lon': 18.8}  # Abisko, Sweden
    }

# PFT names
pft_names = [
    "broadleaf_evergreen_tropical_tree",
    "needleleaf_evergreen_extratrop_tree",
    "needleleaf_colddecid_extratrop_tree",
    "broadleaf_evergreen_extratrop_tree",
    "broadleaf_hydrodecid_tropical_tree",
    "broadleaf_colddecid_extratrop_tree",
    "broadleaf_evergreen_extratrop_shrub",
    "broadleaf_hydrodecid_extratrop_shrub",
    "broadleaf_colddecid_extratrop_shrub",
    "broadleaf_evergreen_arctic_shrub",
    "broadleaf_colddecid_arctic_shrub",
    "arctic_c3_grass",
    "cool_c3_grass",
    "c4_grass"
]

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
    files_to_read = timeseries_files#[:-1]
else:
    files_to_read = timeseries_files

for filename in files_to_read: 
    print(filename)   
    try:        
        with xr.open_dataset(filename, engine='netcdf4') as data:
            if first_file:        
                #print('Reading latitude and longitude')
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
#print(f"{loc}: The index of the closest point is: {closest_idx}")
                    #print(f"{loc}: Closest point is at: {lats[closest_idx]}, {lons[closest_idx]}")
                    if 'landfrac' in data:                        
                        landfrac = data['landfrac'][closest_idx].values
                        #print(f"{loc}: Land fraction (landfrac) is: {landfrac}")
                    if 'PCT_LANDUNIT' in data:
                        pct_landunit = data['PCT_LANDUNIT'][:,:, closest_idx].values
                        #print(f"{loc}: Percentage of each landunit (PCT_LANDUNIT) is: {pct_landunit}")
                    if 'FATES_NOCOMP_PATCHAREA_PF' in data:
                        pct_pft= data['FATES_NOCOMP_PATCHAREA_PF'][:,:, closest_idx].values
                        pct_pft_str = ', '.join([f"{p:.2f}" for p in pct_pft.flatten()])
                        #print(f"{loc}: Percentage PFT is: {pct_pft_str}")
                        #print(f'{loc}: Dominant PFT is: {pft_names[np.argmax(pct_pft[-1])]}, {np.max(pct_pft[-1]*100):.2f} %')
                        for i, pft in enumerate(pft_names):
                            if pct_pft[-1,i] > 0.1:
                                print(f'{loc}: {pft}: {pct_pft[-1,i]*100:.2f} %')


            for loc in locations:
                closest_idx = locations[loc]['closest_idx']
                for var in variables:
                    if var in data:  # Check if the variable exists in the dataset
                        var_data = data[var]
                        # Check if the variable has more than 2 dimensions
                        if var_data.ndim > 2:
                            # Sum over the second dimension                            
                            var_data = var_data.sum(axis=1)
                        var_data = var_data[:, closest_idx]
                        results[loc][var].append(var_data.values)
                    #else:
                        #print(f"Warning: Variable '{var}' not found in file {filename}. Skipping.")
                    
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

# Add CUE if NPP and GPP are available
if all(var in results[loc] for var in ["FATES_NPP", "FATES_GPP"]):
    for loc in locations:
        # Check for division by zero
        results[loc]["CUE"] = np.divide(results[loc]["FATES_NPP"], results[loc]["FATES_GPP"], out=np.zeros_like(results[loc]["FATES_NPP"]), where=results[loc]["FATES_GPP"] != 0)
        # Set negative values to zero
        results[loc]["CUE"] = np.clip(results[loc]["CUE"], 0, None)
    # Append "CUE" to plot_variables only once
    if "CUE" not in plot_variables:
        plot_variables.append("CUE")

if calc_annual:
    # Calculate annual means
    for loc in locations:
        for var in plot_variables:
            if var in results[loc]:
                results[loc][var] = np.mean(results[loc][var].reshape(-1, 12), axis=1)
                time[loc] = time[loc][::12]  # Adjust time to match the annual data

# Create the suptitle with location names and coordinates
location_titles = ', '.join([f"{loc} (lat: {coords['lat']}, lon: {coords['lon']})" for loc, coords in locations.items()])
if process_last_10_years:
    fig_title = f"Timeseries Data at {location_titles} (last 10 years)"
else:
    fig_title = f"Timeseries Data at {location_titles}"

# Plotting
fig, axes = plt.subplots(len(plot_variables), len(locations), figsize=(30, 25), sharex='col')
fig.suptitle(fig_title, fontsize=20)

# Print the specified variables if they exist
for loc in locations:
    print(f"\nLocation: {loc}")
    if all(var in results[loc] for var in ["FATES_VEGC", "FATES_TOTVEGC"]):
        print(f"FATES_VEGC: {results[loc]['FATES_VEGC']}")
        print(f"FATES_TOTVEGC: {results[loc]['FATES_TOTVEGC']}")
        print(f"FATES_VEGC-FATES_TOTVEGC: {results[loc]['FATES_VEGC'] - results[loc]['FATES_TOTVEGC']}")

# Define the observational datasets and their file paths
obs_datasets = {
    'AVH15C1': f'{obs_dir}AVH15C1/lai.nc',
    'AVHRR': f'{obs_dir}AVHRR/lai_0.5x0.5.nc',
    'GIMMS_LAI4g': f'{obs_dir}GIMMS_LAI4g/cao2023_lai.nc',
    'MODIS': f'{obs_dir}MODIS/lai_0.5x0.5.nc'
}

# Read and process all observational LAI datasets
obs_lai_results = {obs_name: {} for obs_name in obs_datasets.keys()}

for obs_name, obs_file in obs_datasets.items():
    try:
        with xr.open_dataset(obs_file, engine='netcdf4') as obs_data:
            obs_lai = obs_data['lai']  # Variable lai(time, lat, lon)
            obs_lats = obs_data['lat'].values
            obs_lons = obs_data['lon'].values

            # Convert longitudes to the range [0, 360] for consistency
            obs_lons = np.where(obs_lons < 0, obs_lons + 360, obs_lons)

            # Average LAI into a climatological year (monthly mean over all years)
            obs_lai_clim = obs_lai.groupby('time.month').mean(dim='time')

            # Extract LAI for each location
            for loc, coords in locations.items():
                lat = coords['lat']
                lon = coords['lon']

                # Find the closest grid cell in the observational data
                lat_diff = np.abs(obs_lats - lat)
                lon_diff = np.abs(obs_lons - lon)
                closest_idx = np.unravel_index(np.argmin(lat_diff[:, None] + lon_diff), (len(obs_lats), len(obs_lons)))

                # Extract the climatological LAI for the closest grid cell
                obs_lai_results[obs_name][loc] = obs_lai_clim[:, closest_idx[0], closest_idx[1]].values
                #print(f"{loc}: {obs_name} LAI extracted for lat: {lat}, lon: {lon}, LAI: {obs_lai_results[obs_name][loc]}")

    except FileNotFoundError:
        print(f"Observational LAI file not found: {obs_file}")
    except Exception as e:
        print(f"Error processing {obs_name}: {e}")

# Update the plotting section to include all observational LAI datasets
for i, var in enumerate(plot_variables):
    for j, loc in enumerate(locations):
        if var in results[loc] and len(results[loc][var]) > 0:  # Check if the variable has data
            axes[i, j].plot(results[loc][var], label='Model')
            if var == 'TLAI' and not calc_annual:
                # Plot all observational LAI datasets
                for obs_name, obs_data in obs_lai_results.items():
                    if loc in obs_data:
                        # Repeat the 12 months of obs for every model year
                        obs_data[loc] = np.tile(obs_data[loc], int(len(results[loc][var]) / 12))
                        axes[i, j].plot(obs_data[loc], label=obs_name, linestyle='--')
            if var in ["FATES_LITTERC", "FATES_TOTVEGC"]:
                axes[i, j].set_title(f"{var} (kg C m-2)")
            elif var == "CUE":
                axes[i, j].set_title(f"{var} (unitless)")  # Handle CUE separately
            else:   
                axes[i, j].set_title(f"{var} ({data[var].units if var in data else 'unknown units'})")
            axes[i, j].grid(True)
        else:
            print(f"Warning: No data available for variable '{var}' at location '{loc}'. Skipping plot.")
        if i == len(plot_variables) - 1:
            axes[i, j].set_xlabel('Time')
        if j == 0:
            axes[i, j].set_ylabel('')
        if var == 'TLAI' and plot_region != 'Global':# and j == 0:
            # Add legend on right side of the plot
            axes[i, j].legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)            

plt.tight_layout(rect=[0, 0, 1, 0.96])
output_filename = f"figs/Point_Timeseries_{case_name}"
#if not global, add region to filename
if plot_region != 'Global': 
    output_filename += f"_{plot_region}"
if process_last_10_years:
    output_filename += "_last10years"
if plot_structure:
    output_filename += "_structure" 
else:
    output_filename += "_survival"
output_filename += ".png"
plt.savefig(output_filename)

print('Finished CheckFates_Points')