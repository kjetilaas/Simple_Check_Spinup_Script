import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define the directory containing the data files
#case_name = 'i2000.CRUJRA0.5_Scandinavia_nocomp'
#case_dir = f'/nird/datalake/NS9560K/kjetisaa/{case_name}/lnd/hist/'

case_name = 'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S2_TRENDY2025_pt3.202508021'
case_dir = f'/nird/datalake/NS9560K/kjetisaa/TRENDY25/{case_name}/lnd/hist/'

print(f"Data directory: {case_dir}")

# Define the file pattern for the NetCDF files
#file_pattern = f'{case_name}.clm2.h0.207'

# Define the region for Scandinavia (latitude and longitude bounds)
lat_min, lat_max = 50, 75
lon_min, lon_max = 3, 30

syear = 2015
eyear = 2024

# Define the variable to process
var = 'FATES_GPP'  # Change to 'FATES_LAI', 'FATES_NEP', or 'FATES_NPP'

# Define color maps and limits for each variable
var_settings = {
    'FATES_LAI': {'cmap': 'YlGn', 'clim_min': 0, 'clim_max': 5, 'is_flux': False},
    'FATES_NEP': {'cmap': 'RdYlGn', 'clim_min': -0.5, 'clim_max': 0.5, 'is_flux': True},
    'FATES_NPP': {'cmap': 'RdYlGn', 'clim_min': -1.5, 'clim_max': 1.5, 'is_flux': True},
    'FATES_GPP': {'cmap': 'RdYlGn', 'clim_min': -2, 'clim_max': 2, 'is_flux': True},
}

if var not in var_settings:
    raise ValueError(f"Variable {var} is not supported. Choose from {list(var_settings.keys())}.")

cmap = var_settings[var]['cmap']
clim_min = var_settings[var]['clim_min']
clim_max = var_settings[var]['clim_max']
is_flux = var_settings[var]['is_flux']


# Function to get the number of days in a month
def days_in_month(month):
    """Return the number of days in a month."""
    if month in [1, 3, 5, 7, 8, 10, 12]:
        return 31
    elif month in [4, 6, 9, 11]:
        return 30
    elif month == 2:
        return 28
    else:
        raise ValueError("Invalid month: {}".format(month)) 


# Create a figure for the plot
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# Loop through the files in the data directory
for year in range(syear, eyear + 1):
    for month in range(1, 13):
        
        # Construct the file name based on the year and month
        file_name = f'{case_name}.clm2.h0.{year}-{month:02d}.nc'
        file_path = os.path.join(case_dir, file_name)
        
        # Print the file being processed
        print(f"Processing file: {file_path}")

        # Open the NetCDF file
        ds = xr.open_dataset(file_path)
        
        # Extract the variable
        if var in ds:
            print(f"{var} found in {file_name}")
            
            # Process flux variables (sum over time) or state variables (average over time)
            if is_flux:
                var_current = ds[var].mean(dim='time') * 60 * 60 * 24 * days_in_month(month)
            else:
                var_current = ds[var].mean(dim='time') * days_in_month(month) / 365  
            
            print(f"{var} data shape: {var_current.shape}")

            # Initialize or accumulate data
            if 'var_data' not in locals():
                var_data = var_current                
            else:
                var_data += var_current
    
        # Close the dataset after processing
        ds.close()

# Average over the years for state variables or keep the sum for flux variables
var_data = var_data / (eyear - syear + 1)

print(var_data)
print(f"{var} data shape: {var_data.shape}")

# Subset the data for Scandinavia
var_scandinavia = var_data.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

# Print the shape of the subsetted data
print(f"Subsetted {var} data shape: {var_scandinavia.shape}")

# Plot the data
im = var_scandinavia.plot(ax=ax, cmap=cmap, add_colorbar=True, transform=ccrs.PlateCarree())
ax.set_title(f"{var} for CLM-FATES Model", fontsize=20)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
im.set_clim(clim_min, clim_max)
ax.coastlines(resolution='50m', color='black', linewidth=1)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)

# Save the figure
output_file = os.path.join("figs/Trendy25/", f"{var.lower()}_scandinavia.png")
plt.tight_layout()
plt.savefig(output_file, dpi=300)
print(f"Figure saved as {output_file}")

