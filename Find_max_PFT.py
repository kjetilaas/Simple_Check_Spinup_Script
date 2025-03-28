import xarray as xr
import numpy as np
import glob

print('Starting Find_Max_PFT')

# File names. Modify manually!
case_name = 'i1850.FATES-NOCOMP-noLU_cal15_twosteam.ne30pg3_tn14.noresm3_0_alpha01.20250326'
case_dir = f'/cluster/work/users/kjetisaa/archive/{case_name}/lnd/hist/'


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

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/{case_name}.clm2.h0.*-*.nc'))

if not timeseries_files:
    print("No files found. Please check the directory and file naming.")
    exit()

# Read the first file
first_file = timeseries_files[0]
print(f"Reading file: {first_file}")

try:
    with xr.open_dataset(first_file, engine='netcdf4') as data:
        # Extract latitude and longitude
        lats = data['lat'].values
        lons = data['lon'].values

        # Check if the PFT fraction variable exists
        if 'FATES_NOCOMP_PATCHAREA_PF' not in data:
            print("Variable 'FATES_NOCOMP_PATCHAREA_PF' not found in the dataset.")
            exit()

        # Extract PFT fractions
        pft_fractions = data['FATES_NOCOMP_PATCHAREA_PF'].values  # Shape: (time, fates_levpft, lndgrid)

        # Use the last time step for analysis
        pft_fractions_last = pft_fractions[-1, :, :]  # Shape: (fates_levpft, lndgrid)

        # Find the grid cell with the highest fraction for each PFT
        for pft_idx in range(pft_fractions_last.shape[0]):  
            # First calculate mean fraction over all the grid cells
            mean_fraction = np.nanmean(pft_fractions_last[pft_idx, :])

            # Only seach Northern  hemisphere
            pft_fractions_last[pft_idx, lats < 0] = np.nan
            # For actic C3 grass, only search in the Arctic
            if pft_idx == 11:  # Arctic C3 grass
                pft_fractions_last[pft_idx, lats < 55] = np.nan
            max_fraction_idx = np.nanargmax(pft_fractions_last[pft_idx, :])
            max_fraction = pft_fractions_last[pft_idx, max_fraction_idx]
            max_lat = lats[max_fraction_idx]
            max_lon = lons[max_fraction_idx]

            # Convert longitude to the range [-180, 180] if necessary
            if max_lon > 180:
                max_lon -= 360

            print(f"PFT {pft_names[pft_idx]}: Max Fraction = {max_fraction:.2f}, Lat = {max_lat:.2f}, Lon = {max_lon:.2f}, Mean Fraction = {mean_fraction:.4f}")

except FileNotFoundError:
    print(f"File not found: {first_file}")
except ValueError as e:
    print(f"Error reading the file: {e}")

print('Finished Find_Max_PFT')