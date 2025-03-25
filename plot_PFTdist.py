import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

print('Starting plot_PFTdist')

# File names. Modify manually!
case_name = 'i1850.FATES-NOCOMP-LUH2-constant.f45_f45_mg37.alpha09a.20250319'
#case_name = 'n1850.FATES-NOCOMP-postAD.ne30pg3_tn14.alpha08d.20250206_fixFincl1'
case_dir = f'/cluster/work/users/kjetisaa/archive/{case_name}/lnd/hist/'

# Find the timeseries file
filename = f'{case_dir}/{case_name}.clm2.h0.0001-01.nc'

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

# Initialize a dictionary to store the accumulated results for each PFT
results = {pft: [] for pft in pft_names}
latitudes = []

print(filename)
try:
    with xr.open_dataset(filename, engine='netcdf4') as data:
        if 'FATES_NOCOMP_PATCHAREA_PF' in data:
            patch_area = data['FATES_NOCOMP_PATCHAREA_PF']
            lat = data['lat'].values
            mean_patch_area = patch_area.mean(dim='lon')
            for i, pft in enumerate(pft_names):
                results[pft] = mean_patch_area[:, i].values.squeeze()
            latitudes = lat
except FileNotFoundError:
    print(f"File not found: {filename}")
except ValueError as e:
    print(f"Error reading the file: {e}")


# Calculate the dominant PFT for each latitude
dominant_pft = np.argmax([results[pft] for pft in pft_names], axis=0).astype(float)+1
# if the maximum value is less than 0.20, set the dominant PFT to np.nan
dominant_pft[np.max([results[pft] for pft in pft_names], axis=0) < 0.20] = np.nan

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={'width_ratios': [3, 1]})
for i, pft in enumerate(pft_names):
    if np.nanmax(results[pft]) >= 0.20:
        ax1.plot(results[pft], latitudes, label=pft.strip())
ax1.set_title("PFT Distribution vs Latitude", fontsize=16)
ax1.set_xlabel('Mean Patch Area', fontsize=14)
ax1.set_ylabel('Latitude', fontsize=12)
ax1.legend(loc='upper right')
ax1.grid(True)  
ax1.tick_params(axis='both', which='major', labelsize=12)

# Plot dominant PFT number vs latitude
ax2.plot(dominant_pft, latitudes, 'o', markersize=5)
ax2.set_title("Dominant PFT vs Latitude", fontsize=16)
ax2.set_xlabel('PFT Number', fontsize=14)
ax2.set_ylabel('Latitude', fontsize=12)
ax2.grid(True)
ax2.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()

output_filename = f"figs/PFT_Distribution_vs_Latitude_{case_name}.png"
plt.savefig(output_filename)

#Print counts of dominant PFTs
for i, pft in enumerate(pft_names):
    count = np.count_nonzero(dominant_pft == i+1)
    if count > 0:
        print(f"{pft}: {count}")

print('Finished plot_PFTdist')