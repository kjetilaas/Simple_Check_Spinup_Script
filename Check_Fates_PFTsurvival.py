import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

print('Starting CheckFates_PFT_Variables')

# File names. Modify manually!
case_name = 'i1850.FATES-NOCOMP-LUH2-evolve.ne30pg3_tn14.alpha09a.20250321'
#case_name = 'n1850.FATES-NOCOMP-postAD.ne30pg3_tn14.alpha08d.20250206_fixFincl1'
case_dir = f'/cluster/work/users/kjetisaa/archive/{case_name}/lnd/hist/'

# Option to process only the last 10 years, only works for monthly data
process_last_10_years = True

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/{case_name}.clm2.h0.*-*.nc'))

if process_last_10_years:
    # Filter files to include only the last 10 years
    timeseries_files = timeseries_files[-120:]  # Assuming monthly data, 10 years = 120 months

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
time = []
fates_leaf_slamax = np.array([0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.012, 
                              0.03, 0.03, 0.012, 0.032, 0.05, 0.05, 0.05])

fates_leaf_slatop = np.array([0.012, 0.005, 0.024, 0.009, 0.03, 0.03, 0.012, 0.03, 
                              0.03, 0.01, 0.032, 0.027, 0.05, 0.05 ])

# Loop through all timeseries files
for filename in timeseries_files:
    print(filename)
    try:
        with xr.open_dataset(filename, engine='netcdf4') as data:
            if 'FATES_NOCOMP_PATCHAREA_PF' in data and 'FATES_LEAFC_PF' in data:
                patch_area = data['FATES_NOCOMP_PATCHAREA_PF']                
                leafc = data['FATES_LEAFC_PF']

                # Apply the SLA and filter by patch area
                adjusted_leafc = np.zeros_like(leafc)
                for i in range(len(fates_leaf_slatop)):
                    mask = patch_area[:, i] > 0
                    adjusted_leafc[:, i] = np.where(mask, (leafc[:, i] * fates_leaf_slatop[i] * 1000), 0)

                weighted_mean_leafc = (adjusted_leafc * patch_area).sum(dim='lndgrid') / patch_area.sum(dim='lndgrid')
                for i, pft in enumerate(pft_names):
                    results[pft].append(weighted_mean_leafc[:, i].values)
                time.append(data['time'].values)
    except FileNotFoundError:
        print(f"File not found: {filename}")
    except ValueError as e:
        print(f"Error reading the file: {e}")

# Convert lists to numpy arrays
for pft in pft_names:
    if results[pft]:  # Check if the list is not empty
        results[pft] = np.concatenate(results[pft])
    else:
        results[pft] = np.array([])

if time:  # Check if the list is not empty
    time = np.concatenate(time)
else:
    time = np.array([])

# Plotting
fig, axes = plt.subplots(5, 3, figsize=(20, 30), sharex='col')
fig.suptitle(f'Mean PFT level C_leaf*SLATOP, weighted by FATES_NOCOMP_PATCHAREA_PF', fontsize=20)

# Determine the y-axis limits from the first subplot
#y_min, y_max = np.inf, -np.inf
#for i, pft in enumerate(pft_names):
#    if results[pft].size > 0:  # Check if the array is not empty
#        y_min = min(y_min, results[pft].min())
#        y_max = max(y_max, results[pft].max())
y_min, y_max = 0, 10

# Plot each subplot with the same y-axis limits
for i, pft in enumerate(pft_names):
    ax = axes[i // 3, i % 3]
    ax.plot(results[pft])
    ax.set_title(pft.strip(), fontsize=14)
    ax.grid(True)
    ax.set_ylim(y_min, y_max)  # Set the same y-axis limits for all subplots
    if i // 3 == 4:
        ax.set_xlabel('Time', fontsize=12)
    if i % 3 == 0:
        ax.set_ylabel('Mean Product', fontsize=12)

plt.tight_layout(rect=[0, 0, 1, 0.96])
output_filename = f"figs/PFT_Timeseries_{case_name}"
if process_last_10_years:
    output_filename += "_last10years"
output_filename += ".png"
plt.savefig(output_filename)

print('Finished CheckFates_PFT_Variables')
