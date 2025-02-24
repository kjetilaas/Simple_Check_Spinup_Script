import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

print('Starting CheckStability_Fates_spinup')

# File names. Modify manually!
#case_name = 'n1850.FATES-NOCOMP-AD.ne30_tn14.alpha08d.20250127_fixFincl1'
#case_name = 'n1850.FATES-NOCOMP-postAD.ne30pg3_tn14.alpha08d.20250129_fixFincl1'
case_name = 'i1850.FATES-NOCOMP.ne30pg3_tn14.alpha08d.20250218_coldstart_adrianas_tuning'

#case_dir =f'/nird/datalake/NS9560K/kjetisaa/{case_name}/lnd/hist'
case_dir =f'/cluster/work/users/kjetisaa/archive/{case_name}/lnd/hist/'
#case_dir =f'/cluster/work/users/kjetisaa/noresm/{case_name}/run/'

# Conversion factors
km2_to_m2 = 1e6
g_to_Gt = 1e-15
seconds_per_year = 365 * 24 * 60 * 60

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/{case_name}.clm2.h0.*-*.nc'))
variables = ["TOTSOMC", "TOTECOSYSC","FATES_VEGC", "FATES_GPP", "TWS", "TLAI"]

# Initialize a dictionary to store the accumulated results
results = {var: [] for var in variables}
time = []

# Loop through all timeseries files except the last one
first_file=True
for filename in timeseries_files[:-1]: 
    print(filename)   
    try:        
        with xr.open_dataset(filename, engine='netcdf4') as data:
            if first_file:        
                print('Reading area and landfrac')
                area = data["area"].fillna(0)
                landfrac = data["landfrac"].fillna(0)
                first_file=False
                print(filename)
            for var in variables:
                if var in data:
                    var_data = data[var]
                    # Carbon pools: multiply by area and landfrac (convert m^2 to km^2)
                    if  data[var].units=="gC/m^2":
                        # Carbon pools: multiply by area and landfrac (convert m^2 to km^2)
                        total = (var_data * area * km2_to_m2 * landfrac * g_to_Gt).sum(dim=("lndgrid"))
                        results[var].append(total.values)
                    elif data[var].units=="kg m-2 s-1":
                        # Fluxes: multiply by area, landfrac, and seconds per year
                        total = (var_data * area * km2_to_m2 * landfrac * g_to_Gt * seconds_per_year).sum(dim=("lndgrid"))                         
                        results[var].append(total.values)
                    else:
                        total = (var_data * area * landfrac).sum(dim=("lndgrid")) /(area * landfrac).sum(dim=("lndgrid"))                        
                        results[var].append(total.values)
                elif var == "FATES_VEGC":
                        var_data = data["TOTECOSYSC"] - data["TOTSOMC"]
                        total = (var_data * area * km2_to_m2 * landfrac * g_to_Gt).sum(dim=("lndgrid"))
                        results["FATES_VEGC"].append(total.values)
                        print("WARNING: calculating FATES_VEGC from TOTECOSYSC - TOTSOMC")
            time.append(data['time'].values)
    except FileNotFoundError:
        print(f"File not found: {filename}")
    except ValueError as e:
        print(f"Error reading the file: {e}")

# Convert lists to numpy arrays
for var in variables:
    results[var] = np.concatenate(results[var])

time = np.concatenate(time)


# Calculate deltas
deltas = {var: np.diff(results[var], prepend=np.nan) for var in variables}

# Plotting
fig, axes = plt.subplots(len(variables), 2, figsize=(20, 15), sharex='col')
fig.suptitle('Timeseries Data and Changes')


for i, var in enumerate(variables):
    # Original timeseries
    axes[i, 0].plot(results[var])
    axes[i, 0].set_title(var)
    axes[i, 0].set_ylabel(var)
    axes[i, 0].grid(True)
    
    #Todo: use time variable in plots

    # Delta timeseries
    axes[i, 1].plot(deltas[var])
    axes[i, 1].set_title(f"Delta {var}")
    axes[i, 1].set_ylabel(f"Delta {var}")
    axes[i, 1].grid(True)

axes[-1, 0].set_xlabel('Time')
axes[-1, 1].set_xlabel('Time')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(f"figs/Stability_{case_name}.png")

print('Finished CheckStability_Fates_spinup.py')