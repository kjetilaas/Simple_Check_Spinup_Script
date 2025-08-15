import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

print('Starting CheckStability_Fates_spinup')

#Modify case name
case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.ADspinup_TRENDY2025.20250813'
#case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.PostADspinup_TRENDY2025.20250814'

#Modify case directory
#case_dir =f'/nird/datalake/NS9560K/kjetisaa/{case_name}/lnd/hist'
#case_dir =f'/cluster/work/users/kjetisaa/archive/{case_name}/lnd/hist/'
case_dir =f'/cluster/work/users/kjetisaa/noresm/{case_name}/run/'
#case_dir =f'/cluster/work/users/kjetisaa/archive/{case_name}/'

out_dir = 'figs/Trendy25/'

# if outpath does not exist, create it
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

calc_annual = False
read_last_n = False
nfiles = 240

# Conversion factors
km2_to_m2 = 1e6
g_to_Gt = 1e-15
kg_to_Gt = 1e-12
CO2_to_C = 12 / 44  # Conversion factor from CO2 to C
seconds_per_year = 365 * 24 * 60 * 60

#variables = ["TOTSOMC", "TOTECOSYSC","FATES_VEGC", "FATES_GPP", "TWS", "TLAI",'FCO2'] #"TOTCOLC", "TOTSOMC_1m", 
variables = ["TOTSOMC", "TOTECOSYSC", "FATES_NONSTRUCTC", "FATES_STRUCTC", "TWS", "TLAI", "FATES_LITTER_AG_CWD_EL","FATES_LITTER_AG_FINE_EL","FATES_LITTER_BG_CWD_EL","FATES_LITTER_BG_FINE_EL", "FCO2", "FATES_NEP", "FATES_NPP"]

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/*.clm2.h0.*-*.nc'))

# Initialize a dictionary to store the accumulated results
results = {var: [] for var in variables}
time = []

if calc_annual:
    timeseries_files = timeseries_files[:-1]
if read_last_n:
    timeseries_files = timeseries_files[-nfiles:]


# Loop through all timeseries files except the last one
first_file=True
for filename in timeseries_files[:]:
    print(f"Processing file: {filename}")
    try:
        with xr.open_dataset(filename, engine='netcdf4') as data:            
            if first_file:
                area = data["area"].fillna(0)
                landfrac = data["landfrac"].fillna(0)
                first_file = False
            for var in variables:
                if var in data:
                    var_data = data[var]
                    # Determine spatial dimensions
                    if "lndgrid" in var_data.dims:
                        spatial_dims = ("lndgrid",)
                    elif "lat" in var_data.dims and "lon" in var_data.dims:
                        spatial_dims = ("lat", "lon")
                    else:
                        spatial_dims = None

                    # Carbon pools: multiply by area and landfrac (convert m^2 to km^2)
                    if data[var].units == "gC/m^2":
                        if spatial_dims:
                            total = (var_data * area * km2_to_m2 * landfrac * g_to_Gt).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg m-2":
                        if spatial_dims:
                            total = (var_data * area * km2_to_m2 * landfrac * kg_to_Gt).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg m-2 s-1" or data[var].units == "kgCO2/m2/s":
                        if spatial_dims:
                            total = (var_data * area * km2_to_m2 * landfrac * kg_to_Gt * CO2_to_C * seconds_per_year).sum(dim=spatial_dims)

                        results[var].append(total.values)
                    else:
                        if spatial_dims:
                            total = (var_data * area * landfrac).sum(dim=spatial_dims) / (area * landfrac).sum(dim=spatial_dims)

                        results[var].append(total.values)

            time.append(data['time'].values)
    except FileNotFoundError:
        print(f"File not found: {filename}")
    except ValueError as e:
        print(f"Error reading the file: {e}")

# Convert lists to numpy arrays
for var in variables:
    results[var] = np.concatenate(results[var])

time = np.concatenate(time)

if calc_annual:
    annual_means = {var: [] for var in variables}
    annual_time = []
    for year in np.unique(time.astype('datetime64[Y]').astype(int)):
        mask = (time.astype('datetime64[Y]').astype(int) == year)
        annual_time.append(time[mask][0])
        for var in variables:
            annual_means[var].append(results[var][mask].mean())
    time = np.array(annual_time)
    for var in variables:
        results[var] = np.array(annual_means[var])



deltas = {var: np.diff(results[var], prepend=np.nan) for var in variables}

# Plotting
fig, axes = plt.subplots(len(variables), 2, figsize=(20, 15), sharex='col')
fig.suptitle('Timeseries Data and Changes')



for i, var in enumerate(variables):
    # Determine correct units for title based on calculation
    if var == "FATES_VEGC" or (var in data and getattr(data[var], "units", None) == "gC/m^2"):
        title_unit = "GtC"
    elif var in data and getattr(data[var], "units", None) == "kg m-2":
        title_unit = "GtC"
    elif var in data and getattr(data[var], "units", None) == "kg m-2 s-1":
        title_unit = "GtC/yr"
    elif var in data and getattr(data[var], "units", None) == "kgCO2/m2/s":
        title_unit = "GtC/yr"        
    elif var in data and getattr(data[var], "units", None) is not None:
        title_unit = getattr(data[var], "units")
    else:
        title_unit = "unitless"
    # Original timeseries
    axes[i, 0].plot(results[var])
    axes[i, 0].set_title(f"{var} ({title_unit})")
    axes[i, 0].set_ylabel("")
    axes[i, 0].grid(True)
    axes[i, 1].plot(deltas[var])
    axes[i, 1].set_title(f"Delta {var} ({title_unit})")
    axes[i, 1].set_ylabel("")
    axes[i, 1].grid(True)

axes[-1, 0].set_xlabel('Time')
axes[-1, 1].set_xlabel('Time')

plt.tight_layout(rect=[0, 0, 1, 0.96])
if calc_annual:
    plt.savefig(f"{out_dir}/Stability_{case_name}_Annual.png")
    print(f'Wrote file: {out_dir}/Stability_{case_name}_Annual.png')
else:
    plt.savefig(f"{out_dir}/Stability_{case_name}.png")
    print(f'Wrote file: {out_dir}/Stability_{case_name}.png')

print('Finished CheckStability_Fates_spinup.py')
