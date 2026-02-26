import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import cftime
import sys
import os

# Input arguments:
#   case_name_in (mandatory)
#   out_dir (optional, default='figs')
#   case_path (optional, overrides auto case directory lookup)
out_dir = 'figs'
case_path = None

if len(sys.argv) < 2:
    print("Usage: python Plot_NBP.py <case_name_in> [out_dir] [case_path]")
    sys.exit(1)

case_name_in = sys.argv[1]
if len(sys.argv) > 2:
    out_dir = sys.argv[2]
if len(sys.argv) > 3:
    case_path = sys.argv[3]

print('Starting Plot NBP')

# if outpath does not exist, create it
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

calc_annual = False
read_last_n = False
read_first_n = False
nfiles = 120

#Additional plots
plot_total = True

# Conversion factors
km2_to_m2 = 1e6
g_to_Gt = 1e-15
kg_to_Gt = 1e-12
CO2_to_C = 12.011 / 44.01  # Conversion factor from CO2 to C
seconds_per_year = 365 * 24 * 60 * 60
days_per_year = 365
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
sec_per_month = [days * 24 * 60 * 60 for days in days_per_month]

# Obs values (from https://www.ilamb.org/CMIP6/historical/):
obs_data = {
    'biomass': {
        'ESACCI': 626.0,     # PgC
        'GLOBCARBON': 364.0, # PgC
    },
    'cSoil': {
        'HWSD': 1330.0,      # PgC
    }
}

# Mapping from simulation variable names to ILAMB observation variable names
sim_to_ilamb = {
    "FATES_VEGC": "biomass",
    "TOTSOMC": "cSoil"
}



variables = ["TWS", "TLAI", "TOTSOMC", "TOT_WOODPRODC","TOT_WOODPRODC_LOSS", "FATES_SEEDS_IN_EXTERN_EL", "TOTLITC", "SOM_C_LEACHED",
            "FATES_VEGC", "FATES_LITTER_AG_CWD_EL","FATES_LITTER_AG_FINE_EL","FATES_LITTER_BG_CWD_EL","FATES_LITTER_BG_FINE_EL",  
            "FCO2", "FATES_NEP", "FATES_GRAZING", "FATES_FIRE_CLOSS", "FATES_SEED_BANK"]
stabilitycheck_vars = ["TLAI", "TWS", "TOTSOMC", "FATES_VEGC", "TOT_WOODPRODC", "TOTLITC",
                    "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", 
                    "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"
                    ]
stock_vars = ["TOTSOMC", "FATES_VEGC", "TOT_WOODPRODC", "TOTLITC", "FATES_SEED_BANK",
            "FATES_LITTER_AG_CWD_EL","FATES_LITTER_BG_CWD_EL"]

total_flux_vars = ["FATES_NEP", "FATES_GRAZING", "FATES_FIRE_CLOSS", "TOT_WOODPRODC_LOSS", "SOM_C_LEACHED"]

case_name = case_name_in

if case_path:
    #case_dir =  f"{case_path}/{case_name}/lnd/hist"
    case_dir =  f"{case_path}/{case_name}/run/"
    print(f'Using user-provided case directory: {case_dir}')
else:
    # starts with n, use NorESM folder, if it starts with i, use CTSM folder
    if case_name.__contains__('TRENDY2025'):
        case_dir = f"/nird/datalake/NS9560K/kjetisaa/TRENDY25/{case_name}/lnd/hist"
        print(f'Using TRENDY2025 case directory: {case_dir}')
    elif case_name.startswith('n'):
        case_dir = f"/nird/datalake/NS9560K/noresm3/cases/{case_name}/lnd/hist"
        print(f'Using NorESM case directory: {case_dir}')
    elif case_name.startswith('i'):
        case_dir = f"/nird/datalake/NS9560K/kjetisaa/{case_name}/lnd/hist"
        print(f'Using CTSM case directory: {case_dir}')
    else:
        print("Could not infer case directory from case_name. Provide case_path as third argument.")
        print("Usage: python Plot_NBP.py <case_name_in> [out_dir] [case_path]")
        sys.exit(1)

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/*.clm2.h0a*-*.nc'))
if not timeseries_files:
    print(f'No timeseries files found in directory: {case_dir}. Testing h0 files instead.')
    timeseries_files = sorted(glob.glob(f'{case_dir}/*.clm2.h0.*.nc'))
    if not timeseries_files:
        print(f'No timeseries files found in directory: {case_dir}. Exiting.')    
        sys.exit(1)

# Initialize a dictionary to store the accumulated results
results = {var: [] for var in variables}
time = []
co2_conc = []
fates_LU_area = []
if calc_annual:
    timeseries_files = timeseries_files[:-1]
if read_last_n:
    timeseries_files = timeseries_files[-nfiles:]
elif read_first_n:
    timeseries_files = timeseries_files[:nfiles]

# Loop through all timeseries files
first_file=True
for filename in timeseries_files[:]:
    print(f"Processing file: {filename}")
    mon=filename.split('.')[-2][-2:]
    mon_index=int(mon)-1
    print(f"Month index: {mon_index}")
    try:
        with xr.open_dataset(filename, engine='netcdf4') as data:            
            if first_file:
                area = data["area"].fillna(0)
                landfrac = data["landfrac"].fillna(0)               
                fates_fraction = data["FATES_FRACTION"].fillna(0) 
                # squeeze to remove any singleton dimensions
                fates_fraction = fates_fraction.squeeze()
            for var in variables:
                if var in data:
                    var_data_in = data[var]                                        
                    var_data = var_data_in
                    if var.startswith("FATES_") :
                        var_data = var_data * fates_fraction
                    # Determine spatial dimensions
                    if "lndgrid" in var_data.dims:
                        spatial_dims = ("lndgrid",)
                    elif "lat" in var_data.dims and "lon" in var_data.dims:
                        spatial_dims = ("lat", "lon")
                    else:
                        spatial_dims = None
                    if first_file:
                        print(f'Reading {var}, with units: {data[var].units}')
                    # Carbon pools: multiply by area and landfrac (convert m^2 to km^2)
                    area_land_m2 = area * km2_to_m2 * landfrac
                    if data[var].units == "gC/m^2":
                        if spatial_dims:
                            total = (var_data * area_land_m2 * g_to_Gt).sum(dim=spatial_dims)
                        if 'fates_levelem' in var_data.dims:
                            #print(f"Summing over fates_levelem for {var}")
                            total = total.sum(dim='fates_levelem')
                        results[var].append(total.values)
                    elif data[var].units == "kg m-2":
                        #print(f"Processing {var} with units {data[var].units} and dimensions {var_data.dims}")                        
                        if spatial_dims:
                            total = (var_data * area_land_m2 * kg_to_Gt).sum(dim=spatial_dims)
                        if 'fates_levelem' in var_data.dims:
                            #print(f"Summing over fates_levelem for {var}")
                            total = total.sum(dim='fates_levelem')
                        results[var].append(total.values)
                    elif data[var].units == "gC/m^2/s":
                        #print(f"Processing {var} with units {data[var].units} and dimensions {var_data.dims}")
                        if spatial_dims:
                            total = (var_data * area_land_m2 * g_to_Gt * sec_per_month[mon_index]).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg m-2 yr-1":
                        if spatial_dims:
                            total = (var_data * area_land_m2 * kg_to_Gt).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg m-2 s-1":
                        if spatial_dims:
                            total = (var_data * area_land_m2 * kg_to_Gt * sec_per_month[mon_index]).sum(dim=spatial_dims)
                        if 'fates_levelem' in var_data.dims:
                            total = total.sum(dim='fates_levelem')
                        results[var].append(total.values)                    
                    elif data[var].units == "kgCO2/m2/s":
                        if spatial_dims:
                            total = (var_data * area_land_m2 * kg_to_Gt * CO2_to_C * sec_per_month[mon_index]).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg s-1" :
                        if spatial_dims:
                            total = (var_data * kg_to_Gt * sec_per_month[mon_index] ).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    else:
                        if spatial_dims:
                            total = (var_data * area * landfrac).sum(dim=spatial_dims) / (area * landfrac).sum(dim=spatial_dims)
                        results[var].append(total.values)
            first_file = False
            time.append(data['time'].values)
    except FileNotFoundError:
        print(f"File not found: {filename}")
    except ValueError as e:
        print(f"Error reading the file: {e}")

# Convert lists to numpy arrays
for var in variables:
    print(f"Processing {var}")
    if len(results[var]) > 0:
        results[var] = np.concatenate(results[var])
    else:
        results[var] = np.array([])

if len(time) > 0:
    time = np.concatenate(time)
else:
    time = np.array([])

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

# Calculate deltas
deltas = {var: np.diff(results[var], prepend=np.nan) for var in variables}

# Plotting
fig, axes = plt.subplots(len(stabilitycheck_vars), 2, figsize=(20, 15), sharex='col')
fig.suptitle('Timeseries Data and Changes')

for i, var in enumerate(stabilitycheck_vars):
    # Determine correct units for title based on calculation
    if var in data and getattr(data[var], "units", None) == "gC/m^2":
        title_unit = "GtC"
    elif var in data and getattr(data[var], "units", None) == "gC/m^2/s":
        title_unit = "GtC/month"
    elif var in data and getattr(data[var], "units", None) == "kg m-2":
        title_unit = "GtC"
    elif var in data and getattr(data[var], "units", None) == "kg m-2 s-1":
        title_unit = "GtC/month"
    elif var in data and getattr(data[var], "units", None) == "kg m-2 yr-1":
        title_unit = "GtC/yr"
    elif var in data and getattr(data[var], "units", None) == "kgCO2/m2/s":
        title_unit = "GtC/month"        
    elif var in data and getattr(data[var], "units", None) == "kg s-1":
        title_unit = "GtC/month"
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

#axes[-1, 0].set_xlabel('Time')
#axes[-1, 1].set_xlabel('Time')

plt.tight_layout(rect=[0, 0, 1, 0.96])
if calc_annual:
    plt.savefig(f"{out_dir}/Stability_{case_name}_Annual.png")
    print(f'Wrote file: {out_dir}/Stability_{case_name}_Annual.png')
else:
    plt.savefig(f"{out_dir}/Stability_{case_name}.png")
    print(f'Wrote file: {out_dir}/Stability_{case_name}.png')


if plot_total:
    #calculate NBP as sum of nep
    results["NBP"] = results["FATES_NEP"] + results["FATES_SEEDS_IN_EXTERN_EL"] - results["FATES_GRAZING"] - results["FATES_FIRE_CLOSS"] - results["TOT_WOODPRODC_LOSS"]  - results["SOM_C_LEACHED"]  #TODO: shouldn't add this to NBP, but for balance check we do it here
    results["NBP"][0] = 0.0 #To be compareable with delta C stocks 
    results["FCO2"][0] = 0.0 #To be compareable with delta C stocks 
    results["TotalCarbon"] = sum(results[var] for var in stock_vars)
    deltas["TotalCarbon"] = np.diff(results["TotalCarbon"], prepend=np.nan)
  
    nyears = results["FCO2"].shape[0] / 12
    print(f'Number of years in FCO2: {nyears}')

    # Plotting
    figTot, axes = plt.subplots(2, 1, figsize=(25, 15), sharex='col')
    figTot.suptitle('Total carbon and cumulative changes')
    #axes[0].plot(np.sum([results[var] for var in stock_vars], axis=0))
    for var in stock_vars:
        axes[0].plot(results[var], label=var, linewidth=2) 
    axes[0].set_title('Carbon stocks')
    axes[0].set_ylabel('Carbon (GtC)')
    axes[0].grid(True)
    axes[0].legend()
 
    axes[1].plot(-np.cumsum([results["FCO2"]]), label="FCO2", linestyle='-', linewidth=2)
    axes[1].plot(np.cumsum([results["NBP"]]), label="NBP", linestyle='-')

    for obs_var in obs_data:
        if obs_var in sim_to_ilamb.values():
            for dataset, obs_val in obs_data[obs_var].items():
                # Find the corresponding simulation variable to get the correct x-axis length
                sim_var = [k for k, v in sim_to_ilamb.items() if v == obs_var]
                if sim_var and sim_var[0] in results:
                    x = np.arange(len(results[sim_var[0]]))
                else:
                    x = np.arange(len(results[stock_vars[0]]))
                axes[0].plot(x, np.full_like(x, obs_val, dtype=float), linestyle='--', linewidth=1.5, label=f"{obs_var} ({dataset})")
    axes[0].legend()

    axes[1].set_title('Cumulative NBP and FCO2')
    axes[1].set_ylabel('Carbon (GtC)')
    axes[1].grid(True)
    axes[1].legend()
    plt.savefig(f"{out_dir}/Stability_{case_name}_Total.png")
    print(f'Wrote file: {out_dir}/Stability_{case_name}_Total.png')    

print('Finished Check_Land_Carbon.py')
