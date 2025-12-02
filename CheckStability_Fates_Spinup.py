import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

print('Starting CheckStability_Fates_spinup')

#Modify case name
#case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.ADspinup_TRENDY2025.20250813'
#case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.PostADspinup_TRENDY2025.20250814'
#case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.PostADspinup_TRENDY2025.20250816'

#case_name = 'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S0_TRENDY2025_pt1.202508011'
#case_name = 'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt1.202508018'

#case_name = 'i1850.ne16pg3_tn14.ctsm5.3.045_noresm_v14.CPLHIST_postADspinup_fullOutput.2025-11-18'
#case_name = 'i1850.ne16pg3_tn14.ctsm5.3.045_noresm_v15.CPLHIST_nograze.2025-11-20'
#case_name = 'i1850.ne16pg3_tn14.ctsm5.3.045_noresm_v15.CPLHIST_GrassTest.2025-11-19'
case_name = 'i1850.f45_f45_mg37.fatesnocomp.noresm3_0_beta06.CRUJRA2024_grazingtest_mods.2025-12-01'
#Modify case directory

#case_dir =f'/nird/datalake/NS9560K/kjetisaa/{case_name}/lnd/hist/'
case_dir =f'/nird/datalake/NS9560K/rosief/{case_name}/lnd/hist/'


out_dir = 'figs/'
annual_output=False

# if outpath does not exist, create it
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

calc_annual = False
read_last_n = True
read_first_n = False
nfiles = 120
plot_runmean = True  # Plot running mean of 12 months

#Additional plots
plot_total = True
plot_trendy_checks = True


if annual_output:
    calc_annual = False
    read_last_n = False
    read_first_n = False
    nfiles = 0
    plot_runmean = False  # No running mean for annual data
    plot_total = False
    plot_trendy_checks = False


# Conversion factors
km2_to_m2 = 1e6
g_to_Gt = 1e-15
kg_to_Gt = 1e-12
CO2_to_C = 12.011 / 44.01  # Conversion factor from CO2 to C
seconds_per_year = 365 * 24 * 60 * 60
if annual_output:
    variables = ["TOTSOMC", "FATES_NEP", "TWS", "H2OSNO", "TLAI"] #"TOTCOLC", "TOTSOMC_1m",
    stabilitycheck_vars = variables
else:
    variables = ["TOTSOMC","TWS", "TLAI", "TOTLITC","FATES_SEEDS_IN_EXTERN_EL", "FATES_LITTER_IN", "FATES_LITTER_OUT", 
                "TOT_WOODPRODC","FATES_VEGC", "FATES_LITTER_AG_CWD_EL","FATES_LITTER_AG_FINE_EL","FATES_LITTER_BG_CWD_EL","FATES_LITTER_BG_FINE_EL",  
                "FCO2", "FATES_NEP","CROPPROD1C_LOSS", "TOT_WOODPRODC_LOSS","FATES_GRAZING", "FATES_FIRE_CLOSS", 
                "FATES_HARVEST_WOODPROD_C_FLUX", "DWT_WOODPRODC_GAIN", "FATES_CBALANCE_ERROR"]
    stabilitycheck_vars = ["TLAI", "TWS", "TOTSOMC", "FATES_VEGC", "TOT_WOODPRODC",
                        "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", 
                        "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"
                        ]
    stock_vars = ["TOTSOMC", "FATES_VEGC", "TOT_WOODPRODC",
                "FATES_LITTER_AG_CWD_EL","FATES_LITTER_AG_FINE_EL","FATES_LITTER_BG_CWD_EL","FATES_LITTER_BG_FINE_EL"]

    trendy_vars = ["PBOT", "PCO2", "FATES_PATCHAREA_LU"]

#variables = ["TOTSOMC", "TOTECOSYSC","FATES_VEGC", "FATES_GPP", "TWS", "TLAI",'FCO2'] #"TOTCOLC", "TOTSOMC_1m", 
#variables = ["TOTSOMC", "TOTECOSYSC", "FATES_NONSTRUCTC", "FATES_STRUCTC", "TWS", "TLAI", "FATES_LITTER_AG_CWD_EL","FATES_LITTER_AG_FINE_EL","FATES_LITTER_BG_CWD_EL","FATES_LITTER_BG_FINE_EL", "FCO2", "FATES_NEP", "FATES_NPP"]

# Find all timeseries files
timeseries_files = sorted(glob.glob(f'{case_dir}/*.clm2.h0.*-*.nc'))

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



fates_landuseclass_name =["primaryland", "secondaryland", "rangeland", "pastureland", "cropland"]

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

            if plot_trendy_checks:
                if "PBOT" in data and "PCO2" in data:
                    pbot_data = data["PBOT"]
                    pco2_data = data["PCO2"]
                    co2_conc_fld = pco2_data / pbot_data
                    co2_conc.append(np.nanmean(co2_conc_fld.values)*1e6)
                    #print(f"CO2 concentration (ppm): {co2_conc[-1]}")
                if "FATES_PATCHAREA_LU" in data:
                    # Determine spatial dimensions
                    if "lndgrid" in data["FATES_PATCHAREA_LU"].dims:
                        spatial_dims = ("lndgrid",)
                    elif "lat" in data["FATES_PATCHAREA_LU"].dims and "lon" in data["FATES_PATCHAREA_LU"].dims:
                        spatial_dims = ("lat", "lon")
                    else:
                        spatial_dims = None
                    fates_patcharea = data["FATES_PATCHAREA_LU"]
                    fates_fraction = data["FATES_FRACTION"].fillna(0)
                    total = (fates_patcharea * area * landfrac * fates_fraction).sum(dim=spatial_dims)
                    fates_LU_area.append(total.values)
                    #print(f"FATES_PATCHAREA_LU total area: {total.values}")

            for var in variables:
                if var in data:
                    var_data_in = data[var]                                        
                    var_data = var_data_in
                    if not annual_output: #TODO: this is a hack because FATES_FRACTION is currently missing in annual files
                        fates_fraction = data["FATES_FRACTION"].fillna(0)
                        if var.startswith("FATES_") and var_data_in.dims == fates_fraction.dims:
                            if first_file:
                                print(f'Multiplying {var} by FATES_FRACTION')
                            #var_data = var_data_in * fates_fraction
                            var_data.attrs = var_data_in.attrs  # Preserve metadata/attributes
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
                            total = (var_data * area_land_m2 * g_to_Gt * seconds_per_year).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg m-2 yr-1":
                        if spatial_dims:
                            total = (var_data * area_land_m2 * kg_to_Gt).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg m-2 s-1":
                        if spatial_dims:
                            total = (var_data * area_land_m2 * kg_to_Gt * seconds_per_year).sum(dim=spatial_dims)
                        if 'fates_levelem' in var_data.dims:
                            #print(f"Summing over fates_levelem for {var}")
                            total = total.sum(dim='fates_levelem')
                        results[var].append(total.values)                    
                    elif data[var].units == "kgCO2/m2/s":
                        if spatial_dims:
                            total = (var_data * area_land_m2 * kg_to_Gt * CO2_to_C * seconds_per_year).sum(dim=spatial_dims)
                        results[var].append(total.values)
                    elif data[var].units == "kg s-1" :
                        if spatial_dims:
                            total = (var_data * kg_to_Gt * seconds_per_year).sum(dim=spatial_dims)
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
    #print(f"Processing {var}")
    results[var] = np.concatenate(results[var])
if plot_trendy_checks:
    fates_LU_area  = np.concatenate(fates_LU_area)
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
fig, axes = plt.subplots(len(stabilitycheck_vars), 2, figsize=(20, 15), sharex='col')
fig.suptitle('Timeseries Data and Changes')

for i, var in enumerate(stabilitycheck_vars):
    # Determine correct units for title based on calculation
    if var in data and getattr(data[var], "units", None) == "gC/m^2":
        title_unit = "GtC"
    elif var in data and getattr(data[var], "units", None) == "gC/m^2/s":
        title_unit = "GtC/yr"
    elif var in data and getattr(data[var], "units", None) == "kg m-2":
        title_unit = "GtC"
    elif var in data and getattr(data[var], "units", None) == "kg m-2 s-1":
        title_unit = "GtC/yr"
    elif var in data and getattr(data[var], "units", None) == "kg m-2 yr-1":
        title_unit = "GtC/yr"
    elif var in data and getattr(data[var], "units", None) == "kgCO2/m2/s":
        title_unit = "GtC/yr"        
    elif var in data and getattr(data[var], "units", None) == "kg s-1":
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


if plot_total:
    #calculate NBP as sum of nep
    results["NBP"] = results["FATES_NEP"] + results["FATES_SEEDS_IN_EXTERN_EL"] - results["FATES_GRAZING"] - results["FATES_FIRE_CLOSS"] - results["TOT_WOODPRODC_LOSS"]  
    results["NBP"][0] = 0.0 #To be compareable with delta C stocks 
    results["FCO2"][0] = 0.0 #To be compareable with delta C stocks 
    results["TotalCarbon"] = sum(results[var] for var in stock_vars)
    deltas["TotalCarbon"] = np.diff(results["TotalCarbon"], prepend=np.nan)
    if annual_output:
        #Calculate number of years from length of NBP
        nyears=results["NBP"].shape[0]  
    else:   
        nyears = results["NBP"].shape[0] / 12
    print(f'Number of years in NBP: {nyears}')

    # Plotting
    figTot, axes = plt.subplots(3, 1, figsize=(20, 20), sharex='col')
    figTot.suptitle('Total carbon stocks and changes')
    axes[0].plot(np.sum([results[var] for var in stock_vars], axis=0))
    for var in stock_vars:
        axes[0].plot(results[var], label=var)
    axes[0].set_title('Total Carbon Stock (GtC)')
    axes[0].set_ylabel('Carbon Stock (GtC)')
    axes[0].grid(True)
    axes[0].legend()

    axes[1].plot(np.sum([results[var] for var in stock_vars], axis=0)-np.sum([results[var][0] for var in stock_vars]), label="Total C Change")
    axes[1].plot(np.cumsum([results["NBP"]])/12, label="NBP", linestyle='-')
    axes[1].plot(-np.cumsum([results["FCO2"]])/12, label="FCO2", linestyle='-')
    for var in stock_vars:
        axes[1].plot(results[var]-results[var][0], label=var, linestyle='--')
    axes[1].set_title('Total Carbon Stock Change (GtC)')
    axes[1].set_ylabel('Carbon Stock (GtC)')
    axes[1].grid(True)
    axes[1].legend()



    #plot FCO2 in same plot:
    if plot_runmean:
        import pandas as pd        
        axes[2].plot(-pd.Series(results["FCO2"]).rolling(window=12).mean(), label="FCO2", linestyle='--')
        axes[2].plot(pd.Series(results["FATES_NEP"]).rolling(window=12).mean(), label="FATES_NEP", linestyle='--')
        axes[2].plot(pd.Series(results["FATES_GRAZING"]).rolling(window=12).mean(), label="FATES_GRAZING", linestyle='--')
        axes[2].plot(pd.Series(results["TOT_WOODPRODC_LOSS"]).rolling(window=12).mean(), label="TOT_WOODPRODC_LOSS", linestyle='--')
        axes[2].plot(pd.Series(results["CROPPROD1C_LOSS"]).rolling(window=12).mean(), label="CROPPROD1C_LOSS", linestyle='--')
        axes[2].plot(pd.Series(results["FATES_FIRE_CLOSS"]).rolling(window=12).mean(), label="FATES_FIRE_CLOSS", linestyle='--')
        axes[2].plot(pd.Series(results["FATES_NEP"]-results["FATES_GRAZING"]-results["FATES_FIRE_CLOSS"]-results["TOT_WOODPRODC_LOSS"]).rolling(window=12).mean(), label="NEP-GRAZ-FIRE-PRODLOSS", linestyle='--')
        axes[2].plot(pd.Series(results["NBP"]).rolling(window=12).mean(), label="NBP", linestyle='-', color='black', linewidth=2)
    else:
        axes[2].plot(np.sum([deltas[var] for var in stock_vars], axis=0)*12, label="Total Change")
        axes[2].plot(-results["FCO2"], label="FCO2", linestyle='--')
        axes[2].plot(results["FATES_NEP"], label="FATES_NEP", linestyle='--')
        axes[2].plot(results["FATES_GRAZING"], label="FATES_GRAZING", linestyle='--')
        axes[2].plot(results["TOT_WOODPRODC_LOSS"], label="TOT_WOODPRODC_LOSS", linestyle='--')
        axes[2].plot(results["CROPPROD1C_LOSS"], label="CROPPROD1C_LOSS", linestyle='--')
        axes[2].plot(results["FATES_FIRE_CLOSS"], label="FATES_FIRE_CLOSS", linestyle='--')
        axes[2].plot(results["FATES_NEP"]-results["FATES_GRAZING"]-results["FATES_FIRE_CLOSS"]-results["TOT_WOODPRODC_LOSS"], label="NEP-GRAZ-FIRE-PRODLOSS", linestyle='--')    
        axes[2].plot(results["NBP"], label="NBP", linestyle='-', color='black', linewidth=2)
    axes[2].legend()
    axes[2].set_title('Total Change in Carbon Stock (GtC/yr)')
    axes[2].set_ylabel('Change in Carbon Stock (GtC/yr)')
    axes[2].grid(True)
    plt.savefig(f"{out_dir}/Stability_{case_name}_Total.png")
    print(f'Wrote file: {out_dir}/Stability_{case_name}_Total.png')    
    print(f'Total Carbon Stock at the start: {results["TotalCarbon"][0]} GtC')
    print(f'Total Carbon Stock at the end: {results["TotalCarbon"][-1]} GtC')
    print(f'Total Change in Carbon Stock: {results["TotalCarbon"][-1] - results["TotalCarbon"][0]} GtC')
    for var in stock_vars:
        print(f'Carbon stock change in {var}: {(results[var][-1]) - (results[var][0])} GtC')
    #print(f'Accumulated NBP (NEP-fire-grazing-prodloss): {np.sum(results["FATES_NEP"]-results["FATES_GRAZING"]-results["FATES_FIRE_CLOSS"]-results["TOT_WOODPRODC_LOSS"])/12} GtC')
    print(f'Accumulated NBP (NEP-fire-grazing-prodloss): {np.mean(results["NBP"])*nyears} GtC')
    print(f'Accumulated NBP from FCO2: {np.mean(-results["FCO2"])*nyears} GtC')
    print(f'Accumulated FATES_NEP: {np.mean(results["FATES_NEP"])*nyears} GtC')
    print(f'Accumulated FATES_GRAZING: {np.mean(results["FATES_GRAZING"])*nyears} GtC')
    print(f'Accumulated FATES_FIRE_CLOSS: {np.mean(results["FATES_FIRE_CLOSS"])*nyears} GtC')
    print(f'Accumulated TOT_WOODPRODC_LOSS: {np.mean(results["TOT_WOODPRODC_LOSS"])*nyears} GtC')
    print(f'Accumulated FATES_LITTER_IN: {np.mean(results["FATES_LITTER_IN"])*nyears} GtC')
    print(f'Accumulated FATES_LITTER_OUT: {np.mean(results["FATES_LITTER_OUT"])*nyears} GtC')
    print(f'Change in litter from in-out: {np.mean(results["FATES_LITTER_IN"])*nyears - np.mean(results["FATES_LITTER_OUT"])*nyears} GtC')

if plot_trendy_checks:
    import matplotlib.ticker as mticker

    #Read CO2 from forcing file
#    co2_forcing_file = "/cluster/work/users/kjetisaa/Trendy_2025_forcing/CO2field/global_co2_ann_1700_2024.txt"
#    co2_forcing_data = np.loadtxt(co2_forcing_file, skiprows=1)
    # Convert co2_conc to numpy array if not already
    co2_conc = np.array(co2_conc)
    co2_conc_annual = np.array([co2_conc[i*12:(i+1)*12].mean() for i in range(int(nyears))])

    figTrendy, axes = plt.subplots(3, 1, figsize=(20, 20), sharex='col')
    
    #Plot CO2
    axes[0].plot(co2_conc, label="CO2 Concentration FATES", linestyle='-')
 #   axes[0].plot(np.arange(6, co2_conc.shape[0], 12), co2_forcing_data[0:co2_conc.shape[0]//12, 1], label="CO2 Concentration Forcing", marker='*', linestyle='None')
    axes[0].set_title('CO2 Concentration Over Time')
    axes[0].set_ylabel('CO2 Concentration (ppm)')
    axes[0].grid(True)
    axes[0].legend()
    axes[0].yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.2f}'))

    #Read LUH land use data
    print('Reading LUH land use data')
    #luh_landuse_states = "/cluster/work/users/kjetisaa/Trendy_2025_forcing/luh/states4.nc"
    #luh_landuse_static = "/cluster/work/users/kjetisaa/Trendy_2025_forcing/luh/staticData_quarterdeg.nc"

    #luh_landuse_static_data = xr.open_dataset(luh_landuse_static, engine='netcdf4')
    #luh_carea = luh_landuse_static_data["carea"].values
    #luh_icwtr = luh_landuse_static_data["icwtr"].values

    #luh_landuse_data = xr.open_dataset(luh_landuse_states, engine='netcdf4', decode_times=False)
    #luh_primforest = np.nansum(luh_landuse_data["primf"].values * luh_carea * luh_icwtr, axis=(1, 2)) 
    #luh_secforest = np.nansum(luh_landuse_data["secdf"].values * luh_carea * luh_icwtr, axis=(1, 2))
    #luh_rangeland = np.nansum(luh_landuse_data["range"].values * luh_carea * luh_icwtr, axis=(1, 2))
    #luh_pastr = np.nansum(luh_landuse_data["pastr"].values * luh_carea * luh_icwtr, axis=(1, 2))
    #luh_cropland = np.nansum((luh_landuse_data["c3ann"].values + luh_landuse_data["c4ann"].values + \
    #    luh_landuse_data["c3per"].values + luh_landuse_data["c4per"].values + luh_landuse_data["c3nfx"].values) * luh_carea * luh_icwtr, axis=(1, 2))
    #luh_urban = np.nansum(luh_landuse_data["urban"].values * luh_carea * luh_icwtr, axis=(1, 2))
    #luh_startyrnr = 1701-850

    #plot land use
    for lunr in range(len(fates_landuseclass_name)):
        axes[1].plot(fates_LU_area[:, lunr], label=fates_landuseclass_name[lunr], linestyle='-')
    #axes[1].plot(np.arange(6, co2_conc.shape[0], 12), luh_primforest[luh_startyrnr:luh_startyrnr+nyears], label='LUH PrimForest', linestyle='--')
    #axes[1].plot(np.arange(6, co2_conc.shape[0], 12), luh_secforest[luh_startyrnr:luh_startyrnr+nyears], label='LUH SecForest', linestyle='--')
    #axes[1].plot(np.arange(6, co2_conc.shape[0], 12), luh_rangeland[luh_startyrnr:luh_startyrnr+nyears], label='LUH Range', linestyle='--')
    #axes[1].plot(np.arange(6, co2_conc.shape[0], 12), luh_pastr[luh_startyrnr:luh_startyrnr+nyears], label='LUH Pasture', linestyle='--')
    #axes[1].plot(np.arange(6, co2_conc.shape[0], 12), luh_cropland[luh_startyrnr:luh_startyrnr+nyears], label='LUH Cropland', linestyle='--')
    #axes[1].plot(np.arange(6, co2_conc.shape[0], 12), luh_urban[luh_startyrnr:luh_startyrnr+nyears], label='LUH Urban', linestyle='--')

    axes[1].set_title('FATES_PATCHAREA_LU')
    axes[1].set_ylabel('FATES_PATCHAREA_LU (km2)')
    axes[1].grid(True)
    axes[1].legend()

    for lunr in range(len(fates_landuseclass_name)):
        axes[2].plot(fates_LU_area[:, lunr]-fates_LU_area[0, lunr], label=fates_landuseclass_name[lunr], linestyle='-')
    axes[2].set_title('Change in FATES_PATCHAREA_LU')
    axes[2].set_ylabel('FATES_PATCHAREA_LU (km2)')
    axes[2].grid(True)
    axes[2].legend()

    plt.savefig(f"{out_dir}/Stability_{case_name}_TrendyChecks.png")
    print(f'Wrote file: {out_dir}/Stability_{case_name}_TrendyChecks.png')
print('Finished CheckStability_Fates_spinup.py')
