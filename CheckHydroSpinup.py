import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

print('Starting CheckHydroSpinup.py')

#Modify case name
#case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.ADspinup_TRENDY2025.20250813'
#case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.PostADspinup_TRENDY2025.20250814'
#case_name = 'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.PostADspinup_TRENDY2025.20250816'

#case_name = 'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S0_TRENDY2025_pt1.202508011'
#case_name = 'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt1.202508018'

case_name = 'i1850.ne30pg3_tn14.fatesnocomp.ctsm5.3.045_noresm_v10.CPLHIST_noLU_coldstart_v25u.20250905'
#Modify case directory

case_dir =f'/nird/datalake/NS9560K/noresm3/cases/{case_name}/lnd/hist/'

out_dir = 'figs/NorESM_key_experiments/'

nc_outdir=f'/nird/datalake/NS9560K/kjetisaa/PostProcessed/{case_name}/'

# if outpath does not exist, create it
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
if not os.path.exists(nc_outdir):
    os.makedirs(nc_outdir)

plot_last_n_years = 50
convert_to_Pg = True

# Conversion factors
kg_to_Pg = 1e-12  # Conversion factor from kg to Pg
km2_to_m2 = 1e6
g_to_Gt = 1e-15
kg_to_Gt = 1e-12
CO2_to_C = 12.011 / 44.01  # Conversion factor from CO2 to C
seconds_per_year = 365 * 24 * 60 * 60

variables = ["TWS", "H2OSNO", "H2OSFC", "H2OCAN", "TOTEXICE_VOL", "TOTSOILICE", "TOTSOILLIQ"] #"TOTCOLC", "TOTSOMC_1m",

#Read landarea_m2, landareafates_m2 and nbedrock
landarea_file = f'{nc_outdir}/{case_name}_landarea_m2.nc'
landareafates_file = f'{nc_outdir}/{case_name}_landareafates_m2.nc'
nbedrock_file = f'{nc_outdir}/{case_name}_nbedrock.nc'

nbedrock=xr.open_dataarray(nbedrock_file)
landarea=xr.open_dataarray(landarea_file)
landareafates_data=xr.open_dataset(landareafates_file)
landarea=landareafates_data['landareafates_m2']

results = {var: [] for var in variables}

for var in variables:
    print(f'Reading variable: {var}')
    filename=f'{nc_outdir}/{case_name}_{var}.nc'
    ds=xr.open_dataset(filename)
    data=ds[var]
    data=data[6:-1:12,:] #Select only January data
    results[var] = data    
    print(f'Read variable {var} with shape {data.shape}')

#Plot results
fig, axes = plt.subplots(len(variables), 1, figsize=(10, 5*len(variables)))
for i, var in enumerate(variables):
    print(f'Plotting variable: {var}')
    data = results[var]
    plotyears=np.arange(len(data))
    if plot_last_n_years > 0:
        data = data[-plot_last_n_years:, :]
        plotyears = plotyears[-plot_last_n_years:]
    #plot sum of data times landarea over all grid cells
    if convert_to_Pg:
        if var == 'TOTEXICE_VOL':
            axes[i].plot(plotyears,np.sum(data * landarea, axis=(1))*8600*kg_to_Pg, label=var)  # Convert to Pg
        else:
            axes[i].plot(plotyears,np.sum(data * landarea, axis=(1))*kg_to_Pg, label=var)  # Convert to Pg
    else:
        if var == 'TOTEXICE_VOL':
            axes[i].plot(plotyears,np.sum(data * landarea, axis=(1))*8600/np.sum(landarea), label=var)  # Convert to kg
        else:
            axes[i].plot(plotyears,np.sum(data * landarea, axis=(1))/np.sum(landarea), label=var)  # Area-weighted average
    axes[i].legend()        
    axes[i].set_xlabel('Time (years)')    
    if convert_to_Pg:
        unit='Pg'
    else:
        unit='mm'
    axes[i].set_title(f'{var} ({unit})')
    axes[i].set_ylabel(f'Total ({unit})')        
plt.tight_layout()

if plot_last_n_years > 0:
    plt.suptitle(f'{case_name} - Last {plot_last_n_years} years', y=1.02)
    plt.savefig(f'{out_dir}/{case_name}_hydro_spinup_last_{plot_last_n_years}_years_{unit}.png')
    print(f'Saved figure to {out_dir}/{case_name}_hydro_spinup_last_{plot_last_n_years}_years_{unit}.png')
else:
    plt.suptitle(f'{case_name}', y=1.02)
    plt.savefig(f'{out_dir}/{case_name}_hydro_spinup_{unit}.png')
    print(f'Saved figure to {out_dir}/{case_name}_hydro_spinup_{unit}.png')
plt.close()
