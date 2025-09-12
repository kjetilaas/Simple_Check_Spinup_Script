import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

#case_name= 'i1850.f45_f45_mg37.fatesnocomp.noresm3_0_alpha03c.potveg_coldstart.20250602'
#case_name= 'i1850.f45_f45_mg37.fatesnocomp.noresm3_0_alpha03c.constLU_fromPotveg.20250603'
#case_name= 'ihist.f45_f45_mg37.fatesnocomp.noresm3_0_alpha03c.hist_fromconstLU.20250605'
case_name= 'ihist.f45_f45_mg37.fatesnocomp.noresm3_0_alpha03c.hist_fromconstLU_pt2.20250610'

case_dir = f'/nird/datalake/NS9560K/kjetisaa/FATES_LU_FullTest1_4x5/{case_name}/lnd/hist/'

syear = 1901

# List your NetCDF files (adjust the path as needed)
nc_files = sorted(glob.glob(f'{case_dir}*.nc'))

# Option to process only the last 10 years, only works for monthly data
process_last_2_years = False

if process_last_2_years:
    # Filter files to include only the last 2 years
    nc_files = nc_files[-24:]  # Assuming monthly data, 2 years = 24 months

# Initialize lists to store timeseries
totsomc_ts = []
fates_vegc_ts = []
fates_litter_cwd_eldc_ts = []
times_ts = []

for f in nc_files:
    print(f"Processing file: {f}")
    ds = xr.open_dataset(f)
    # Read area (convert from km^2 to m^2)
    area = ds['area'] * 1e6  # km^2 to m^2
    landfrac = ds['landfrac'] if 'landfrac' in ds else 1.0

    # TOTSOMC: gC/m^2, convert to kgC (optional)
    totsomc = ds['TOTSOMC'] * area * landfrac  # [time, lat, lon]
    totsomc_global = totsomc.sum(dim=['lat', 'lon']) / 1e15 # Convert from gC to PgC

    # FATES_VEGC: kgC/m^2
    fates_vegc = ds['FATES_VEGC'] * area * landfrac  # [time, lat, lon]
    fates_vegc_global = fates_vegc.sum(dim=['lat', 'lon']) / 1e12  # Convert from kgC to PgC

    # FATES_LITTER_CWD_ELDC: sum over fates_levelcwd, then area/landfrac
    fates_litter = ds['FATES_LITTER_CWD_ELDC'].sum(dim='fates_levelcwd') * area * landfrac  # [time, lat, lon]
    fates_litter_global = fates_litter.sum(dim=['lat', 'lon'])/ 1e12  # Convert from kgC to PgC

    # times
    times = ds['time'].values if 'time' in ds else None
    
    
    # Store timeseries
    totsomc_ts.append(totsomc_global.values)
    fates_vegc_ts.append(fates_vegc_global.values)
    fates_litter_cwd_eldc_ts.append(fates_litter_global.values)
    times_ts.append(times)

    ds.close()

# Concatenate timeseries if multiple files
totsomc_ts = np.concatenate(totsomc_ts)
fates_vegc_ts = np.concatenate(fates_vegc_ts)
fates_litter_cwd_eldc_ts = np.concatenate(fates_litter_cwd_eldc_ts)
times_ts = np.concatenate(times_ts)

#print(times_ts[0], times_ts[-1])
#print(times_ts.shape)
#print(totsomc_ts)

plot_times = np.arange(syear, syear + len(totsomc_ts)/12 , 1/12 )  # Assuming monthly data, this will be the time index
print(f"Plotting times from {plot_times[0]} to {plot_times[-1]}")

# Plotting
fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharex=True)
axs[0].plot(plot_times, totsomc_ts)
axs[0].set_title("Global TOTSOMC")
axs[0].set_ylabel("kg C")

axs[1].plot(plot_times, fates_vegc_ts)
axs[1].set_title("Global FATES_VEGC")
axs[1].set_ylabel("kg C")

axs[2].plot(plot_times, fates_litter_cwd_eldc_ts)
axs[2].set_title("Global FATES_LITTER_CWD_ELDC")
axs[2].set_ylabel("kg C")

for ax in axs:
    ax.set_xlabel("Time")

plt.tight_layout()

output_filename = f"figs/Global_Timeseries_{case_name}"
output_filename += ".png"
plt.savefig(output_filename)

print('Finished processing and plotting global timeseries.')