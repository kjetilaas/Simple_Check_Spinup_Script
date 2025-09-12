import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import os


# --- Main function and main guard ---
def main():
    print('Starting CheckFates_Points')
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
    for case_name in CASE_NAMES:
        print(f"Processing case: {case_name}")
        case_dir = get_case_dir(case_name)
        timeseries_files = get_timeseries_files(case_dir, case_name)
        timeseries_files = filter_timeseries_files(timeseries_files)
        variables = get_variables()
        obs_datasets, obs_results = get_obs_datasets_and_results()
        locations = get_locations()
        obs_results = process_obs_datasets(obs_datasets, locations, obs_results)
        results, units, time = process_model_files(timeseries_files, variables, locations)
        plot_variables = get_plot_variables(variables, results, locations)
        results = add_cue_if_possible(results, locations)
        if CALC_ANNUAL:
            results, time = calc_annual_means(results, time, plot_variables, locations)
        plot_results(results, units, time, plot_variables, locations, obs_results, case_name)
    print('Finished CheckFates_Points')

if __name__ == "__main__":
    main()

# --- Main function and main guard moved to end ---

# --- All helper functions moved above main() and main guard ---

# =====================
# CONFIGURATION SECTION
# =====================
CASE_NAMES = [
    'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S2_TRENDY2025_pt1.202508011'
]
OUTPATH = 'figs/Trendy25/'
PLOT_REGION = 'Boreal'  # Options: 'Norway', 'Nordic', 'Global', 'Biased', 'Boreal'
PLOT_ILAMB = False
PLOT_DUST = False
PLOT_STRUCTURE = False
PLOT_DIM1 = False
PLOT_SEED = False
PLOT_NBP = False
PLOT_BY_PFT = True
PFT_NUMBER = 1  # zero indexed
PFT_NAMES = [
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
PROCESS_LAST_10_YEARS = False
PROCESS_FIRST_N_YEARS = False
FIRST_N_YEARS = 5
PROCESS_SELECTED_YEARS = True
SELECT_YR_RANGE = [98, 100]
CALC_ANNUAL = True

# --- Helper: Observation datasets and results dict ---
def get_obs_datasets_and_results():
    obs_dir_ilamb = '/nird/datalake/NS9560K/diagnostics/ILAMB-Data/DATA/'
    if PLOT_ILAMB:
        obs_datasets = {
            'tas': f'{obs_dir_ilamb}tas/CRU4.07-1900/tas.nc',
            'pr': f'{obs_dir_ilamb}pr/CMAPv1904/pr.nc',
            'lai': f'{obs_dir_ilamb}lai/AVHRR/lai_0.5x0.5.nc',
            'gpp': f'{obs_dir_ilamb}gpp/FLUXCOM/gpp.nc',
            'biomass': f'{obs_dir_ilamb}biomass/Saatchi2011/biomass_0.5x0.5.nc',
            'cSoil': f'{obs_dir_ilamb}cSoil/HWSD2/hwsd2_cSoil.nc',
            'hfss': f'{obs_dir_ilamb}hfss/FLUXCOM/hfss.nc',
            'hfls': f'{obs_dir_ilamb}hfls/FLUXCOM/hfls.nc',
            'swe': f'{obs_dir_ilamb}swe/CanSISE/swe.nc',
        }
    else:
        obs_datasets = {
            'AVH15C1': f'{obs_dir_ilamb}lai/AVH15C1/lai.nc',
            'AVHRR': f'{obs_dir_ilamb}lai/AVHRR/lai_0.5x0.5.nc',
            'GIMMS_LAI4g': f'{obs_dir_ilamb}lai/GIMMS_LAI4g/cao2023_lai.nc',
            'MODIS': f'{obs_dir_ilamb}lai/MODIS/lai_0.5x0.5.nc'
        }
    obs_results = {obs_name: {} for obs_name in obs_datasets.keys()}
    return obs_datasets, obs_results

# --- Helper: Locations by region ---
def get_locations():
    if PLOT_DUST:
        return {'Sahara': {'lat': 23, 'lon': -5}, 'Mid Austr.': {'lat': -25.0, 'lon': 136.0}, 'SW Austr.': {'lat': -29.0, 'lon': 127.0}, 'Turkmenist.': {'lat': 40, 'lon': 58}, 'Atacama': {'lat': -24.5, 'lon': -69.25}}
    if PLOT_REGION == 'Global':
        return {'Amazon': {'lat': -0.5, 'lon': -65}, 'C4_Gr': {'lat': 13.0, 'lon': 16.5}, 'BL_CD': {'lat': 39.0, 'lon': -80.5}, 'Boreal': {'lat': 57.5, 'lon': -121.5}, 'Cool_C3': {'lat': 53, 'lon': -7.5}, 'Larch': {'lat': 61, 'lon': 122}, 'Arc_grass': {'lat': 68, 'lon': 120.0}}
    if PLOT_REGION == 'Biased':
        return {'C4_Gr': {'lat': 13.0, 'lon': 16.5}, 'Boreal': {'lat': 57.5, 'lon': -121.5}, 'Cool_C3': {'lat': 53, 'lon': -7.5}, 'Larch': {'lat': 61, 'lon': 122}, 'GUICHOU': {'lat': 27.0, 'lon': 106.0}}
    if PLOT_REGION == 'Boreal':
        return {'Boreal Can': {'lat': 57.5, 'lon': -121.5}, 'Hurdal': {'lat': 60.4, 'lon': 11.1}, 'Hyytiälä': {'lat': 61.5, 'lon': 24.0}, 'Boreal Fin': {'lat': 64.0, 'lon': 26.0}, 'Boreal W-Rus': {'lat': 62.0, 'lon': 55.0}, 'Boreal Sib': {'lat': 58.0, 'lon': 100.0}}
    if PLOT_REGION == 'Norway':
        return {'Hurdal': {'lat': 60.4, 'lon': 11.1}, 'Iskoras': {'lat': 69.1, 'lon': 25.0}, 'Finse': {'lat': 60.6, 'lon': 7.5}}
    if PLOT_REGION == 'Nordic':
        return {'Hyytiälä': {'lat': 61.5, 'lon': 24.0}, 'Hurdal': {'lat': 60.4, 'lon': 11.1}, 'Abisko': {'lat': 68.4, 'lon': 18.8}}
    return {}

# --- Helper: Process obs datasets ---
def process_obs_datasets(obs_datasets, locations, obs_results):
    import warnings
    for obs_name, obs_file in obs_datasets.items():
        try:
            with xr.open_dataset(obs_file, engine='netcdf4') as obs_data:
                if PLOT_ILAMB:
                    obs_var = obs_name
                    obs_field = obs_data[obs_var] if obs_var in obs_data else list(obs_data.data_vars.values())[0]
                else:
                    obs_field = obs_data['lai']
                obs_lats = obs_data['lat'].values
                obs_lons = obs_data['lon'].values
                if hasattr(obs_field, 'units'):
                    if obs_name in ['hfss', 'hfls'] and obs_field.units == 'MJ m-2 day-1':
                        obs_field = obs_field * (1e6 / 86400)
                        obs_field.attrs['units'] = 'W/m2'
                    elif obs_name == 'pr' and obs_field.units == 'kg m-2 s-1':
                        obs_field = obs_field * 86400
                        obs_field.attrs['units'] = 'mm/d'
                    elif obs_name == 'tas' and obs_field.units == 'K':
                        obs_field = obs_field - 273.15
                        obs_field.attrs['units'] = 'degC'
                    elif obs_name == 'gpp' and obs_field.units == 'g m-2 day-1':
                        obs_field = obs_field * 0.001 * 365
                        obs_field.attrs['units'] = 'kg/m2/yr'
                    elif obs_name == 'swe' and obs_field.units == 'm':
                        obs_field = obs_field * 1000
                        obs_field.attrs['units'] = 'mm'
                        obs_field = obs_field.where(obs_field < 1e30, 0)
                        obs_field = obs_field.where(obs_field > 1e-5, 0)
                if 'time' in obs_field.dims:
                    obs_clim = obs_field.groupby('time.month').mean(dim='time')
                else:
                    obs_clim = obs_field
                for loc, coords in locations.items():
                    lat = coords['lat']
                    lon = coords['lon']
                    lat_diff = np.abs(obs_lats - lat)
                    lon_diff = np.abs(obs_lons - lon)
                    closest_idx = np.unravel_index(np.argmin(lat_diff[:, None] + lon_diff), (len(obs_lats), len(obs_lons)))
                    if 'time' in obs_field.dims:
                        obs_results[obs_name][loc] = obs_clim[:, closest_idx[0], closest_idx[1]].values
                    else:
                        obs_results[obs_name][loc] = obs_clim[closest_idx[0], closest_idx[1]].values
        except FileNotFoundError:
            warnings.warn(f"Obs file not found: {obs_file}")
        except Exception as e:
            warnings.warn(f"Error processing obs {obs_name}: {e}")
    return obs_results

# --- Helper: Model file processing (core data extraction) ---
def process_model_files(timeseries_files, variables, locations):
    results = {loc: {var: [] for var in variables} for loc in locations}
    units = {loc: {var: None for var in variables} for loc in locations}
    time = {loc: [] for loc in locations}
    first_file = True
    SE_grid = None
    for filename in timeseries_files:
        try:
            with xr.open_dataset(filename, engine='netcdf4') as data:
                if first_file:
                    lats = data['lat'].values
                    lons = data['lon'].values
                    area = data['area'].values
                    SE_grid = area.ndim == 1
                    for loc, coords in locations.items():
                        lat = coords['lat']
                        lon = coords['lon']
                        lat_diff = np.abs(lats - lat)
                        lon_diff = np.abs(lons - lon)
                        if SE_grid:
                            closest_idx = np.nanargmin(lat_diff + lon_diff)
                        else:
                            closest_idx = np.unravel_index(np.nanargmin(lat_diff[:, None] + lon_diff), (len(lats), len(lons)))
                        locations[loc]['closest_idx'] = closest_idx
                    first_file = False
                for loc in locations:
                    closest_idx = locations[loc]['closest_idx']
                    for var in variables:
                        if var == 'RAIN+SNOW':
                            rain_data = data['RAIN'][:, closest_idx] if SE_grid else data['RAIN'][:, closest_idx[0], closest_idx[1]]
                            snow_data = data['SNOW'][:, closest_idx] if SE_grid else data['SNOW'][:, closest_idx[0], closest_idx[1]]
                            var_data = rain_data + snow_data
                            if hasattr(rain_data, 'units') and rain_data.units == 'mm/s':
                                var_data = var_data * 86400
                                units[loc][var] = 'mm/d'
                            else:
                                units[loc][var] = getattr(rain_data, 'units', 'unknown units')
                            results[loc][var].append(var_data.values)
                            continue
                        if var in data:
                            var_data = data[var]
                            if SE_grid:
                                if var_data.ndim > 2:
                                    if PLOT_BY_PFT and var.endswith('_PF'):
                                        var_data = var_data[:, PFT_NUMBER, :]
                                    else:
                                        var_data = var_data.sum(axis=1)
                                var_data = var_data[:, closest_idx]
                            else:
                                if var_data.ndim > 3:
                                    if PLOT_BY_PFT and var.endswith('_PF'):
                                        var_data = var_data[:, PFT_NUMBER, :]
                                    else:
                                        var_data = var_data.sum(axis=1)
                                var_data = var_data[:, closest_idx[0], closest_idx[1]]
                            if hasattr(var_data, 'units'):
                                if var_data.units == 'mm/s':
                                    var_data = var_data * 86400
                                    units[loc][var] = 'mm/d'
                                elif var_data.units == 'K':
                                    var_data = var_data - 273.15
                                    units[loc][var] = 'C'
                                elif var_data.units == 'kg m-2 s-1':
                                    var_data = var_data * 86400 * 365
                                    units[loc][var] = 'kg m-2 yr-1'
                                elif var_data.units == 'gC/m^2/s':
                                    var_data = var_data * 86400 * 365 / 1000
                                    units[loc][var] = 'kgC/m2/yr'
                                elif var_data.units == 'gC/m^2':
                                    var_data = var_data / 1000
                                    units[loc][var] = 'kgC/m2'
                                else:
                                    units[loc][var] = var_data.units
                            else:
                                units[loc][var] = 'unknown units'
                            results[loc][var].append(var_data.values)
                    time[loc].append(data['time'].values)
        except FileNotFoundError:
            print(f"File not found: {filename}")
        except ValueError as e:
            print(f"Error reading the file: {e}")
    for loc in locations:
        for var in variables:
            if results[loc][var]:
                results[loc][var] = np.concatenate(results[loc][var])
            else:
                results[loc][var] = np.array([])
        if time[loc]:
            time[loc] = np.concatenate(time[loc])
        else:
            time[loc] = np.array([])
    return results, units, time

# --- Helper: Plot variables (optionally add CUE) ---
def get_plot_variables(variables, results, locations):
    plot_variables = variables.copy()
    if all(var in results[loc] for var in ["FATES_NPP", "FATES_GPP"] for loc in locations):
        if "CUE" not in plot_variables:
            plot_variables.append("CUE")
    return plot_variables

# --- Helper: Add CUE if possible ---
def add_cue_if_possible(results, locations):
    for loc in locations:
        if "FATES_NPP" in results[loc] and "FATES_GPP" in results[loc]:
            npp = results[loc]["FATES_NPP"]
            gpp = results[loc]["FATES_GPP"]
            results[loc]["CUE"] = np.divide(npp, gpp, out=np.zeros_like(npp), where=gpp != 0)
            results[loc]["CUE"] = np.clip(results[loc]["CUE"], 0, None)
    return results

# --- Helper: Calculate annual means ---
def calc_annual_means(results, time, plot_variables, locations):
    for loc in locations:
        for var in plot_variables:
            if var in results[loc] and results[loc][var].size > 0:
                results[loc][var] = np.mean(results[loc][var].reshape(-1, 12), axis=1)
        if time[loc].size > 0:
            time[loc] = time[loc][::12]
    return results, time

# --- Helper: Plotting ---
def plot_results(results, units, time, plot_variables, locations, obs_results, case_name):
    fig, axes = plt.subplots(len(plot_variables), len(locations), figsize=(30, 25), sharex='col')
    location_titles = ', '.join([f"{loc} (lat: {coords['lat']}, lon: {coords['lon']})" for loc, coords in locations.items()])
    fig_title = f"Timeseries Data at {location_titles}"
    fig.suptitle(fig_title, fontsize=20)
    sim_to_ilamb = {"TSA": "tas", "RAIN+SNOW": "pr", "TLAI": "lai", "FATES_GPP": "gpp", "FATES_VEGC": "biomass", "TOTSOMC_1m": "cSoil", "FSH": "hfss", "EFLX_LH_TOT": "hfls", "H2OSNO": "swe"}
    for i, var in enumerate(plot_variables):
        for j, loc in enumerate(locations):
            if var in results[loc] and len(results[loc][var]) > 0:
                var_data = results[loc][var]
                var_units = 'unknown' if var == 'CUE' else units[loc][var]
                axes[i, j].plot(var_data, label='Model', linewidth=2)
                if PLOT_ILAMB:
                    obs_var = sim_to_ilamb.get(var, None)
                    if obs_var and obs_var in obs_results and loc in obs_results[obs_var]:
                        obs_data = obs_results[obs_var][loc]
                        if obs_data is not None:
                            if obs_data.ndim == 1 and len(var_data) > 0:
                                if len(obs_data) < 12:
                                    axes[i, j].axhline(np.mean(obs_data), color='C1', linestyle='--', label=f'ILAMB {obs_var}', linewidth=2)
                                else:
                                    obs_repeats = int(len(var_data) / 12)
                                    obs_plot = np.tile(obs_data, obs_repeats)
                                    axes[i, j].plot(obs_plot, label=f'ILAMB {obs_var}', linestyle='--', linewidth=2)
                            elif obs_data.ndim == 0:
                                axes[i, j].axhline(obs_data, color='C1', linestyle='--', label=f'ILAMB {obs_var}', linewidth=2)
                else:
                    if var == 'TLAI':
                        for obs_name, obs_data_dict in obs_results.items():
                            if loc in obs_data_dict:
                                obs_data = obs_data_dict[loc]
                                if len(var_data) > 120 and not CALC_ANNUAL:
                                    obs_repeats = 10
                                    obs_plot = np.tile(obs_data, obs_repeats)
                                    x_obs = np.arange(len(var_data)-120, len(var_data))
                                    axes[i, j].plot(x_obs, obs_plot, label=obs_name, linestyle='--', linewidth=2)
                                elif CALC_ANNUAL:
                                    obs_repeats = len(var_data) // 12
                                    obs_plot = np.mean(obs_data)
                                    obs_plot = np.tile(obs_plot, obs_repeats)
                                    axes[i, j].plot(obs_plot, label=obs_name, linestyle='--', linewidth=2)
                                else:
                                    obs_repeats = int(len(var_data) / 12)
                                    obs_plot = np.tile(obs_data, obs_repeats)
                                    axes[i, j].plot(obs_plot, label=obs_name, linestyle='--', linewidth=2)
                axes[i, j].set_title(f"{var} ({var_units})")
                axes[i, j].grid(True)
            if i == len(plot_variables) - 1:
                axes[i, j].set_xlabel('Time')
            if j == 0:
                axes[i, j].set_ylabel('')
            if var == 'TLAI' and PLOT_REGION != 'Global':
                axes[i, j].legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    output_filename = f"{OUTPATH}/Point_Timeseries_{case_name}"
    if PLOT_REGION != 'Global':
        output_filename += f"_{PLOT_REGION}"
    if PROCESS_LAST_10_YEARS:
        output_filename += "_last10years"
    elif PROCESS_FIRST_N_YEARS:
        output_filename += f"_first{FIRST_N_YEARS}years"
    elif PROCESS_SELECTED_YEARS:
        output_filename += f"_selected{SELECT_YR_RANGE[0]}-{SELECT_YR_RANGE[1]}years"
    if PLOT_STRUCTURE:
        output_filename += "_structure"
    elif PLOT_DUST:
        output_filename += "_dust"
    elif PLOT_NBP:
        output_filename += "_NBP"
    elif PLOT_SEED:
        output_filename += "_seed"
    elif PLOT_ILAMB:
        output_filename += "_ILAMB"
    elif PLOT_DIM1:
        output_filename += "_dim1"
    else:
        output_filename += "_survival"
    output_filename += ".png"
    plt.savefig(output_filename)
    print(f"Plots saved to {output_filename}")
