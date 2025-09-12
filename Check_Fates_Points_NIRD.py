import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

print('Starting CheckFates_Points')

user='kjetisaa'  # Change to your username
#user='rosief'  # Change to your username

case_names = [
            'i1850.ne30pg3_tn14.fatesnocomp.ctsm5.3.045_noresm_v10.CPLHIST_noLU_coldstart_v25u.20250905'
            ]

outpath = f'/datalake/NS9560K/www/diagnostics/noresm/{user}/{case_names[0]}'

#if outpath dir excist 
if not os.path.exists(outpath):
    print(f'Could not find outpath: {outpath}, use figs/ instead')
    outpath = 'figs/'  # Directory to save figures
#outpath = f'/datalake/NS9560K/www/diagnostics/noresm/kjetisaa/TRENDY25/{case_names[0]}'

# Locations
plot_region = 'Arctic'# 'Norway', 'Nordic', 'Global', 'Biased', 'Boreal', 'Arctic'

# Set plotting options
plot_varset = 'ilamb' # 'ilamb', 'dust', 'structure', 'dim1', 'seed', 'NBP', 'PFT'

# Option to process only the last 10 years, only works for monthly data
process_last_10_years = True
process_first_n_years = False
first_n_years = 1
process_selected_years = False
select_yr_range = [98, 100]

calc_annual = False

year_label = "AllYears"
if process_last_10_years:
    year_label = "Last10years"
elif process_first_n_years:
    year_label = f"First{first_n_years}years"
elif process_selected_years:
    year_label = f"Years_{select_yr_range[0]}-{select_yr_range[1]}"

if plot_varset == 'ilamb':
    plot_ilamb = True
else:
    plot_ilamb = False

plot_by_pft=False
if plot_varset == 'PFT':
    plot_by_pft=True
    pft_number=1  #zero indexed.

# PFT names with 3-letter short-names
pft_names = {
    "broadleaf_evergreen_tropical_tree":      "BET",
    "needleleaf_evergreen_extratrop_tree":    "NET",
    "needleleaf_colddecid_extratrop_tree":    "NDT",
    "broadleaf_evergreen_extratrop_tree":     "BEET",
    "broadleaf_hydrodecid_tropical_tree":     "BHT",
    "broadleaf_colddecid_extratrop_tree":     "BDT",
    "broadleaf_evergreen_extratrop_shrub":    "BEES",
    "broadleaf_hydrodecid_extratrop_shrub":   "BHS",
    "broadleaf_colddecid_extratrop_shrub":    "BDS",
    "broadleaf_evergreen_arctic_shrub":       "BEAS",
    "broadleaf_colddecid_arctic_shrub":       "BDAS",
    "arctic_c3_grass":                        "AC3",
    "cool_c3_grass":                          "C3G",
    "c4_grass":                               "C4G",
}

# Mapping from simulation variable names to ILAMB observation variable names
sim_to_ilamb = {
    "TSA": "tas",
    "RAIN+SNOW": "pr",
    "TLAI": "lai",
    "FATES_GPP": "gpp",
    "FATES_VEGC": "biomass",
    "TOTSOMC_1m": "cSoil",
    "FSH": "hfss",
    "EFLX_LH_TOT": "hfls",
    "H2OSNO": "swe",
    "FSR": "rsus",
    "FSDS": "rsds"
}

# if outpath does not exist, create it    
if not os.path.exists(outpath):
    os.makedirs(outpath)

for case_name in case_names:
    print(f"Processing case: {case_name}")  

    # If 'i' case, use kjetisaa, else use noresm3
    if case_name.startswith('i'):
        case_dir = f'/nird/datalake/NS9560K/{user}/{case_name}/lnd/hist/'
#        case_dir = f'/nird/datalake/NS9560K/kjetisaa/TRENDY25/{case_name}/lnd/hist/'
    else:
        case_dir = f'/nird/datalake/NS9560K/noresm3/cases/{case_name}/lnd/hist/'

    print(f"Case directory: {case_dir}")

    # Find all timeseries files
    timeseries_files = sorted(glob.glob(f'{case_dir}/{case_name}.clm2.h0.*-*.nc'))
    print(f"Found {len(timeseries_files)} timeseries files.")
    if len(timeseries_files) == 0:
        print(f"No timeseries files found for case {case_name} in directory {case_dir}. Skipping to next case.")
        continue

    if process_last_10_years:
        # Filter files to include only the last 10 years
        timeseries_files = timeseries_files[-120:]  # Assuming monthly data, 10 years = 120 months
    elif process_first_n_years:
        # Filter files to include only the first n years
        timeseries_files = timeseries_files[:first_n_years * 12]  # Assuming monthly data
    elif process_selected_years:
        timeseries_files = timeseries_files[select_yr_range[0] * 12:select_yr_range[1] * 12]  # Assuming monthly data
    elif calc_annual:
        # Filter files to include only full years
        timeseries_files = timeseries_files[0:len(timeseries_files) - (len(timeseries_files) % 12)]  # Ensure we have complete years
        print(f"Processing {len(timeseries_files)} files for annual means.")

    if plot_varset == 'ilamb': #,"QSOIL","QVEGT","QVEGE"
        variables = ["TLAI", "FATES_GPP", "FATES_VEGC", "TOTSOMC_1m", "TSA", "RAIN+SNOW", "FSH","QRUNOFF","EFLX_LH_TOT","FSR","FSDS", "H2OSNO"] #, "ALBEDO", "H2OSNO"
    elif plot_varset == 'dust':
        variables = ["TLAI","SOILWATER_10CM", 
        "U10_DUST", "DSTFLXT", "DSTDEP", "FATES_CA_WEIGHTED_HEIGHT"]
    elif plot_varset == 'structure':
        variables = ["FATES_NPLANT_SZ","FATES_NCOHORTS", "FATES_NPATCHES", "TLAI","TOTECOSYSC", 
        "FATES_NONSTRUCTC", "FATES_STRUCTC","FATES_BA_WEIGHTED_HEIGHT","FATES_CA_WEIGHTED_HEIGHT"]
    elif plot_varset == 'dim1':
        variables = ["TLAI", "RAIN", "SNOW", "TBOT", "FATES_GPP", "FATES_NPP","FATES_NEP", "BTRAN", "SOILWATER_10CM", "TOTSOMC", 
                    "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_COLD_STATUS", "FATES_STOREC_TF","FATES_CA_WEIGHTED_HEIGHT", 
                    "FATES_MORTALITY_CFLUX_CANOPY", "FATES_MORTALITY_CFLUX_USTORY"]                   
    elif plot_varset == 'NBP':
        variables = ["FCO2", "TLAI", "FATES_GPP", "FATES_NPP", "HR", "FATES_NPLANT_PF", "FATES_MORTALITY_PF", 
        "FATES_MORTALITY_CSTARV_SZ", 
        "FATES_MORTALITY_BACKGROUND_SZ", "FATES_MORTALITY_FREEZING_SZ", "FATES_MORTALITY_HYDRAULIC_SZ", "FATES_MORTALITY_AGESCEN_SZ", "FATES_MORTALITY_SENESCENCE_SZ",
        "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_STOREC_TF", "FATES_NONSTRUCTC", "FATES_STRUCTC", "TOTSOMC_1m", "FATES_CA_WEIGHTED_HEIGHT"]
    elif plot_varset == 'PFT':
        variables = ["TLAI", "FATES_GPP_PF", "FATES_NPP_PF", "FATES_NPLANT_PF", "FATES_VEGC_PF", "FATES_MORTALITY_PF", 
        "FATES_MORTALITY_CFLUX_PF", "FATES_MORTALITY_CSTARV_CFLUX_PF", "FATES_MORTALITY_FIRE_CFLUX_PF",
        "FATES_MORTALITY_HYDRO_CFLUX_PF", "FATES_MORTALITY_BACKGROUND_SZ", "FATES_MORTALITY_AGESCEN_SZ", "FATES_MORTALITY_SENESCENCE_SZ",
        "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_STOREC_TF", "FATES_CA_WEIGHTED_HEIGHT"]
    elif plot_varset == 'seed':
        variables = ["TLAI", "FATES_NPP", "FATES_SEEDLING_POOL", "FATES_SEEDS_IN", 
        "FATES_SEED_ALLOC", "FATES_SEED_BANK", "FATES_SEED_DECAY_EL",
        "FATES_SEED_GERM_EL", "FATES_AREA_PLANTS", "FATES_CA_WEIGHTED_HEIGHT"]
    else:
        variables = ["TLAI", "FATES_GPP", "FATES_NPP", "BTRAN", "SOILWATER_10CM", "TOTSOMC", "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_MORTALITY_BACKGROUND_SZ", 
        "FATES_MORTALITY_CSTARV_SZ", "FATES_MORTALITY_FREEZING_SZ", "FATES_MORTALITY_HYDRAULIC_SZ", 
        "FATES_MORTALITY_IMPACT_SZ", "FATES_MORTALITY_CANOPY_SZ", "FATES_MORTALITY_USTORY_SZ"] #"FATES_MORTALITY_FIRE_SZ"

    # --- Observational datasets setup (after plot_ilamb and locations are defined) ---
    obs_dir_ilamb = '/nird/datalake/NS9560K/diagnostics/ILAMB-Data/DATA/'
    if plot_ilamb:
        obs_datasets = {
            'tas': f'{obs_dir_ilamb}tas/CRU4.07-1900/tas.nc',
            'pr': f'{obs_dir_ilamb}pr/CMAPv1904/pr.nc',                
            #'lai': f'{obs_dir_ilamb}lai/AVHRR/lai_0.5x0.5.nc',
            'lai': f'{obs_dir_ilamb}lai/MODIS/lai_0.5x0.5.nc',
            'gpp': f'{obs_dir_ilamb}gpp/FLUXCOM/gpp.nc',
            #'biomass': f'{obs_dir_ilamb}biomass/Saatchi2011/biomass_0.5x0.5.nc',
            'biomass': f'{obs_dir_ilamb}biomass/ESACCI/biomass.nc',
            'cSoil': f'{obs_dir_ilamb}cSoil/HWSD2/hwsd2_cSoil.nc',
            'hfss': f'{obs_dir_ilamb}hfss/FLUXCOM/hfss.nc',
            'hfls': f'{obs_dir_ilamb}hfls/FLUXCOM/hfls.nc',
            'swe': f'{obs_dir_ilamb}swe/CanSISE/swe.nc',
            'rsds': f'{obs_dir_ilamb}rsds/CERESed4.2/rsds.nc',
            'rsus': f'{obs_dir_ilamb}rsus/CERESed4.2/rsus.nc'
        }
    else:
        obs_datasets = {
            #'AVH15C1': f'{obs_dir_ilamb}lai/AVH15C1/lai.nc',
            'AVHRR': f'{obs_dir_ilamb}lai/AVHRR/lai_0.5x0.5.nc',
            #'GIMMS_LAI4g': f'{obs_dir_ilamb}lai/GIMMS_LAI4g/cao2023_lai.nc',
            'MODIS': f'{obs_dir_ilamb}lai/MODIS/lai_0.5x0.5.nc'
        }

    obs_results = {obs_name: {} for obs_name in obs_datasets.keys()}


    # Initialize results dictionary for observations
    if plot_region == 'Global':
        locations = {
            'Amazon': {'lat': -0.5, 'lon': -65},  # Tropical rainforest
            'C4_Gr': {'lat': 13.0, 'lon': 16.5},  # C4 grassland in Africa
            'BL_CD': {'lat': 39.0, 'lon': -80.5},  # Broadleaf temperate forest in Eastern US
            'Boreal': {'lat': 57.5, 'lon': -121.5},  # Boreal forest
            'Cool_C3': {'lat': 53, 'lon': -7.5},  # Cool C3 grassland
            'Larch': {'lat': 61, 'lon': 122},  # Larch forest
            'Arc_grass': {'lat': 68, 'lon': 120.0},  # Siberian grassland
        }
    elif plot_region == 'Arctic': #Mid Alaska, East Canada, Northern Sweden, West Siberia, East Siberia
        locations = {
            'Mid Alaska': {'lat': 64.5, 'lon': -150.0},  # Mid Alaska
            'East Canada': {'lat': 56.0, 'lon': -75.0},  # East Canada
            'Northern Sweden': {'lat': 68.0, 'lon': 18.0},  # Northern Sweden
            'West Siberia': {'lat': 60.0, 'lon': 80.0},  # West Siberia
            'East Siberia': {'lat': 65.0, 'lon': 100.0},  # East Siberia
        }
    elif plot_region == 'Biased':    
        locations = {        
            'C4_Gr': {'lat': 13.0, 'lon': 16.5},  # C4 grassland in Africa        
            'Boreal': {'lat': 57.5, 'lon': -121.5},  # Boreal forest
            'Cool_C3': {'lat': 53, 'lon': -7.5},  # Cool C3 grassland
            'Larch': {'lat': 61, 'lon': 122},  # Larch forest        
            'GUICHOU': {'lat': 27.0, 'lon': 106.0}  # Guichou, China
        }
    elif plot_region == 'Boreal':
        locations = {
            'W-Can': {'lat': 57.5, 'lon': -121.5},  # Boreal forest in Canada
            'Hurdal': {'lat': 60.4, 'lon': 11.1},  # Hurdal, Norway
            'Hyytiälä': {'lat': 61.5, 'lon': 24.0},  # Hyytiälä, Finland                        
            'Cen Fin': {'lat': 64.0, 'lon': 26.0},  # Boreal forest in Finland
            'W-Rus': {'lat': 62.0, 'lon': 55.0},  # Boreal forest in Western Russia
            'Bor Sib': {'lat': 58.0, 'lon': 100.0}  # Boreal forest in Siberia
        }
    elif plot_region == 'Dust':    
        locations = {        
            'Sahara': {'lat': 23, 'lon': -5}, 
            'Mid Austr.': {'lat': -25.0, 'lon': 136.0}, 
            'SW Austr.': {'lat': -29.0, 'lon': 127.0}, 
            'Turkmenist.': {'lat': 40, 'lon': 58}, 
            'Atacama': {'lat': -24.5, 'lon': -69.25} 
        }    
    elif plot_region == 'Norway':
        # Plot locations in Norway   
        locations = {
            'Hurdal': {'lat': 60.4, 'lon': 11.1},  # Hurdal, Norway        
            'Iskoras': {'lat': 69.1, 'lon': 25.0},  # Iskoras, Norway
            'Finse': {'lat': 60.6, 'lon': 7.5}  # Finse, Norway
        }
    elif plot_region == 'Nordic':
        locations = {
            'Hyytiälä': {'lat': 61.5, 'lon': 24.0},  # Hyytiälä, Finland
            'Hurdal': {'lat': 60.4, 'lon': 11.1},  # Hurdal, Norway
            'Abisko': {'lat': 68.4, 'lon': 18.8}  # Abisko, Sweden
        }

    for obs_name, obs_file in obs_datasets.items():
        try:
            with xr.open_dataset(obs_file, engine='netcdf4') as obs_data:
                # For ILAMB, variable name matches obs_name; for LAI, use 'lai'
                if plot_ilamb:
                    obs_var = obs_name
                    if obs_var in obs_data:
                        obs_field = obs_data[obs_var]
                    else:
                        obs_field = list(obs_data.data_vars.values())[0]
                else:
                    obs_field = obs_data['lai']
                obs_lats = obs_data['lat'].values
                obs_lons = obs_data['lon'].values
                #obs_lons = np.where(obs_lons < 0, obs_lons + 360, obs_lons)
                # --- Unit conversions for obs ---
                if hasattr(obs_field, 'units'):
                    if obs_name in ['hfss', 'hfls'] and obs_field.units == 'MJ m-2 day-1':
                        # Convert MJ m-2 day-1 to W/m2: 1 MJ m-2 day-1 = 1e6 J / 86400 s = 11.574 W/m2
                        obs_field = obs_field * (1e6 / 86400)
                        obs_field.attrs['units'] = 'W/m2'
                    elif obs_name == 'pr' and obs_field.units == 'kg m-2 s-1':
                        # Convert kg m-2 s-1 to mm/d: 1 kg m-2 s-1 = 86400 mm/d
                        obs_field = obs_field * 86400
                        obs_field.attrs['units'] = 'mm/d'
                    elif obs_name == 'tas' and obs_field.units == 'K':
                        # Convert K to deg C
                        obs_field = obs_field - 273.15
                        obs_field.attrs['units'] = 'degC'
                    elif obs_name == 'biomass' and obs_field.units == 'Mg/ha':
                        # Convert from Mg/ha to kg/m2: 1 Mg/ha = 0.1 kg/m2
                        obs_field = obs_field * 0.1
                        obs_field.attrs['units'] = 'kg/m2'
                    elif obs_name == 'gpp' and obs_field.units == 'g m-2 day-1':
                        # Convert g m-2 day-1 to kg m-2 yr-1: 1 g m-2 day-1 = 0.001 kg * 365
                        obs_field = obs_field * 0.001 * 365
                        obs_field.attrs['units'] = 'kg/m2/yr'    
                    elif obs_name == 'swe' and obs_field.units == 'm':
                        # Convert m to mm
                        obs_field = obs_field * 1000
                        obs_field.attrs['units'] = 'mm'
                        #set missing values (>1e30) to zero 
                        obs_field = obs_field.where(obs_field < 1e30, 0)
                        #set tiny values (< 1e-5) to zero
                        obs_field = obs_field.where(obs_field > 1e-5, 0)

                # --- end unit conversions ---
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
            print(f"Obs file not found: {obs_file}")
        except Exception as e:
            print(f"Error processing obs {obs_name}: {e}")


    #convert negative longitudes to positive
    for loc in locations:
        if locations[loc]['lon'] < 0:
            locations[loc]['lon'] = 360 + locations[loc]['lon']

    # Initialize dictionaries to store results and units
    results = {loc: {var: [] for var in variables} for loc in locations}
    units = {loc: {var: None for var in variables} for loc in locations}  # New dictionary for units

    # Initialize a dictionary to store the accumulated results for each location
    time = {loc: [] for loc in locations}

    # Loop through all timeseries files except the last one
    first_file = True
    if len(timeseries_files) > 1:
        files_to_read = timeseries_files#[:-1]
    else:
        files_to_read = timeseries_files

    for filename in files_to_read: 
        print(filename)   
        try:        
            with xr.open_dataset(filename, engine='netcdf4') as data:
                if first_file:        
                    #print('Reading latitude and longitude')
                    lats = data['lat'].values
                    lons = data['lon'].values   
                    area = data['area'].values #hacky way to get dimensions
                    first_file = False

                    #Check if grid is 1D or 2D
                    if area.ndim == 1:
                        SE_grid = True
                        print('SE grid')
                    elif area.ndim == 2:
                        SE_grid = False
                        print('2D grid')
                    

                    for loc, coords in locations.items():
                        lat = coords['lat']
                        lon = coords['lon']
                        # Find the absolute difference between the given lat/lon and the lat/lon values in the dataset
                        lat_diff = np.abs(lats - lat)
                        lon_diff = np.abs(lons - lon)

                        if SE_grid:
                            # Find the index of the minimum difference for lat and lon, ignoring NaNs
                            closest_idx = np.nanargmin(lat_diff + lon_diff)
                            locations[loc]['closest_idx'] = closest_idx
                        else:
                            # Find the x,y indices of the closest grid cell
                            closest_idx = np.unravel_index(np.nanargmin(lat_diff[:, None] + lon_diff), (len(lats), len(lons)))
                            locations[loc]['closest_idx'] = closest_idx

                        # Add dominant PFT info to locations dict for later use in title
                        if 'FATES_NOCOMP_PATCHAREA_PF' in data:                        
                            if SE_grid:
                                pct_pft= data['FATES_NOCOMP_PATCHAREA_PF'][:,:, closest_idx].values
                            else:
                                pct_pft= data['FATES_NOCOMP_PATCHAREA_PF'][:, :, closest_idx[0], closest_idx[1]].values

                            dom_idx = int(np.argmax(pct_pft[-1]))
                            dom_frac = float(pct_pft[-1, dom_idx])
                            dom_longname = list(pft_names.keys())[dom_idx]
                            dom_shortname = pft_names[dom_longname]
                            locations[loc]['dom_pft_short'] = dom_shortname
                            locations[loc]['dom_pft_frac'] = dom_frac
                            # Print for debug
                            print(f"{loc}: Dominant PFT is: {dom_longname} ({dom_shortname}), {dom_frac*100:.2f} %")


                for loc in locations:
                    closest_idx = locations[loc]['closest_idx']
                    for var in variables:
                        # Special handling for RAIN+SNOW
                        if var == 'RAIN+SNOW':
                            rain_data = data['RAIN'][:, closest_idx] if SE_grid else data['RAIN'][:, closest_idx[0], closest_idx[1]]
                            snow_data = data['SNOW'][:, closest_idx] if SE_grid else data['SNOW'][:, closest_idx[0], closest_idx[1]]
                            var_data = rain_data + snow_data
                            # Unit conversion (assume same units for both)
                            if hasattr(rain_data, 'units') and rain_data.units == 'mm/s':
                                var_data = var_data * 86400  # Convert mm/s to mm/d
                                units[loc][var] = 'mm/d'
                            else:
                                units[loc][var] = getattr(rain_data, 'units', 'unknown units')
                            results[loc][var].append(var_data.values)
                            continue
                        # Check if the variable exists in the dataset
                        if var in data:  
                            var_data = data[var]

                            #TODO: rewrite this part, to not hardcode the dimensions to sum over!
                            if SE_grid:
                                # Check if the variable has more than 2 dimensions
                                if var_data.ndim > 2:
                                    if plot_by_pft and var.endswith('_PF'):
                                        var_data = var_data[:, pft_number, :]
                                    else:
                                        # Sum over the second dimension                            
                                        var_data = var_data.sum(axis=1)
                                var_data = var_data[:, closest_idx]
                            else:
                                # Check if the variable has more than 3 dimensions
                                if var_data.ndim > 3:
                                    if plot_by_pft and var.endswith('_PF'):
                                        var_data = var_data[:, pft_number, :]
                                    else:
                                        # Sum over the third dimension                            
                                        var_data = var_data.sum(axis=1)                                
                                var_data = var_data[:, closest_idx[0], closest_idx[1]]

                            # Unit conversion
                            if hasattr(var_data, 'units'):
                                if var_data.units == 'mm/s':
                                    var_data = var_data * 86400  # Convert mm/s to mm/d
                                    units[loc][var] = 'mm/d'
                                elif var_data.units == 'K':
                                    var_data = var_data - 273.15  # Convert K to Celsius
                                    units[loc][var] = 'C'
                                elif var_data.units == 'kg m-2 s-1':
                                    var_data = var_data * 86400 * 365  # Convert kg m-2 s-1 to kg m-2 yr-1
                                    units[loc][var] = 'kg m-2 yr-1'
                                elif var_data.units == 'gC/m^2/s':
                                    var_data = var_data * 86400 * 365 / 1000  # Convert gC/m2/s to kgC/m2/yr
                                    units[loc][var] = 'kgC/m2/yr'
                                elif var_data.units == 'gC/m^2':
                                    var_data = var_data / 1000  # Convert gC/m2 to kgC/m2
                                    units[loc][var] = 'kgC/m2'
                                else:
                                    units[loc][var] = var_data.units  # Keep original units if no conversion is applied
                            else:
                                units[loc][var] = 'unknown units'

                            results[loc][var].append(var_data.values)  # Store only the data
                        #else:
                            #print(f"Warning: Variable '{var}' not found in file {filename}. Skipping.")
                        
                    time[loc].append(data['time'].values)
        except FileNotFoundError:
            print(f"File not found: {filename}")
        except ValueError as e:
            print(f"Error reading the file: {e}")

    # Convert lists to numpy arrays
    for loc in locations:
        for var in variables:
            if results[loc][var]:  # Check if the list is not empty
                results[loc][var] = np.concatenate(results[loc][var])
            else:
                results[loc][var] = np.array([])

        if time[loc]:  # Check if the list is not empty
            time[loc] = np.concatenate(time[loc])
        else:
            time[loc] = np.array([])

    plot_reduced_list = False
    if plot_reduced_list:
        #Collect variables to total carbon pools
        fates_litter_vars = ["FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"]
        #fates_veg_vars = ["FATES_STRUCTC", "FATES_NONSTRUCTC", "FATES_STOREC", "FATES_REPROC"]
        for loc in locations:   
            if all(var in results[loc] for var in fates_litter_vars):
                results[loc]["FATES_LITTERC"] = results[loc]["FATES_LITTER_CWD_ELDC"] + results[loc]["FATES_LITTER_AG_CWD_EL"] + results[loc]["FATES_LITTER_AG_FINE_EL"] + results[loc]["FATES_LITTER_BG_CWD_EL"] + results[loc]["FATES_LITTER_BG_FINE_EL"]            
        #    if all(var in results[loc] for var in fates_veg_vars):
        #        results[loc]["FATES_TOTVEGC"] = results[loc]["FATES_STRUCTC"] + results[loc]["FATES_NONSTRUCTC"] + results[loc]["FATES_STOREC"] + results[loc]["FATES_REPROC"]  
        #plot_variables = [var for var in variables if var not in ["FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL", "FATES_LEAFC", "FATES_FROOTC", "FATES_STOREC", "FATES_SAPWOODC", "FATES_REPROC"]] + ["FATES_LITTERC", "FATES_TOTVEGC"]
        plot_variables = [var for var in variables if var not in ["FATES_LITTER_CWD_ELDC", "FATES_LITTER_AG_CWD_EL", "FATES_LITTER_AG_FINE_EL", "FATES_LITTER_BG_CWD_EL", "FATES_LITTER_BG_FINE_EL"]] + ["FATES_LITTERC"]
    else:
        plot_variables = variables

    # Add CUE if NPP and GPP are available
    if all(var in results[loc] for var in ["FATES_NPP", "FATES_GPP"]) :
        for loc in locations:
            # Check for division by zero
            results[loc]["CUE"] = np.divide(results[loc]["FATES_NPP"], results[loc]["FATES_GPP"], out=np.zeros_like(results[loc]["FATES_NPP"]), where=results[loc]["FATES_GPP"] != 0)
            # Set negative values to zero
            results[loc]["CUE"] = np.clip(results[loc]["CUE"], 0, None)
        # Append "CUE" to plot_variables only once
        if "CUE" not in plot_variables:
            plot_variables.append("CUE")

    if calc_annual:
        # Calculate annual means
        for loc in locations:
            for var in plot_variables:
                if var in results[loc]:
                    results[loc][var] = np.mean(results[loc][var].reshape(-1, 12), axis=1)
                    time[loc] = time[loc][::12]  # Adjust time to match the annual data

    #----- Plotting -----
    fig, axes = plt.subplots(len(plot_variables), len(locations), figsize=(5*len(locations), 25), sharex='col')

    # Update the plotting section to use the updated units
    for i, var in enumerate(plot_variables):
        for j, loc in enumerate(locations):
            if var in results[loc] and len(results[loc][var]) > 0:  # Check if the variable has data
                var_data = results[loc][var]  # Extract data
                if var == 'CUE':
                    var_units = 'unknown'  # Default if no units found
                else:
                    var_units = units[loc][var]  # Extract units
                axes[i, j].plot(var_data, label='Model', linewidth=2)

                # Plot obs if available for this variable/location
                if plot_ilamb:
                    # Use mapping to get correct obs variable name
                    obs_var = sim_to_ilamb.get(var, None)
                    if obs_var and obs_var in obs_results and loc in obs_results[obs_var]:
                        obs_data = obs_results[obs_var][loc]
                        if obs_data is not None and np.all(obs_data < 1e30):  # Check that obs values not > 1e30
                            if obs_data.ndim == 1 and len(var_data) > 0:
                                if len(obs_data) < 12:
                                    # Plot as a straight line if not monthly
                                    axes[i, j].axhline(np.mean(obs_data), color='C1', linestyle='--', label=f'ILAMB {obs_var}', linewidth=2)
                                else:
                                    obs_repeats = int(len(var_data) / 12)
                                    obs_plot = np.tile(obs_data, obs_repeats)
                                    axes[i, j].plot(obs_plot, label=f'ILAMB {obs_var}', linestyle='--', linewidth=2)
                            elif obs_data.ndim == 0:
                                axes[i, j].axhline(obs_data, color='C1', linestyle='--', label=f'ILAMB {obs_var}', linewidth=2)
                            if obs_var == 'cSoil':
                                axes[i, j].set_ylim(0, 50) #Hardcode ylim for cSoil
                else:
                    # Only plot obs for TLAI
                    if var == 'TLAI':
                        for obs_name, obs_data_dict in obs_results.items():
                            if loc in obs_data_dict:
                                obs_data = obs_data_dict[loc]
                                if not process_last_10_years and len(var_data) > 120 and not calc_annual:
                                    obs_repeats = 10
                                    obs_plot = np.tile(obs_data, obs_repeats)
                                    x_obs = np.arange(len(var_data)-120, len(var_data))
                                    axes[i, j].plot(x_obs, obs_plot, label=obs_name, linestyle='--', linewidth=2)                                
                                elif calc_annual:
                                    obs_repeats = len(var_data) // 12
                                    obs_plot = np.mean(obs_data)
                                    obs_plot = np.tile(obs_plot, obs_repeats)
                                    axes[i, j].plot(obs_plot, label=obs_name, linestyle='--', linewidth=2)
                                else:
                                    obs_repeats = int(len(var_data) / 12)
                                    obs_plot = np.tile(obs_data, obs_repeats)
                                    axes[i, j].plot(obs_plot, label=obs_name, linestyle='--', linewidth=2)
                axes[i, j].set_title(f"{var} ({var_units})")  # Use updated units
                axes[i, j].grid(True)
            else:
                print(f"Warning: No data available for variable '{var}' at location '{loc}'. Skipping plot.")
            if i == len(plot_variables) - 1:
                axes[i, j].set_xlabel('Time')
            if j == 0:
                axes[i, j].set_ylabel('')

            # Only add legend for the first subplot (i==0, j==0), and move it inside
            if i == 0 and j == 0:
                axes[i, j].legend(loc='upper left', fontsize=12, frameon=True)


    # Create the suptitle with location names, coordinates, and dominant PFT info
    def loc_title(loc, coords):
        dom_pft = coords.get('dom_pft_short', '?')
        dom_frac = coords.get('dom_pft_frac', None)
        if dom_frac is not None:
            dom_str = f", {dom_pft} {dom_frac*100:2.0f}%"
        else:
            dom_str = ''
        return f"{loc} (lat: {coords['lat']}, lon: {coords['lon']}{dom_str})"

    location_titles = ';   '.join([loc_title(loc, coords) for loc, coords in locations.items()])
    if process_last_10_years:
        fig_title = f"{location_titles} (last 10 years)"
    elif process_first_n_years:
        fig_title = f"{location_titles} (first {first_n_years} years)"
    else:
        fig_title = f"{location_titles}"

    # Add case name as a bold suptitle line above the main title
    fig.suptitle(f"{case_name} ({year_label})", fontsize=22, fontweight='bold', y=0.99)
    # Add the original title below the case name
    fig.text(0.5, 0.96, fig_title, ha='center', va='top', fontsize=15)

    plt.tight_layout(rect=[0, 0, 1, 0.96])


    output_filename = f"{outpath}/Point_Timeseries_{case_name}"
    output_filename += f"_{plot_region}"
    output_filename += f'_{plot_varset}'
    output_filename += f'_{year_label}'
    output_filename += ".png"

    plt.savefig(output_filename)
    print(f"Plots saved to {output_filename}")

    print('Finished CheckFates_Points')