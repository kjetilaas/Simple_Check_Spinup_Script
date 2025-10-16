import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import glob
import os

print('Starting CheckFates_Points (Refactored Version)')

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

user = 'kjetisaa'
case_names = [
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_secm2secy_SHORT.2025-10-07'
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_secm2zero_SHORT.2025-10-07'
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_secYandM2zero_SHORT.2025-10-07'
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_AllSecHarv2zero_SHORT.2025-10-08'
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_AllHarvestToZero_SHORT.2025-10-08'    
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_corrected.202508024'
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_loggingFix.20251001'
    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S0_TRENDY2025_pt3.202508021'
    #'i1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.PostADspinup_paramV26j_fatesparamV26i.2025-09-24'

    #'iHIST1700.f19_g17.fatesnocomp.ctsm5.3.045_noresm_v10.S3_TRENDY2025_pt3_SecAge20_SHORT.2025-10-13'

    'n1850.ne16pg3_tn14.hybrid_fates-nocomp.noresm3_0_beta03a.2025-10-10'
    #'n1850.ne16_tn14.noresm3_0_beta03_rafsipmasstonum.20250930'
    #'i1850.ne16pg3_tn14.fatesnocomp.ctsm5.3.045_noresm_v13.CPLHIST_noLU_paramV26j_fatesparamV26i.2025-10-08'
]

trendy_flag = False  # If True, use TRENDY file naming and paths
key_noresmflag = True  # If True, use noresm file naming and paths

# Land Use forcing file configuration
lu_forcing_dir = '/nird/datalake/NS9560K/kjetisaa/LU_forcing_copy/'
lu_forcing_file = 'LUH2_timeseries_to_surfdata_1.9x2.5_250723_cdf5.nc'

# Output path
if trendy_flag:
    outpath = f'/datalake/NS9560K/www/diagnostics/noresm/{user}/TRENDY25/{case_names[0]}'
elif key_noresmflag:
    outpath = f'/datalake/NS9560K/www/diagnostics/noresm/{user}/NorESM_Key_Simulations/{case_names[0]}'
else:
    outpath = f'/datalake/NS9560K/www/diagnostics/noresm/{user}/{case_names[0]}'

if not os.path.exists(outpath):
    print(f'Could not find outpath: {outpath}, use figs/ instead')
    outpath = 'figs/'

# Plotting options
plot_region = 'South_America'  # 'Norway', 'Nordic', 'Global', 'Biased', 'Boreal', 'Arctic', 'Dust', 'Spikes'
plot_varset = 'ilamb'      # 'ilamb', 'dust', 'structure', 'dim1', 'seed', 'NBP', 'PFT', 'LU', 'mortality', 'size_class'

# Time filtering options
process_last_10_years = False
process_first_n_years = False
first_n_years = 1
process_selected_years = False
select_yr_range = [0, 3]  # Assumes 12 files pr year, and where zero is first year of simulation (regarless of start year)
calc_annual = False

# PFT options
plot_by_pft = False
pft_number = 1  # zero indexed

# Size class and age bin options
plot_by_size_and_age = True  # If True, plot size class variables (_SZ) and age bin variables (_AP) as separate lines instead of summing

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def get_year_label(model_time=None):
    """Generate year label based on processing options and actual model years"""
    if model_time is not None and len(model_time) > 0:
        start_year = model_time[0].year
        end_year = model_time[-1].year
        
        if process_last_10_years:
            return f"Last10years_{start_year}-{end_year}"
        elif process_first_n_years:
            return f"First{first_n_years}years_{start_year}-{end_year}"
        elif process_selected_years:
            return f"Years_{start_year}-{end_year}"
        else:
            return f"AllYears_{start_year}-{end_year}"
    else:
        # Fallback to original labels if no model time available
        if process_last_10_years:
            return "Last10years"
        elif process_first_n_years:
            return f"First{first_n_years}years"
        elif process_selected_years:
            return f"Years_{select_yr_range[0]}-{select_yr_range[1]}"
        else:
            return "AllYears"

def get_variable_list(varset):
    """Get list of variables based on plot_varset"""
    variable_sets = {
        'ilamb': ["TLAI", "FATES_GPP", "FATES_VEGC", "TOTSOMC_1m", "TSA", "RAIN+SNOW", 
                 "FSH", "EFLX_LH_TOT", "FSR", "FSDS", "H2OSNO"],
        'dust': ["TLAI", "SOILWATER_10CM", "U10_DUST", "DSTFLXT", "DSTDEP", "FATES_CA_WEIGHTED_HEIGHT"],
        'structure': ["FATES_NPLANT_SZ", "FATES_NCOHORTS", "FATES_NPATCHES", "TLAI", "TOTECOSYSC",
                     "FATES_NONSTRUCTC", "FATES_STRUCTC", "FATES_BA_WEIGHTED_HEIGHT", "FATES_CA_WEIGHTED_HEIGHT"],
        'dim1': ["TLAI", "RAIN", "SNOW", "TBOT", "FATES_GPP", "FATES_NPP", "FATES_NEP", "BTRAN", 
                "SOILWATER_10CM", "TOTSOMC", "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_COLD_STATUS",
                "FATES_STOREC_TF", "FATES_CA_WEIGHTED_HEIGHT", "FATES_MORTALITY_CFLUX_CANOPY", 
                "FATES_MORTALITY_CFLUX_USTORY"],
        'NBP': ["TLAI", "FATES_FRACTION", "FATES_VEGC", "FATES_NPP", "HR", "FATES_NEP", "FATES_SEEDS_IN_EXTERN_EL", 
                "FATES_GRAZING", "FATES_FIRE_CLOSS", "TOT_WOODPRODC_LOSS", "TOTSOMC_1m",
               "FATES_CA_WEIGHTED_HEIGHT", "FCO2"],
        'mortality': ["TLAI", "FATES_NPLANT_PF", "FATES_MORTALITY_PF", "FATES_MORTALITY_CSTARV_SZ", "FATES_MORTALITY_LOGGING_SZ", "FATES_MORTALITY_FIRE_SZ",
               "FATES_MORTALITY_BACKGROUND_SZ", "FATES_MORTALITY_FREEZING_SZ", "FATES_MORTALITY_HYDRAULIC_SZ",
               "FATES_MORTALITY_AGESCEN_SZ", "FATES_MORTALITY_SENESCENCE_SZ" ],
        'PFT': ["TLAI", "FATES_GPP_PF", "FATES_NPP_PF", "FATES_NPLANT_PF", "FATES_VEGC_PF", "FATES_MORTALITY_PF",
               "FATES_MORTALITY_CFLUX_PF", "FATES_MORTALITY_CSTARV_CFLUX_PF", "FATES_MORTALITY_FIRE_CFLUX_PF",
               "FATES_MORTALITY_HYDRO_CFLUX_PF", "FATES_MORTALITY_BACKGROUND_SZ", "FATES_MORTALITY_AGESCEN_SZ",
               "FATES_MORTALITY_SENESCENCE_SZ", "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_STOREC_TF",
               "FATES_CA_WEIGHTED_HEIGHT"],
        'seed': ["TLAI", "FATES_NPP", "FATES_SEEDLING_POOL", "FATES_SEEDS_IN", "FATES_SEED_ALLOC",
                "FATES_SEED_BANK", "FATES_SEED_DECAY_EL", "FATES_SEED_GERM_EL", "FATES_AREA_PLANTS",
                "FATES_MORTALITY_LOGGING_SZ", "FATES_CA_WEIGHTED_HEIGHT"],
        'size_class': ["TLAI", "FATES_NPLANT_SZ", "FATES_VEGC_ABOVEGROUND_SZ", "FATES_LAI_CANOPY_SZ", "FATES_LAI_USTORY_SZ",
                      "FATES_MORTALITY_LOGGING_SZ", "FATES_BASALAREA_SZ", "FATES_BA_WEIGHTED_HEIGHT", "FATES_PATCHAREA_AP", 
                      "FATES_PRIMARY_AREA_AP", "FATES_SECONDARY_AREA_AP", "FATES_LAI_AP"],
        #'LU': ["TLAI", "FATES_VEGC", "TOT_WOODPRODC", "FATES_WOOD_PRODUCT", "FATES_HARVEST_WOODPROD_C_FLUX", "FATES_LUCHANGE_WOODPROD_C_FLUX", 
        #       "FATES_LITTER_IN", "FATES_LITTER_OUT", "FATES_DISTURBANCE_RATE_LOGGING", "FATES_MORTALITY_LOGGING_SZ","FATES_CA_WEIGHTED_HEIGHT"],
        'LU': ["FATES_HARVEST_WOODPROD_C_FLUX", "FATES_LUCHANGE_WOODPROD_C_FLUX", "TLAI", "FATES_VEGC", "FATES_DISTURBANCE_RATE_LOGGING", "FATES_MORTALITY_LOGGING_SZ"], 
        'LU_forcing': ['primf', 'secdf', 'primn', 'secdn', 'urban', 'pastr', 'range', 'c3ann', 'c4ann',
               'c3per', 'c4per', 'c3nfx', 'primf_harv', 'primn_harv', 'secmf_harv', 'secyf_harv', 'secnf_harv' ],  # LU forcing variables
        'LU_forcing_transitions': [
               # Primary forest transitions
               'primf_to_secdn', 'primf_to_urban', 'primf_to_c3ann', 'primf_to_c4ann', 'primf_to_c3per', 'primf_to_c4per', 
               'primf_to_c3nfx', 'primf_to_pastr', 'primf_to_range',
               # Primary non-forest transitions  
               'primn_to_secdf', 'primn_to_urban', 'primn_to_c3ann', 'primn_to_c4ann', 'primn_to_c3per', 'primn_to_c4per',
               'primn_to_c3nfx', 'primn_to_pastr', 'primn_to_range',
               # Secondary forest transitions
               'secdf_to_secdn', 'secdf_to_urban', 'secdf_to_c3ann', 'secdf_to_c4ann', 'secdf_to_c3per', 'secdf_to_c4per',
               'secdf_to_c3nfx', 'secdf_to_pastr', 'secdf_to_range',
               # Secondary non-forest transitions
               'secdn_to_secdf', 'secdn_to_urban', 'secdn_to_c3ann', 'secdn_to_c4ann', 'secdn_to_c3per', 'secdn_to_c4per',
               'secdn_to_c3nfx', 'secdn_to_pastr', 'secdn_to_range',
               # Urban transitions
               'urban_to_secdf', 'urban_to_secdn', 'urban_to_c3ann', 'urban_to_c4ann', 'urban_to_c3per', 'urban_to_c4per',
               'urban_to_c3nfx', 'urban_to_pastr', 'urban_to_range',
               # C3 annual crop transitions
               'c3ann_to_secdf', 'c3ann_to_secdn', 'c3ann_to_urban', 'c3ann_to_c4ann', 'c3ann_to_c3per', 'c3ann_to_c4per',
               'c3ann_to_c3nfx', 'c3ann_to_pastr', 'c3ann_to_range',
               # C4 annual crop transitions
               'c4ann_to_secdf', 'c4ann_to_secdn', 'c4ann_to_urban', 'c4ann_to_c3ann', 'c4ann_to_c3per', 'c4ann_to_c4per',
               'c4ann_to_c3nfx', 'c4ann_to_pastr', 'c4ann_to_range',
               # C3 perennial crop transitions
               'c3per_to_secdf', 'c3per_to_secdn', 'c3per_to_urban', 'c3per_to_c3ann', 'c3per_to_c4ann', 'c3per_to_c4per',
               'c3per_to_c3nfx', 'c3per_to_pastr', 'c3per_to_range',
               # C4 perennial crop transitions
               'c4per_to_secdf', 'c4per_to_secdn', 'c4per_to_urban', 'c4per_to_c3ann', 'c4per_to_c4ann', 'c4per_to_c3per',
               'c4per_to_c3nfx', 'c4per_to_pastr', 'c4per_to_range',
               # C3 nitrogen-fixing crop transitions
               'c3nfx_to_secdf', 'c3nfx_to_secdn', 'c3nfx_to_urban', 'c3nfx_to_c3ann', 'c3nfx_to_c4ann', 'c3nfx_to_c3per',
               'c3nfx_to_c4per', 'c3nfx_to_pastr', 'c3nfx_to_range',
               # Pasture transitions
               'pastr_to_secdf', 'pastr_to_secdn', 'pastr_to_urban', 'pastr_to_c3ann', 'pastr_to_c4ann', 'pastr_to_c3per',
               'pastr_to_c4per', 'pastr_to_c3nfx', 'pastr_to_range',
               # Rangeland transitions
               'range_to_secdf', 'range_to_secdn', 'range_to_urban', 'range_to_c3ann', 'range_to_c4ann', 'range_to_c3per',
               'range_to_c4per', 'range_to_c3nfx', 'range_to_pastr'
        ],
        'LU_forcing_harvest': [
               # Wood harvest
               'primf_harv', 'primn_harv', 'secmf_harv', 'secyf_harv', 'secnf_harv',
               # Biofuel harvest 
               'primf_bioh', 'primn_bioh', 'secmf_bioh', 'secyf_bioh', 'secnf_bioh',
               # Crop management (fertilizer, irrigation, burning)
               'fertl_c3ann', 'irrig_c3ann', 'crpbf_c3ann', 'fertl_c4ann', 'irrig_c4ann', 'crpbf_c4ann',
               'fertl_c3per', 'irrig_c3per', 'crpbf_c3per', 'fertl_c4per', 'irrig_c4per', 'crpbf_c4per',
               'fertl_c3nfx', 'irrig_c3nfx', 'crpbf_c3nfx',
               # Forage harvest
               'fharv_c3per', 'fharv_c4per',
               # Natural disturbances and other
               'flood', 'rndwd', 'fulwd', 'combf', 'crpbf_total'     
        ]          
    }
    
    return variable_sets.get(varset, ["TLAI", "FATES_GPP", "FATES_NPP", "BTRAN", "SOILWATER_10CM", "TOTSOMC", 
                                     "FATES_GROWTH_RESP", "FATES_MAINT_RESP", "FATES_MORTALITY_BACKGROUND_SZ", 
                                     "FATES_MORTALITY_CSTARV_SZ", "FATES_MORTALITY_FREEZING_SZ", "FATES_MORTALITY_HYDRAULIC_SZ", 
                                     "FATES_MORTALITY_IMPACT_SZ", "FATES_MORTALITY_CANOPY_SZ", "FATES_MORTALITY_USTORY_SZ"])

def get_locations_dict(region):
    """Get locations dictionary based on plot_region"""
    all_locations = {
        'Global': {
            'Amazon': {'lat': -0.5, 'lon': -65},
            'C4_Gr': {'lat': 13.0, 'lon': 16.5},
            'BL_CD': {'lat': 39.0, 'lon': -80.5},
            'Hyytiälä': {'lat': 61.5, 'lon': 24.0},
            'Cool_C3': {'lat': 53, 'lon': -7.5},
            'Larch': {'lat': 61, 'lon': 122},
            'Arc_grass': {'lat': 68, 'lon': 120.0},
        },
        'Arctic': {
            'Mid Alaska': {'lat': 64.5, 'lon': -150.0},
            'East Canada': {'lat': 56.0, 'lon': -75.0},
            'Northern Sweden': {'lat': 68.0, 'lon': 18.0},
            'AC3_YeniseiHead': {'lat': 64.5, 'lon': 95.5},
            'BDAS_eSib': {'lat': 60.5, 'lon': 162.0},
        },
        'Biased': {
            'C4_Gr': {'lat': 13.0, 'lon': 16.5},
            'Boreal': {'lat': 57.5, 'lon': -121.5},
            'Cool_C3': {'lat': 53, 'lon': -7.5},
            'Larch': {'lat': 61, 'lon': 122},
            'GUICHOU': {'lat': 27.0, 'lon': 106.0}
        },
        'Spikes': { 
            'Spike1': {'lat': 20.0, 'lon': 110.0},
            'Spike2': {'lat': 46.5, 'lon': 25.0},
            'Spike3': {'lat': 57.5, 'lon': 60.0},
            'Spike4': {'lat': 57.5, 'lon': 12.5}            
        },
        'Boreal': {
            'W-Can': {'lat': 57.5, 'lon': -121.5},
            'Hurdal': {'lat': 60.4, 'lon': 11.1},
            'Hyytiälä': {'lat': 61.5, 'lon': 24.0},
            'Cen Fin': {'lat': 64.0, 'lon': 26.0},
            'W-Rus': {'lat': 62.0, 'lon': 55.0},
            'Bor Sib': {'lat': 58.0, 'lon': 100.0}
        },
        'South_America': {
            'Amazonas': {'lat': -0.5, 'lon': -65},
            'Suriname': {'lat': 4.0, 'lon': -55.0},
            'NE_Bolivia': {'lat': -12.0, 'lon': -62.0},
            'Porto Alegre': {'lat': -30.0, 'lon': -51.0}
        },
        'Dust': {
            'Sahara': {'lat': 23, 'lon': -5},
            'Mid Austr.': {'lat': -25.0, 'lon': 136.0},
            'SW Austr.': {'lat': -29.0, 'lon': 127.0},
            'Turkmenist.': {'lat': 40, 'lon': 58},
            'Atacama': {'lat': -24.5, 'lon': -69.25}
        },
        'Norway': {
            'Bergen': {'lat': 60.4, 'lon': 5.3},
            'Hurdal': {'lat': 60.4, 'lon': 11.1},
            'Austmarka': {'lat': 60.2, 'lon': 12.5},
            'Finse': {'lat': 60.6, 'lon': 7.5},
            'Trondheim': {'lat': 63.4, 'lon': 10.4},
            'Iskoras': {'lat': 69.1, 'lon': 25.0}            
        },
        'Nordic': {
            'Bergen': {'lat': 60.4, 'lon': 5.3},
            'Hurdal': {'lat': 60.4, 'lon': 11.1},            
            'Finse': {'lat': 60.6, 'lon': 7.5},
            'Norunda': {'lat': 60.1, 'lon': 17.4},            
            'Abisko': {'lat': 68.4, 'lon': 18.8},
            'Hyytiälä': {'lat': 61.5, 'lon': 24.0}
        }
    }
    
    return all_locations.get(region, all_locations['Global'])

def convert_model_units(var_data, var_name):
    """Convert model variable units"""
    if not hasattr(var_data, 'units'):
        return var_data, 'unknown units'
        
    original_units = var_data.units

    if original_units == 'mm/s':
        var_data = var_data * 86400
        return var_data, 'mm/d'
    elif original_units == 'K':
        var_data = var_data - 273.15
        return var_data, 'C'
    elif original_units == 'kg m-2 s-1':
        var_data = var_data * 86400 * 365
        return var_data, 'kg m-2 yr-1'
    elif original_units == 'gC/m^2/s':
        var_data = var_data * 86400 * 365 / 1000
        return var_data, 'kgC/m2/yr'
    elif original_units == 'gC/m^2':
        var_data = var_data / 1000
        return var_data, 'kgC/m2'
    elif original_units == 'kgCO2/m2/s':
        var_data = var_data * 86400 * 365 / 44 * 12
        return var_data, 'kgC/m2/yr'
    else:
        return var_data, original_units

def convert_obs_units(obs_field, obs_name):
    """Convert observational data units"""
    if not hasattr(obs_field, 'units'):
        return obs_field
        
    if obs_name in ['hfss', 'hfls'] and obs_field.units == 'MJ m-2 day-1':
        obs_field = obs_field * (1e6 / 86400)
        obs_field.attrs['units'] = 'W/m2'
    elif obs_name == 'pr' and obs_field.units == 'kg m-2 s-1':
        obs_field = obs_field * 86400
        obs_field.attrs['units'] = 'mm/d'
    elif obs_name == 'tas' and obs_field.units == 'K':
        obs_field = obs_field - 273.15
        obs_field.attrs['units'] = 'degC'
    elif obs_name == 'biomass' and obs_field.units == 'Mg/ha':
        obs_field = obs_field * 0.1
        obs_field.attrs['units'] = 'kg/m2'
    elif obs_name == 'gpp' and obs_field.units == 'g m-2 day-1':
        obs_field = obs_field * 0.001 * 365
        obs_field.attrs['units'] = 'kg/m2/yr'
    elif obs_name == 'swe' and obs_field.units == 'm':
        obs_field = obs_field * 1000
        obs_field.attrs['units'] = 'mm'
        obs_field = obs_field.where(obs_field < 1e30, 0)
        obs_field = obs_field.where(obs_field > 1e-5, 0)
        
    return obs_field

def get_multiple_obs_datasets():
    """Parse observational datasets from available data and organize by variable type"""
    obs_datasets_multi = {
        'tas': [
            'tas/CRU4.07-1900/tas.nc',
            'tas/CRU4.02/tas.nc'
        ],
        'pr': [
            'pr/CMAPv1904/pr.nc',
            'pr/CLASS/pr.nc',
            'pr/GPCCv2018/pr.nc',
            'pr/GPCPv2.3/pr.nc'
        ],
        'lai': [
            'lai/MODIS/lai_0.5x0.5.nc',
            'lai/AVHRR/lai_0.5x0.5.nc',
            'lai/AVH15C1/lai.nc',
            'lai/GIMMS_LAI4g/cao2023_lai.nc'
        ],
        'gpp': [
            'gpp/FLUXCOM/gpp.nc',
            'gpp/WECANN/gpp.nc'
        ],
        'biomass': [
            'biomass/ESACCI/biomass.nc',
            'biomass/GEOCARBON/biomass.nc',
            'biomass/NBCD2000/biomass_0.5x0.5.nc',
            'biomass/Saatchi2011/biomass_0.5x0.5.nc',
            'biomass/Thurner/biomass_0.5x0.5.nc',
            'biomass/USForest/biomass_0.5x0.5.nc',
            'biomass/XuSaatchi2021/XuSaatchi.nc'
        ],
        'cSoil': [
            'cSoil/HWSD2/hwsd2_cSoil.nc',
            'cSoil/Wang2024/wang2024_cSoil.nc'
        ],
        'hfss': [
            'hfss/FLUXCOM/hfss.nc',
            'hfss/CLASS/hfss.nc',
            'hfss/WECANN/hfss.nc'
        ],
        'hfls': [
            'hfls/FLUXCOM/hfls.nc',
            'hfls/CLASS/hfls.nc',
            'hfls/WECANN/hfls.nc'
        ],
        'swe': [
            'swe/CanSISE/swe.nc'
        ],
        'rsds': [
            'rsds/CERESed4.2/rsds.nc',
            'rsds/GEWEX.SRB/rsds_0.5x0.5.nc'
        ],
        'rsus': [
            'rsus/CERESed4.2/rsus.nc',
            'rsus/GEWEX.SRB/rsus_0.5x0.5.nc'
        ]
    }
    
    return obs_datasets_multi

def get_bin_info(var_type):
    """Get bin colors and labels based on FATES size classes or age bins"""
    
    if var_type == 'size':
        # FATES size class bin edges in cm
        bins = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        unit = 'cm'
        
        # Generate distinct colors for 13 size classes
        colors1 = cm.tab10(np.linspace(0, 1, 10))  # First 10 colors from tab10
        colors2 = cm.Set3(np.linspace(0, 1, 12))   # Additional colors from Set3
        
        # Select 13 distinct colors
        bin_colors = [
            colors1[0], colors1[1], colors1[2], colors1[3], colors1[4],  # tab10: 0-4
            colors1[5], colors1[6], colors1[7], colors1[8], colors1[9],  # tab10: 5-9
            colors2[2], colors2[5], colors2[8]                           # Set3: additional 3
        ]
    elif var_type == 'age':
        # FATES age bin edges in years
        bins = [0, 1, 2, 5, 10, 20, 50]
        unit = 'yr'
        
        # Generate distinct colors for 7 age bins
        colors1 = cm.tab10(np.linspace(0, 1, 10))  # First 10 colors from tab10
        
        # Select 7 distinct colors
        bin_colors = [
            colors1[0], colors1[1], colors1[2], colors1[3], colors1[4],  # tab10: 0-4
            colors1[5], colors1[6]                                       # tab10: 5-6
        ]
    else:
        raise ValueError(f"Unknown var_type: {var_type}. Use 'size' or 'age'.")
    
    # Generate bin labels with ranges
    bin_labels = []
    for i in range(len(bins)):
        if i == 0:
            label = f"0-{bins[1]}{unit}"
        elif i == len(bins) - 1:
            label = f"{bins[i]}+{unit}"
        else:
            label = f"{bins[i]}-{bins[i+1]}{unit}"
        bin_labels.append(label)
    
    return bin_colors, bin_labels


def load_lu_forcing_data(locations):
    """Load Land Use forcing data and extract at locations"""
    lu_file_path = os.path.join(lu_forcing_dir, lu_forcing_file)
    
    if not os.path.exists(lu_file_path):
        print(f"Warning: LU forcing file not found: {lu_file_path}")
        return {}, []
        
    lu_data = {}
    lu_years = []
    try:
        print(f"Loading LU forcing data from: {lu_file_path}")
        with xr.open_dataset(lu_file_path, engine='netcdf4') as data:
            lats = data['lat'].values
            lons = data['lon'].values
            
            # Get LU years from YEAR variable
            if 'YEAR' in data:
                lu_years = data['YEAR'].values
                print(f"LU data years: {lu_years[0]}-{lu_years[-1]} ({len(lu_years)} years)")
            else:
                print("Warning: No YEAR variable found in LU data")
                return {}, []
            
            # Determine grid type
            area_var = None
            for var_name in ['area', 'primf']:  # Try to find a 2D variable to check grid type
                if var_name in data:
                    area_var = data[var_name]
                    break
            
            if area_var is None:
                print("Warning: Could not determine grid type for LU data")
                return {}, []
                
            se_grid = area_var.values.ndim == 2  # If 2D in file, it's SE grid (time, ncol)
            
            # Get LU forcing variables
            lu_variables = get_variable_list('LU_forcing')
            available_lu_vars = [var for var in lu_variables if var in data]
            
            print(f"Found LU forcing variables: {available_lu_vars}")
            
            # Initialize data structure
            for loc in locations:
                lu_data[loc] = {}
                
            # Find closest grid points and extract data
            for loc, coords in locations.items():
                lat = coords['lat']
                lon = coords['lon']
                
                # Find closest grid point
                lat_diff = np.abs(lats - lat)
                lon_diff = np.abs(lons - lon)
                
                if se_grid:
                    closest_idx = np.nanargmin(lat_diff + lon_diff)
                else:
                    closest_idx = np.unravel_index(np.nanargmin(lat_diff[:, None] + lon_diff), (len(lats), len(lons)))
                
                for var in available_lu_vars:
                    if se_grid:
                        var_data = data[var][:, closest_idx]
                    else:
                        var_data = data[var][:, closest_idx[0], closest_idx[1]]
                    
                    lu_data[loc][var] = var_data.values
                    
            print("LU forcing data loaded successfully")

    except Exception as e:
        print(f"Error loading LU forcing data: {e}")
        return {}, []

    return lu_data, lu_years

# ============================================================================
# SETUP DERIVED VARIABLES
# ============================================================================

plot_ilamb = (plot_varset == 'ilamb')
if plot_varset == 'PFT':
    plot_by_pft = True

# PFT names with 3-letter short-names
pft_names = {
    "broadleaf_evergreen_tropical_tree": "BET",
    "needleleaf_evergreen_extratrop_tree": "NET",
    "needleleaf_colddecid_extratrop_tree": "NDT",
    "broadleaf_evergreen_extratrop_tree": "BEET",
    "broadleaf_hydrodecid_tropical_tree": "BHT",
    "broadleaf_colddecid_extratrop_tree": "BDT",
    "broadleaf_evergreen_extratrop_shrub": "BEES",
    "broadleaf_hydrodecid_extratrop_shrub": "BHS",
    "broadleaf_colddecid_extratrop_shrub": "BDS",
    "broadleaf_evergreen_arctic_shrub": "BEAS",
    "broadleaf_colddecid_arctic_shrub": "BDAS",
    "arctic_c3_grass": "AC3",
    "cool_c3_grass": "C3G",
    "c4_grass": "C4G",
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

# Get variables and locations
variables = get_variable_list(plot_varset)
locations = get_locations_dict(plot_region)

# if outpath does not exist, create it    
if not os.path.exists(outpath):
    os.makedirs(outpath)

# ============================================================================
# LOAD OBSERVATIONAL DATA
# ============================================================================

obs_dir_ilamb = '/nird/datalake/NS9560K/diagnostics/ILAMB-Data/DATA/'
if plot_ilamb:
    # Get multiple observational datasets for each variable
    obs_datasets_multi = get_multiple_obs_datasets()
    
    # Create a single obs_datasets dict with unique names for each dataset
    obs_datasets = {}
    for var_name, file_list in obs_datasets_multi.items():
        for file_path in file_list:
            # Extract dataset name from file path (e.g., "CRU4.07-1900" from "tas/CRU4.07-1900/tas.nc")
            dataset_name = file_path.split('/')[1] if '/' in file_path else file_path.split('.')[0]
            unique_name = f"{var_name}_{dataset_name}"
            obs_datasets[unique_name] = f'{obs_dir_ilamb}{file_path}'
else:
    obs_datasets = {
        'AVHRR': f'{obs_dir_ilamb}lai/AVHRR/lai_0.5x0.5.nc',
        'MODIS': f'{obs_dir_ilamb}lai/MODIS/lai_0.5x0.5.nc'
    }

# Create a mapping from unique obs names back to variable names
obs_name_to_var = {}
if plot_ilamb:
    for var_name, file_list in obs_datasets_multi.items():
        for file_path in file_list:
            dataset_name = file_path.split('/')[1] if '/' in file_path else file_path.split('.')[0]
            unique_name = f"{var_name}_{dataset_name}"
            obs_name_to_var[unique_name] = var_name

obs_results = {obs_name: {} for obs_name in obs_datasets.keys()}

for obs_name, obs_file in obs_datasets.items():
    try:
        with xr.open_dataset(obs_file, engine='netcdf4') as obs_data:
            # For ILAMB, determine variable name from the obs_name
            if plot_ilamb:
                if obs_name in obs_name_to_var:
                    var_name = obs_name_to_var[obs_name]
                else:
                    var_name = obs_name.split('_')[0]  # fallback
                
                # Try different variable name variations
                possible_vars = [var_name, obs_name]
                if var_name == 'cSoil':
                    possible_vars.extend(['soilc', 'csoil', 'cSoil'])
                if var_name == 'gpp':
                    possible_vars.extend(['GPP'])
                
                obs_field = None
                for var_candidate in possible_vars:
                    if var_candidate in obs_data:
                        obs_field = obs_data[var_candidate]
                        break
                
                if obs_field is None:
                    # Take the first data variable if no exact match
                    data_vars = list(obs_data.data_vars.keys())
                    if data_vars:
                        obs_field = obs_data[data_vars[0]]
                        print(f"Warning: Using variable '{data_vars[0]}' for {obs_name}")
                    else:
                        print(f"Warning: No data variables found in {obs_name}")
                        continue
            else:
                obs_field = obs_data['lai']
                
            obs_lats = obs_data['lat'].values
            obs_lons = obs_data['lon'].values
            
            # Convert units
            if plot_ilamb and obs_name in obs_name_to_var:
                obs_field = convert_obs_units(obs_field, obs_name_to_var[obs_name])
            else:
                obs_field = convert_obs_units(obs_field, obs_name)
            
            # Calculate climatology
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


# Convert negative longitudes to positive
for loc in locations:
    if locations[loc]['lon'] < 0:
        locations[loc]['lon'] = 360 + locations[loc]['lon']

# ============================================================================
# LOAD LAND USE FORCING DATA
# ============================================================================

if plot_varset == 'LU' and 2==2:
    print("Loading Land Use forcing data...")
    lu_data, lu_years = load_lu_forcing_data(locations)
else:
    lu_data, lu_years = {}, []

# ============================================================================
# MAIN PROCESSING LOOP
# ============================================================================

for case_name in case_names:
    print(f"Processing case: {case_name}")  

    # If 'i' case, use kjetisaa, else use noresm3
    if case_name.startswith('i') and trendy_flag:
        case_dir = f'/nird/datalake/NS9560K/kjetisaa/TRENDY25/{case_name}/lnd/hist/'
    elif case_name.startswith('i') and not trendy_flag:
        case_dir = f'/nird/datalake/NS9560K/kjetisaa/{case_name}/lnd/hist/'            
    else:
        case_dir = f'/nird/datalake/NS9560K/noresm3/cases/{case_name}/lnd/hist/'

    print(f"Case directory: {case_dir}")

    # Find all timeseries files
    timeseries_files = sorted(glob.glob(f'{case_dir}/{case_name}.clm2.h0.*-*.nc'))
    print(f"Found {len(timeseries_files)} timeseries files.")
    if len(timeseries_files) == 0:
        print(f"No timeseries files found for case {case_name} in directory {case_dir}. Skipping to next case.")
        continue

    # Filter files based on time options
    if process_last_10_years:
        timeseries_files = timeseries_files[-120:]  # 10 years = 120 months
    elif process_first_n_years:
        timeseries_files = timeseries_files[:first_n_years * 12]
    elif process_selected_years:
        timeseries_files = timeseries_files[select_yr_range[0] * 12:select_yr_range[1] * 12]
    elif calc_annual:
        timeseries_files = timeseries_files[0:len(timeseries_files) - (len(timeseries_files) % 12)]
        print(f"Processing {len(timeseries_files)} files for annual means.")

    # ============================================================================
    # PROCESS MODEL DATA FILES
    # ============================================================================
    
    # Initialize dictionaries to store results and units
    results = {loc: {var: [] for var in variables} for loc in locations}
    units = {loc: {var: None for var in variables} for loc in locations}
    time = {loc: [] for loc in locations}

    # Process files
    first_file = True
    if len(timeseries_files) > 1:
        files_to_read = timeseries_files
    else:
        files_to_read = timeseries_files

    for filename in files_to_read: 
        print(f"Processing: {os.path.basename(filename)}")
        try:        
            with xr.open_dataset(filename, engine='netcdf4') as data:
                if first_file:        
                    lats = data['lat'].values
                    lons = data['lon'].values   
                    area = data['area'].values
                    first_file = False

                    # Check if grid is 1D or 2D
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
                            closest_idx = np.nanargmin(lat_diff + lon_diff)
                            locations[loc]['closest_idx'] = closest_idx
                        else:
                            closest_idx = np.unravel_index(np.nanargmin(lat_diff[:, None] + lon_diff), (len(lats), len(lons)))
                            locations[loc]['closest_idx'] = closest_idx

                        # Add dominant PFT info to locations dict for later use in title
                        if 'FATES_NOCOMP_PATCHAREA_PF' in data:                        
                            if SE_grid:
                                pct_pft = data['FATES_NOCOMP_PATCHAREA_PF'][:, :, closest_idx].values
                            else:
                                pct_pft = data['FATES_NOCOMP_PATCHAREA_PF'][:, :, closest_idx[0], closest_idx[1]].values

                            dom_idx = int(np.argmax(pct_pft[-1]))                            
                            dom_frac = float(pct_pft[-1, dom_idx])
                            dom_longname = list(pft_names.keys())[dom_idx]
                            dom_shortname = pft_names[dom_longname]
                            locations[loc]['dom_pft_short'] = dom_shortname
                            locations[loc]['dom_pft_frac'] = dom_frac
                            print(f"{loc}: Dominant PFT is: {dom_longname} ({dom_shortname}), {dom_frac*100:.2f}%")
                            #print(f"{loc}: Top 3 PFT are: {[list(pft_names.keys())[i] for i in np.argsort(pct_pft[-1])[-3:][::-1]]}") 
                            print(f"{loc}: Top 3 PFT are: {[f'{list(pft_names.keys())[i]} ({pct_pft[-1, i]*100:.2f}%)' for i in np.argsort(pct_pft[-1])[-3:][::-1]]}")

                for loc in locations:
                    closest_idx = locations[loc]['closest_idx']
                    for var in variables:
                        # Special handling for RAIN+SNOW
                        if var == 'RAIN+SNOW':
                            if 'RAIN' in data and 'SNOW' in data:
                                rain_data = data['RAIN'][:, closest_idx] if SE_grid else data['RAIN'][:, closest_idx[0], closest_idx[1]]
                                snow_data = data['SNOW'][:, closest_idx] if SE_grid else data['SNOW'][:, closest_idx[0], closest_idx[1]]
                                var_data = rain_data + snow_data
                                var_data.attrs['units'] = rain_data.units
                                var_data, var_units = convert_model_units(var_data, var)
                                results[loc][var].append(var_data.values)
                                units[loc][var] = var_units                                
                            continue
                            
                        # Check if the variable exists in the dataset
                        if var in data:  
                            var_data = data[var]

                            # Handle multi-dimensional variables
                            if SE_grid:
                                if var_data.ndim > 2:
                                    if plot_by_pft and var.endswith('_PF'):
                                        var_data = var_data[:, pft_number, :]
                                    elif plot_by_size_and_age and (var.endswith('_SZ') or var.endswith('_AP')):
                                        # Keep size class or age bin dimension for individual plotting
                                        var_data = var_data[:, :, closest_idx]  # (time, size_class/age_bin)
                                        # Store each bin separately
                                        bin_suffix = '_SC' if var.endswith('_SZ') else '_AB'
                                        for bin_idx in range(var_data.shape[1]):
                                            bin_var_name = f"{var}{bin_suffix}{bin_idx}"
                                            bin_var_data = var_data[:, bin_idx]
                                            bin_var_data, var_units = convert_model_units(bin_var_data, var)
                                            if bin_var_name not in results[loc]:
                                                results[loc][bin_var_name] = []
                                            results[loc][bin_var_name].append(bin_var_data.values)
                                            units[loc][bin_var_name] = var_units
                                        continue  # Skip normal processing for this variable
                                    else:
                                        var_data = var_data.sum(axis=1)
                                var_data = var_data[:, closest_idx]
                            else:
                                if var_data.ndim > 3:
                                    if plot_by_pft and var.endswith('_PF'):
                                        var_data = var_data[:, pft_number, :]
                                    elif plot_by_size_and_age and (var.endswith('_SZ') or var.endswith('_AP')):
                                        # Keep size class or age bin dimension for individual plotting
                                        var_data = var_data[:, :, closest_idx[0], closest_idx[1]]  # (time, size_class/age_bin)
                                        # Store each bin separately
                                        bin_suffix = '_SC' if var.endswith('_SZ') else '_AB'
                                        for bin_idx in range(var_data.shape[1]):
                                            bin_var_name = f"{var}{bin_suffix}{bin_idx}"
                                            bin_var_data = var_data[:, bin_idx]
                                            bin_var_data, var_units = convert_model_units(bin_var_data, var)
                                            if bin_var_name not in results[loc]:
                                                results[loc][bin_var_name] = []
                                            results[loc][bin_var_name].append(bin_var_data.values)
                                            units[loc][bin_var_name] = var_units
                                        continue  # Skip normal processing for this variable
                                    else:
                                        var_data = var_data.sum(axis=1)                                
                                var_data = var_data[:, closest_idx[0], closest_idx[1]]

                            # Unit conversion
                            var_data, var_units = convert_model_units(var_data, var)
                            results[loc][var].append(var_data.values)
                            units[loc][var] = var_units
                        
                    time[loc].append(data['time'].values)
        except FileNotFoundError:
            print(f"File not found: {filename}")
        except ValueError as e:
            print(f"Error reading the file: {e}")

    # Convert lists to numpy arrays
    for loc in locations:
        # Handle regular variables
        for var in variables:
            if len(results[loc][var]) > 0:
                results[loc][var] = np.concatenate(results[loc][var])
            else:
                results[loc][var] = np.array([])
        
        # Handle bin variables (size class and age bin) if they exist
        if plot_by_size_and_age:
            bin_vars = [key for key in results[loc].keys() if '_SC' in key or '_AP' in key]
            for var in bin_vars:
                if len(results[loc][var]) > 0:
                    results[loc][var] = np.concatenate(results[loc][var])
                else:
                    results[loc][var] = np.array([])

        if len(time[loc]) > 0:
            time[loc] = np.concatenate(time[loc])
        else:
            time[loc] = np.array([])

    # Generate year label using actual model time
    first_loc_time = time[list(locations.keys())[0]]
    year_label = get_year_label(first_loc_time if len(first_loc_time) > 0 else None)

    # ============================================================================
    # INTEGRATE LAND USE DATA
    # ============================================================================
    
    # Add LU data to results if available
    if lu_data: #and len(lu_years) > 0:
        print("Integrating Land Use forcing data...")
        
        # Convert model CF time to years for alignment
        model_years = []
        if len(time[list(locations.keys())[0]]) > 0:
            model_time = time[list(locations.keys())[0]]
            # Extract years from CF time objects
            model_years = np.array([t.year for t in model_time])
            print(f"Model data years: {model_years[0]}-{model_years[-1]} ({len(model_years)} time steps)")
        
        # Get LU forcing variables (separate from model variables)
        lu_forcing_vars = get_variable_list('LU_forcing')
        
        # Track alignment progress
        aligned_vars = []
        aligned_locations = 0
        
        for loc in locations:
            if loc in lu_data:
                aligned_locations += 1
                for var in lu_forcing_vars:
                    if var in lu_data[loc]:
                        lu_var_data = lu_data[loc][var]
                        
                        if len(model_years) > 0:
                            # Interpolate LU data to match model years
                            lu_aligned = np.interp(model_years, lu_years, lu_var_data)
                            
                            # Add LU data as new variables (separate from model variables)
                            results[loc][var] = lu_aligned
                            units[loc][var] = 'fraction'
                            
                            # Track which variables were successfully aligned
                            if var not in aligned_vars:
                                aligned_vars.append(var)
                        else:
                            print(f"Warning: No model time data available for {loc}")
        
        # Summary of alignment
        if aligned_vars:
            print(f"Successfully aligned {len(aligned_vars)} LU variables across {aligned_locations} locations")
            print(f"Aligned variables: {', '.join(aligned_vars[:5])}{'...' if len(aligned_vars) > 5 else ''}")
        else:
            print("No LU variables were successfully aligned")

    # ============================================================================
    # ADD DERIVED VARIABLES
    # ============================================================================
    
    # Combine model variables with LU forcing variables for plotting
    plot_variables = variables.copy()
    
    # Handle bin variables (size class and age bin) if plot_by_size_and_age is enabled
    if plot_by_size_and_age:
        # Find all bin variables and replace original variables with bin versions
        bin_variables = []
        for var in variables:
            if var.endswith('_SZ') or var.endswith('_AP'):
                # Find all bin versions of this variable
                for loc in locations:
                    if var.endswith('_SZ'):
                        bin_vars = [key for key in results[loc].keys() if key.startswith(f"{var}_SC")]
                    else:  # var.endswith('_AP')
                        bin_vars = [key for key in results[loc].keys() if key.startswith(f"{var}_AB")]
                    bin_variables.extend([bv for bv in bin_vars if bv not in bin_variables])
                # Remove original variable from plot list
                if var in plot_variables:
                    plot_variables.remove(var)
        
        # Add bin variables to plot list
        plot_variables.extend(bin_variables)
        print(f"Bin mode: plotting {len(bin_variables)} bin variables (size class and age bin)")
    
    if lu_data and len(lu_years) > 0:
        lu_forcing_vars = get_variable_list('LU_forcing')
        # Filter LU variables to only include those with meaningful values
        lu_threshold = 1e-40
        significant_lu_vars = []
        
        for var in lu_forcing_vars:
            # Check if this LU variable has any significant values across all locations
            has_significant_values = False
            for loc in locations:
                if var in results.get(loc, {}) and len(results[loc][var]) > 0:
                    if np.any(np.abs(results[loc][var]) > lu_threshold):
                        has_significant_values = True
                        break
            
            if has_significant_values:
                significant_lu_vars.append(var)
        
        print(f"Filtered LU variables: {len(significant_lu_vars)} of {len(lu_forcing_vars)} variables above threshold {lu_threshold}")
        if significant_lu_vars:
            print(f"Significant LU variables: {', '.join(significant_lu_vars[:5])}{'...' if len(significant_lu_vars) > 5 else ''}")
        
        # Add only significant LU variables to plot list
        for var in significant_lu_vars:
            if var not in plot_variables:
                plot_variables.append(var)
    
    # Add CUE if NPP and GPP are available
    if all(any(var in results[loc] for loc in locations) for var in ["FATES_NPP", "FATES_GPP"]):
        for loc in locations:
            if "FATES_NPP" in results[loc] and "FATES_GPP" in results[loc] and len(results[loc]["FATES_GPP"]) > 0:
                gpp_data = results[loc]["FATES_GPP"]
                npp_data = results[loc]["FATES_NPP"]
                results[loc]["CUE"] = np.divide(npp_data, gpp_data, 
                                              out=np.zeros_like(npp_data), where=gpp_data != 0)
                results[loc]["CUE"] = np.clip(results[loc]["CUE"], 0, None)
                units[loc]["CUE"] = 'dimensionless'
        
        if "CUE" not in plot_variables:
            plot_variables.append("CUE")

    # Add nbp as FATES_NEP+FATES_SEEDS_IN_EXTERN_EL-FATES_GRAZING-FATES_FIRE_CLOSS-TOT_WOODPRODC_LOSS if all are available
    # where all FATES_* components should be multiplied by FATES_FRACTION 
    nbp_components = ["FATES_NEP", "FATES_SEEDS_IN_EXTERN_EL", "FATES_GRAZING", "FATES_FIRE_CLOSS", "TOT_WOODPRODC_LOSS", "FATES_FRACTION"]
    if all(any(var in results[loc] for loc in locations) for var in nbp_components):
        for loc in locations:
            if all(var in results[loc] for var in nbp_components) and len(results[loc]["FATES_NEP"]) > 0:
                nep_data = results[loc]["FATES_NEP"]
                seeds_in_data = results[loc]["FATES_SEEDS_IN_EXTERN_EL"]
                grazing_data = results[loc]["FATES_GRAZING"]
                fire_closs_data = results[loc]["FATES_FIRE_CLOSS"]
                woodprod_loss_data = results[loc]["TOT_WOODPRODC_LOSS"]
                fraction_data = results[loc]["FATES_FRACTION"]
                # Apply fraction to FATES components
                nep_data *= fraction_data
                seeds_in_data *= fraction_data
                grazing_data *= fraction_data
                fire_closs_data *= fraction_data
                # Calculate NBP
                results[loc]["NBP"] = (nep_data + seeds_in_data - grazing_data - fire_closs_data - woodprod_loss_data)
                #results[loc]["NBP"] = (nep_data + seeds_in_data - grazing_data - woodprod_loss_data)
                units[loc]["NBP"] = 'kgC/m2/yr' 
        if "NBP" not in plot_variables:
            plot_variables.append("NBP")

    # Add missmatch (FCO2 + NBP) if FCO2 and NBP are available
    if all(any(var in results[loc] for loc in locations) for var in ["FCO2", "NBP"]):
        for loc in locations:
            if "FCO2" in results[loc] and "NBP" in results[loc] and len(results[loc]["FCO2"]) > 0:
                fco2_data = results[loc]["FCO2"]
                nbp_data = results[loc]["NBP"]
                results[loc]["Mismatch"] = fco2_data + nbp_data
                units[loc]["Mismatch"] = 'kgC/m2/yr'
        if "Mismatch" not in plot_variables:
            plot_variables.append("Mismatch")


    # Calculate annual means if requested
    if calc_annual:
        for loc in locations:
            # Handle regular variables
            for var in plot_variables:
                if var in results[loc] and len(results[loc][var]) > 0:
                    results[loc][var] = np.mean(results[loc][var].reshape(-1, 12), axis=1)
            
            # Handle bin variables (size class and age bin) if they exist
            if plot_by_size_and_age:
                bin_vars = [key for key in results[loc].keys() if '_SC' in key or '_AP' in key]
                for var in bin_vars:
                    if len(results[loc][var]) > 0:
                        results[loc][var] = np.mean(results[loc][var].reshape(-1, 12), axis=1)
            
            if len(time[loc]) > 0:
                time[loc] = time[loc][::12]

    # ============================================================================
    # PLOTTING
    # ============================================================================

    # Create time axis for plotting (use years from model time)
    plot_time = []
    if len(time[list(locations.keys())[0]]) > 0:
        model_time = time[list(locations.keys())[0]]
        if calc_annual:
            # For annual data, use just the years
            plot_time = np.array([t.year + (t.month - 1) / 12.0 for t in model_time])
        else:
            # For monthly data, use fractional years
            plot_time = np.array([t.year + (t.month - 1) / 12.0 for t in model_time])
        print(f"Plot time range: {plot_time[0]:.1f} - {plot_time[-1]:.1f}")

    # Group variables for plotting
    if plot_by_size_and_age:
        # Group bin variables (size class and age bin) by their base variable name
        plot_groups = {}
        regular_vars = []
        
        for var in plot_variables:
            if '_SC' in var or '_AB' in var:
                if '_SC' in var:
                    base_var = var.split('_SC')[0]
                else:  # '_AB' in var
                    base_var = var.split('_AB')[0]
                if base_var not in plot_groups:
                    plot_groups[base_var] = []
                plot_groups[base_var].append(var)
            else:
                regular_vars.append(var)
        
        # Create final plot list: regular variables + grouped bin variables
        final_plot_list = regular_vars + list(plot_groups.keys())
    else:
        final_plot_list = plot_variables
        plot_groups = {}

    fig, axes = plt.subplots(len(final_plot_list), len(locations), figsize=(5*len(locations), 25), sharex='col')

    for i, var_or_group in enumerate(final_plot_list):
        for j, loc in enumerate(locations):
            if plot_by_size_and_age and var_or_group in plot_groups:
                # Plot all bins (size class or age bin) for this variable
                def extract_bin_number(var_name):
                    """Extract bin number from variable name"""
                    if '_SC' in var_name:
                        return int(var_name.split('_SC')[1])
                    elif '_AB' in var_name:
                        return int(var_name.split('_AB')[1])
                    else:
                        return 0
                
                bin_vars = sorted(plot_groups[var_or_group], key=extract_bin_number)
                var_units = 'unknown units'
                
                # Determine if this is size class or age bin based on the first variable
                is_size_class = '_SC' in bin_vars[0]
                bin_type = 'size' if is_size_class else 'age'
                
                # Get bin colors and labels
                bin_colors, bin_labels = get_bin_info(bin_type)
                
                for k, bin_var in enumerate(bin_vars):
                    if bin_var in results[loc] and len(results[loc][bin_var]) > 0:
                        var_data = results[loc][bin_var]
                        var_units = units[loc].get(bin_var, 'unknown units')
                        
                        if is_size_class:
                            bin_num = int(bin_var.split('_SC')[1])
                        else:
                            bin_num = int(bin_var.split('_AB')[1])
                        
                        # Use custom color and bin-based label
                        color = bin_colors[bin_num] if bin_num < len(bin_colors) else 'black'
                        label = bin_labels[bin_num] if bin_num < len(bin_labels) else f'Bin{bin_num}'
                        
                        # Plot with proper time axis if available, otherwise use indices
                        if len(plot_time) == len(var_data):
                            axes[i, j].plot(plot_time, var_data, label=label, linewidth=2, color=color)
                        else:
                            axes[i, j].plot(var_data, label=label, linewidth=2, color=color)
                
                axes[i, j].set_title(f"{var_or_group} ({var_units})")
                axes[i, j].grid(True)
                
                # Show legend only for the first bin variable of each type and first location
                if j == 0:  # First location only
                    # Check if this is the first occurrence of this bin type
                    current_bin_type = 'size' if is_size_class else 'age'
                    is_first_of_type = True
                    
                    # Check if we've already shown a legend for this bin type
                    for prev_i in range(i):
                        prev_var = final_plot_list[prev_i]
                        if prev_var in plot_groups:
                            prev_bin_vars = plot_groups[prev_var]
                            if len(prev_bin_vars) > 0:
                                prev_is_size_class = '_SC' in prev_bin_vars[0]
                                prev_bin_type = 'size' if prev_is_size_class else 'age'
                                if prev_bin_type == current_bin_type:
                                    is_first_of_type = False
                                    break
                    
                    if is_first_of_type:
                        axes[i, j].legend(loc='upper left', fontsize=9, frameon=True, ncol=2)
                
            else:
                # Regular variable plotting (original code)
                if var_or_group in results[loc] and len(results[loc][var_or_group]) > 0:
                    var_data = results[loc][var_or_group]
                    var_units = units[loc].get(var_or_group, 'unknown units')
                    
                    # FIRST: Plot obs if available for this variable/location (so they appear behind model)
                    if plot_ilamb:
                        obs_var = sim_to_ilamb.get(var_or_group, None)
                        if obs_var:
                            # Find all observational datasets for this variable
                            matching_obs = [obs_name for obs_name in obs_results.keys() 
                                          if obs_name.startswith(f"{obs_var}_") and loc in obs_results[obs_name]]
                            
                            # Use different colors and line styles for multiple obs datasets
                            colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
                            linestyles = ['--', '-.', ':', '--', '-.', ':', '--', '-.', ':']
                            
                            for idx, obs_name in enumerate(matching_obs):
                                obs_data = obs_results[obs_name][loc]
                                if obs_data is not None and np.all(obs_data < 1e30):
                                    # Extract dataset name for label (remove variable prefix)
                                    dataset_name = obs_name.replace(f"{obs_var}_", "")
                                    color = colors[idx % len(colors)]
                                    linestyle = linestyles[idx % len(linestyles)]
                                    
                                    if obs_data.ndim == 1 and len(var_data) > 0:
                                        if len(obs_data) < 12:
                                            axes[i, j].axhline(np.mean(obs_data), color=color, linestyle=linestyle, 
                                                             label=f'{dataset_name}', linewidth=2)
                                        else:
                                            obs_repeats = int(len(var_data) / 12)
                                            obs_plot = np.tile(obs_data, obs_repeats)
                                            if len(plot_time) == len(obs_plot):
                                                axes[i, j].plot(plot_time, obs_plot, label=f'{dataset_name}', 
                                                               linestyle=linestyle, linewidth=2, color=color)
                                            else:
                                                axes[i, j].plot(obs_plot, label=f'{dataset_name}', 
                                                               linestyle=linestyle, linewidth=2, color=color)
                                    elif obs_data.ndim == 0:
                                        axes[i, j].axhline(obs_data, color=color, linestyle=linestyle, 
                                                         label=f'{dataset_name}', linewidth=2)
                                    
                            if obs_var == 'cSoil':
                                axes[i, j].set_ylim(0, 50)
                    else:
                        # Only plot obs for TLAI
                        if var_or_group == 'TLAI' and 2==2:
                            for obs_name, obs_data_dict in obs_results.items():
                                if loc in obs_data_dict:
                                    obs_data = obs_data_dict[loc]
                                    if not process_last_10_years and len(var_data) > 120 and not calc_annual:
                                        obs_repeats = 10
                                        obs_plot = np.tile(obs_data, obs_repeats)
                                        if len(plot_time) >= 120:
                                            x_obs = plot_time[-120:]
                                            axes[i, j].plot(x_obs, obs_plot, label=obs_name, linestyle='--', linewidth=2)
                                        else:
                                            x_obs = np.arange(len(var_data)-120, len(var_data))
                                            axes[i, j].plot(x_obs, obs_plot, label=obs_name, linestyle='--', linewidth=2)                               
                                    elif calc_annual:
                                        obs_repeats = len(var_data) // 12
                                        obs_plot = np.mean(obs_data)
                                        obs_plot = np.tile(obs_plot, obs_repeats)
                                        if len(plot_time) == len(obs_plot):
                                            axes[i, j].plot(plot_time, obs_plot, label=obs_name, linestyle='--', linewidth=2)
                                        else:
                                            axes[i, j].plot(obs_plot, label=obs_name, linestyle='--', linewidth=2)
                                    else:
                                        obs_repeats = int(len(var_data) / 12)
                                        obs_plot = np.tile(obs_data, obs_repeats)
                                        if len(plot_time) == len(obs_plot):
                                            axes[i, j].plot(plot_time, obs_plot, label=obs_name, linestyle='--', linewidth=2)
                                        else:
                                            axes[i, j].plot(obs_plot, label=obs_name, linestyle='--', linewidth=2)

                    # SECOND: Plot model data on top (so it's always visible)
                    if len(plot_time) == len(var_data):
                        axes[i, j].plot(plot_time, var_data, label='Model', linewidth=2, color='black', zorder=5)
                    else:
                        axes[i, j].plot(var_data, label='Model', linewidth=2, color='black', zorder=5)
                                        
                    axes[i, j].set_title(f"{var_or_group} ({var_units})")
                    axes[i, j].grid(True)
                    
                else:
                    print(f"Warning: No data available for variable '{var_or_group}' at location '{loc}'. Skipping plot.")
            
            # Format x-axis to show years nicely
            if len(plot_time) > 0:
                axes[i, j].set_xlim(plot_time[0], plot_time[-1])
                # Set x-axis labels to show every few years depending on data length
                if len(plot_time) > 0:
                    year_span = plot_time[-1] - plot_time[0]
                    if year_span > 50:
                        axes[i, j].set_xticks(np.arange(int(plot_time[0]), int(plot_time[-1])+1, 10))
                    elif year_span > 20:
                        axes[i, j].set_xticks(np.arange(int(plot_time[0]), int(plot_time[-1])+1, 5))
                    else:
                        axes[i, j].set_xticks(np.arange(int(plot_time[0]), int(plot_time[-1])+1, 2))
                
            if i == len(final_plot_list) - 1:
                axes[i, j].set_xlabel('Year')
            if j == 0:
                axes[i, j].set_ylabel('')

            # Add legend logic based on plot type and location
            if j == 0:  # First location only
                if plot_by_size_and_age and var_or_group in plot_groups:
                    # For bin variables, legend is handled separately above
                    pass
                elif plot_ilamb:
                    # For ILAMB plots, add legend for all variables in first location
                    axes[i, j].legend(loc='upper left', fontsize=10, frameon=True)
                elif i == 0:
                    # For non-ILAMB plots, only add legend to first subplot
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
    
    fig.suptitle(f"{case_name} ({year_label})", fontsize=22, fontweight='bold', y=0.99)
    fig.text(0.5, 0.96, location_titles, ha='center', va='top', fontsize=15)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save plot
    output_filename = f"{outpath}/Point_Timeseries_{case_name}"
    output_filename += f"_{plot_region}"
    output_filename += f'_{plot_varset}'
    if plot_by_size_and_age:
        output_filename += '_SizeAndAge'
    output_filename += f'_{year_label}'
    output_filename += ".png"

    plt.savefig(output_filename)
    print(f"Plots saved to {output_filename}")

print('Finished CheckFates_Points (Refactored Version)')