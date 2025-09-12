#!/bin/bash

module load CDO/2.2.2-gompi-2023b

variables=("TOTSOMC_1m" "FATES_NEP" "TWS" "H2OSNO" "H2OSFC" "H2OCAN" "TLAI" "TOTEXICE_VOL" "TOTSOILLIQ" "TOTSOILICE" "FATES_NOCOMP_PATCHAREA_PF" "FATES_GPP_PF" "FATES_NPP_PF" "FATES_VEGC_PF" "FATES_FRACTION" "area" "landfrac" "landarea_m2" "landareafates_m2" "nbedrock")

vars_cdo="TOTSOMC_1m,FATES_NEP,TWS,H2OSNO,H2OSFC,H2OCAN,TLAI,TOTEXICE_VOL,FATES_NOCOMP_PATCHAREA_PF,FATES_GPP_PF,FATES_NPP_PF,FATES_VEGC_PF,TOTSOILLIQ,TOTSOILICE,FATES_FRACTION,area,landfrac,nbedrock"

# Conversion factors
km2_to_m2=1e6
g_to_Gt=1e-15
kg_to_Gt=1e-12
CO2_to_C=0.272915 #12.011/44.01  # Conversion factor from CO2 to C
seconds_per_year=31536000  

case=i1850.ne30pg3_tn14.fatesnocomp.ctsm5.3.045_noresm_v10.CPLHIST_noLU_coldstart_v25u.20250905
#results_path=/nird/datalake/NS9560K/noresm3/cases/${case}/lnd/hist/
results_path=/nird/datalake/NS9560K/kjetisaa/${case}/lnd/hist/
temp_path=/nird/datalake/NS9560K/kjetisaa/PostProcessed/Temp_files/
output_path=/nird/datalake/NS9560K/kjetisaa/PostProcessed/${case}/

if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
    echo "Created output directory: $output_path"
fi
if [ ! -d "$temp_path" ]; then
    mkdir -p "$temp_path"
    echo "Created temporary directory: $temp_path"
fi

cd $results_path
pwd
for file in $case.clm2.h0.*.nc; do
    echo "Processing ${file}"
    cdo selvar,${vars_cdo} ${file} ${temp_path}/Temp_${file%.nc}_selected.nc
done 

cd $temp_path
echo "Merging time dimension..."
cdo -O mergetime Temp_*_selected.nc Merged_Temp1.nc
cdo aexpr,"landarea_m2=area*landfrac*$km2_to_m2;landareafates_m2=area*landfrac*FATES_FRACTION*$km2_to_m2" Merged_Temp1.nc Merged_Temp.nc

for var in ${variables[@]}; do
    cdo selvar,${var} Merged_Temp.nc ${output_path}/${case}_${var}.nc
done

rm ${temp_path}/Temp_*_selected.nc
rm ${temp_path}/Merged_Temp*.nc
rm ${temp_path}/temp_*.nc 