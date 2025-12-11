#!/bin/bash

#Simple wrapper to run different diagnostics for NorESM. This script assumes that you run on ipcc.nird.sigma2.no,
#which also requires a minor change to the xesmf diagnostics scripts ("Marit's diagnostics"): 
#https://github.com/NorESMhub/xesmf_clm_fates_diagnostic/pull/49

do_oceanT=false
do_noresm_diagComp=false  #Requires logint to ipcc-node: ipcc.nird.sigma2.no
do_noresm_diagObs=false #Requires logint to ipcc-node: ipcc.nird.sigma2.no
do_xesm_diag=true
do_point_diag=false

#Load python
. /nird/home/kjetisaa/xesmf_clm_fates_diagnostic/scripts/setup.sh

case1=n1850.ne16pg3_tn14.noresm3_0_beta06.PPEbasedV2.2025-12-08
case2=n1850.ne16pg3_tn14.noresm3_0_beta06.Best_c8_4p7.2025-12-01
label1=b06_PPEtunv2
case_path=/nird/datalake/NS9560K/noresm3/cases/
outpath=/nird/datalake/NS9560K/www/diagnostics/noresm/kjetisaa/NorESM_Key_Simulations/

yr_start1=471
yr_end1=490
yr_start2=471
yr_end2=490

#Make output directory if it does not exist
if [ ! -d "$outpath/${case1}" ]; then
    mkdir -p "$outpath/${case1}"
    echo "Created output directory: $outpath/${case1}"
fi

if [ "$do_oceanT" = true ]; then
    echo "Plotting ocean volume temperatures for case: $case1 "
    echo "Input to Plot_ocean_volume_temp.py: $case1 $label1 $outpath/${case1}"
    python /nird/home/kjetisaa/Simple_Check_Spinup_Script/Plot_ocean_volume_temp.py $case1 $label1 $outpath/${case1}
fi

if [ "$do_noresm_diagComp" = true ]; then
    echo "Running NoRESM diagnostics for $case1 vs $case2"
    /diagnostics/noresm-diagnostics-src/bin/diag_run -m blom,hamocc,cice -c1 $case1 -i1 $case_path -s1 $yr_start1 -e1 $yr_end1 -c2 $case2 -i2 $case_path -s2 $yr_start2 -e2 $yr_end2 --type=climo -o /scratch/kjetisaa/diagnostics/noresm/out/ --web-dir=$outpath
fi

if [ "$do_noresm_diagObs" = true ]; then
    echo "Running NoRESM diagnostics for $case1 vs $case2"
    /diagnostics/noresm-diagnostics-src/bin/diag_run -m blom,hamocc,cice -c1 $case1 -i1 $case_path -s1 $yr_start1 -e1 $yr_end1 --type=climo -o /scratch/kjetisaa/diagnostics/noresm/out/ --web-dir=$outpath
fi

if [ "$do_xesm_diag" = true ]; then
    echo "Running XESM diagnostics for $case1 vs $case2"
    # Note: weight file is for remapping from ne16pg3 to 1.9x2.5 grid. Check that case1 contains "ne16pg3"
    if [[ $case1 == *"ne16pg3"* ]]; then
        echo "outpath=$outpath"
        cd /nird/home/kjetisaa/xesmf_clm_fates_diagnostic/scripts/
        echo "python /nird/home/kjetisaa/xesmf_clm_fates_diagnostic/scripts/run_diagnostic_full_from_terminal.py $case_path$case1/lnd/hist compare=$case_path$case2/lnd/hist/ weight=/nird/datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne16pg3_to_1.9x2.5_nomask_scripgrids_c250425.nc outpath=$outpath pamfile=short_pams.json compare_from_end=20"        
        python /nird/home/kjetisaa/xesmf_clm_fates_diagnostic/scripts/run_diagnostic_full_from_terminal.py $case_path$case1/lnd/hist compare=$case_path$case2/lnd/hist/ weight=/nird/datalake/NS9560K/diagnostics/land_xesmf_diag_data/map_ne16pg3_to_1.9x2.5_nomask_scripgrids_c250425.nc outpath=$outpath pamfile=short_pams.json compare_from_end=20
        cd -
    else
        echo "Error: case1 does not contain 'ne16pg3'. Please provide appropriate weight file for remapping."        
    fi    
fi

if [ "$do_point_diag" = true ]; then
    echo "Running point diagnostics for $case1"    
    python /nird/home/kjetisaa/Simple_Check_Spinup_Script/Check_Fates_Points_NIRD_refactored.py $case1 -o $outpath
fi

echo "Wrote files to: " $outpath
echo "Note that online path is accessable at https://ns9560k.web.sigma2.no/datalake/diagnostics/noresm/..."