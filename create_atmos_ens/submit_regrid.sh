#!/bin/bash
#SBATCH --ntasks=1  
#SBATCH --job-name=regrid_monthly
#SBATCH -t 00:30:00
#SBATCH -A gsienkf
#SBATCH -q debug
#SBATCH -J regrid_ensemble
#SBATCH -o log_noahmp.%j 
#SBATCH -e err_noahmp.%j

regrid_gefs_c96_monthly.x 
