Code to prepare GEFS ensemble forcing files. 
1. Reads in GEFS ensemble forcing from individual files for each ensemble member and each day, and collates into a single ensemble monthly file. 
2. Regrids the forcing.

Clara Draper, Nov 2022. Original code by Zhichang Guo.

older_options contains older code to collate monthly files, and collate monthly files with all ensemble members (no re-gridding). 

To-do: generalize current code so don't need above versions. Generalize grid size. 

To compile:

> compile.sh

To run: 
set dates, input and output directory in regrid_ens.nml. 
> sbatch submit_regrid.sh

December 12 - not working, as needs updated weights file.
