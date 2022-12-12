# land-create_ens

Code to support ensemble offline land forecasts.

1. Perturbation of input parameters. 

This uses the stochastic physics unit test to create perturbation fields with assumed spatial correlations. Currently, we perturb the states once at the start of the forecast, by perturbing the values in the static_fields input file. 

Once we switch to using the model as a component of the UFS, we'll be able to use the stochastic physics component to apply the perturbations - which means that we can evolve them through time.

Note: if only perturbing once, and have an assumed spatial scale, this will create spatial variation in the magnitude of applied perturbations (i.e, some locations will have more spread, some less). If we were to retain the current approach of only applying perturbations once at the start, we should remove the spatial correlattions, so that every grid cell has the same perturbation spread (and if we're doing this, then we probably don't need to use the stochastic physics package). 

>run_generate_perts.sh creates files with the perturbations on the FV3 tile space  

input namelist example: stochy_input.C96.vgf.nml 

most of the relevent entries are under nam_sfcperts namelist

also in the fv_core_nml namelist, use layout=1,1 to get one file per tile. 

>run_apply_perts.sh converts the above perturbation patterns into vector format, and applies them to the static file.

input namelist example for comversion to tile: lndp2tile_input.C96.vgf.nml

lndp_layout needs to match fv_core_nml layout from above.

input namelist for application of the perturbations:  lndp_apply_input.C96.vgf.nml

2.  Generation of ensemble forcing from GEFS. This is done in create_atmos_ens. See README there.
