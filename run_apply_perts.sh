#!/bin/sh

for mem in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
    #cp lndp2tile_input.C96.vgf.nml lndp2tile_input.nml
    #sed -i -e "s/XXMEM/$mem/g" lndp2tile_input.nml
    #./land-vector2tile/vector2tile_converter.exe lndp2tile_input.nml

    cp lndp_apply_input.C96.vgf.nml lndp_apply_input.nml
    sed -i -e "s/XXMEM/$mem/g"  lndp_apply_input.nml
    ./ensemble_pert/lndp_apply_pert.exe lndp_apply_input.nml

done
