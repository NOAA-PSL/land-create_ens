#! /usr/bin/env bash
set -eux

# check if part of workflow. If so, use those modules.
source hera_modules

export FCMP=mpiifort

export INCS="-I${NETCDF}/include"          # netcdf modules does not
#export FFLAGS="$INCS -O3 -fp-model precise -r8 -convert big_endian -traceback -g"
export FFLAGS="$INCS"

export LIBSM="-L${NETCDF}/lib -lnetcdff"

file=regrid_gefs_c96_monthly
$FCMP $FFLAGS -c ${file}.f90 -o $file.o
$FCMP -o $file.x $file.o $LIBSM
