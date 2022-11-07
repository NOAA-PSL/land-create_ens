#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=30
#SBATCH --job-name="stoch_unit_tests"

RES=96

cd stochastic_physics/unit_tests/

source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.stoch_gnu

# compile codes
sh compile_standalone.hera_gnu
if [ ! -f standalone_stochy.x ];then
  echo "compilation errors"
  exit 1
fi

# copy input directory
if [ ! -d INPUT ]; then
   cp -r /scratch2/BMC/gsienkf/Tseganeh.Gichamo/stochastic_physics_mod/unit_tests/INPUT  INPUT
fi

export OMP_NUM_THREADS=1
module list

for mem in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
    cp ../../stochy_input.C96.vgf.nml input.nml
    sed -i -e "s/LNDPSEED/$mem/g" input.nml
    time srun --label -n 6 standalone_stochy.x 
    mkdir ../../stochy_out0$mem
    mv workg* ../../stochy_out0${mem}
done
