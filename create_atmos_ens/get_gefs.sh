#!/bin/bash 
#SBATCH --ntasks=1 -p service 
#SBATCH -t 12:00:00 
#SBATCH -A da-cpu
#SBATCH -q batch 
#SBATCH -J gefs_getfiles

script=$0
year_cur=$1
month_cur=$2
day_beg=$3
day_end=$4
iargs=$#
if [ ! $iargs = 4 ] ; then
    echo 'Usage: '$script' year month day_beg day_end'
    exit
else
    echo $script' '$1' '$2' '$3' '$4
fi

module load hpss

WGRIB2="/home/Fanglin.Yang/bin/wgrib2"
variables="(:TMP:|:PRES:|:RH:|:DSWRF:|:DLWRF:|:UGRD:|:VGRD:|:APCP:)"
#rundir=/scratch2/NCEPDEV/stmp3/Zhichang.Guo/GEFS
#grbdir=$rundir/grb2
rundir=/scratch2/NCEPDEV/land/Zhichang.Guo/GEFS/
grbdir=$rundir/grib2
ncdir=$rundir/nc4
extdir=$rundir/extracted
download="NO"
grib2nc="YES"
grib2grib="NO"
if [ ! -d $grbdir ] ; then mkdir $grbdir ; fi
if [ ! -d $ncdir ] ; then mkdir $ncdir ; fi
if [ ! -d $extdir ] ; then mkdir $extdir ; fi

yyyy=$year_cur
year_end=$year_cur
(( year_end += 1 ))
while [[ $yyyy -lt $year_end ]] ; do
   month_end=$month_cur
   (( month_end += 1 ))
   while [[ $month_cur -lt $month_end ]] ; do
       mm=$(printf "%02d" $month_cur)
       #day_beg=1
       #if [ $month_cur == 9 ] && [ $year_cur == 2020 ] ; then day_beg=23 ; fi
       #day_end=31
       #if [ $month_cur == 2 ] ; then day_end=28 ; fi
       #if [ $month_cur == 4 ] || [ $month_cur == 6 ] || [ $month_cur == 9 ] || [ $month_cur == 11 ] ; then day_end=30 ; fi
       day_cur=$day_beg
       (( day_end += 1 ))
       while [[ $day_cur -lt $day_end ]] ; do
           echo $yyyy $month_cur $day_cur
           dd=$(printf "%02d" $day_cur)
           year=$(($yyyy))
           mon=$(($month_cur))
           day=$(($day_cur+1))
           c1=0
           if [ $day == 32 ] && [ $mon == 1 ] ; then c1=1 ; fi
           if [ $day == 29 ] && [ $mon == 2 ] ; then c1=1 ; fi
           if [ $day == 32 ] && [ $mon == 3 ] ; then c1=1 ; fi
           if [ $day == 31 ] && [ $mon == 4 ] ; then c1=1 ; fi
           if [ $day == 32 ] && [ $mon == 5 ] ; then c1=1 ; fi
           if [ $day == 31 ] && [ $mon == 6 ] ; then c1=1 ; fi
           if [ $day == 32 ] && [ $mon == 7 ] ; then c1=1 ; fi
           if [ $day == 32 ] && [ $mon == 8 ] ; then c1=1 ; fi
           if [ $day == 31 ] && [ $mon == 9 ] ; then c1=1 ; fi
           if [ $day == 32 ] && [ $mon == 10 ] ; then c1=1 ; fi
           if [ $day == 31 ] && [ $mon == 11 ] ; then c1=1 ; fi
           if [ $day == 32 ] && [ $mon == 12 ] ; then c1=2 ; fi
           if [ $c1 == 1 ] ; then
             day=1
             mon=$(($mon+1))
           fi
           if [ $c1 == 2 ] ; then
             day=1
             mon=1
             year=$(($yyyy+1))
           fi
           day2=$(printf "%02d" $day)
           mon2=$(printf "%02d" $mon)
           ymd=$yyyy$mm$dd
           y4=$yyyy
           y4m2=$yyyy$mm
           ymd2=$year$mon2$day2
           tag=rh${yyyy}/$yyyy$mm/$ymd
           base=/NCEPPROD/5year/hpssprod/runhistory/$tag/
           dir=$grbdir/$ymd
           dir0=$extdir/$y4/$y4m2/$ymd
           dir1=$ncdir/$ymd
           dir2=$ncdir/$ymd2
           #echo $dir $dir0 $dir1 $dir2
           if [ "${grib2grib}" = "YES" ]; then
               if [ ! -d $dir0 ] ; then mkdir -p $dir0 ; fi
           fi
           if [ ! -d $dir ] ; then mkdir $dir ; fi
           if [ ! -d $dir1 ] ; then mkdir $dir1 ; fi
           if [ ! -d $dir2 ] ; then mkdir $dir2 ; fi
           for hh in {0..18..6} ; do
               hh=$(printf "%02d" $hh)
               archive=$base/com_gefs_prod_gefs.$yyyy$mm${dd}_${hh}.atmos_pgrb2sp25.tar 
               cd $dir
               if [ "${download}" = "YES" ]; then
                   htar -xvf $archive
               fi
               check='Y'
               for ee in {1..30..1} ; do
                   ee=$(printf "%02d" $ee)
                   check1=$dir/atmos/pgrb2sp25/gep${ee}.t${hh}z.pgrb2s.0p25.f003
                   check2=$dir/atmos/pgrb2sp25/gep${ee}.t${hh}z.pgrb2s.0p25.f006
                   if [ -f $check1 ] && [ -f $check2 ] ; then 
                       if [ "${download}" = "YES" ]; then
                           mv $check1 $dir
                           mv $check2 $dir
                       fi
                   fi
                   check3=$dir/gep${ee}.t${hh}z.pgrb2s.0p25.f003
                   check4=$dir/gep${ee}.t${hh}z.pgrb2s.0p25.f006
                   if [ ! -f $check3 ] || [ ! -f $check4 ] ; then
                       check='N'
                   fi
               done
               if [ $check == 'Y' ] ; then
                   old_dir=$dir/atmos
                   if [ -d $old_dir ] ; then
                       if [ "${download}" = "YES" ]; then
                           rm -rf $old_dir
                       fi
                   fi
               fi
           done
           for e in {1..30..1} ; do
               ee=$(printf "%02d" $e)
               for h in {0..18..6} ; do
                   hh=$(printf "%02d" $h)
                   for t in {3..6..3} ; do
                       tt=$(printf "%03d" $t)
                       time=$(($t+$h))
                       time2=$(printf "%02d" $time)
                       ifname=$dir/gep${ee}.t${hh}z.pgrb2s.0p25.f${tt}
                       if [ -f $ifname ] ; then 
                           efname=$dir0/egep${ee}.t${hh}z.pgrb2s.0p25.f${tt}.grib2
                           if [ $time -eq 24 ]; then
                               ofname=$dir2/gefs.ens${ee}.${year}-${mon2}-${day2}_00Z.nc
                           else
                               ofname=$dir1/gefs.ens${ee}.${yyyy}-${mm}-${dd}_${time2}Z.nc
                           fi
                           if [ "${grib2nc}" = "YES" ]; then
                               ${WGRIB2} $ifname | egrep $variables | ${WGRIB2} -i $ifname -netcdf ${ofname}
                           fi
                           if [ "${grib2grib}" = "YES" ]; then
                               ${WGRIB2} $ifname | egrep $variables | ${WGRIB2} -i $ifname -grib ${efname}
                           fi
                       fi
                   done
               done
           done
           (( day_cur += 1 ))
       done
       (( month_cur += 1 ))
   done
   (( yyyy += 1 ))
done
echo 'The script ended!'
