program get_gefs_monthly

  use netcdf
  implicit none

  integer, parameter :: lats      = 721    ! grid dimension for latitudes
  integer, parameter :: lons      = 1440   ! grid dimension for longitudes
  integer, parameter :: vds       = 8      ! variable size
  integer, parameter :: tds       = 8      ! 8 times (3-hourly) daily
  integer, parameter :: eds       = 30     ! ensemble size
  real,    parameter :: missing   = -999.  ! missing value for ensemble mean
!
  real               :: v3d(lons,lats,vds)
  real               :: v2d(lons,lats)
  real               :: q(lons,lats)
  real               :: avg(lons,lats,vds)
  real               :: cnt(lons,lats,vds)
  character(len=30)  :: ivnames(vds), ovnames(vds-1)
  character(len=256) :: dir, ifname, ifname_last, ofname, exefile
  character(len=4)   :: cyear
  character(len=10)  :: cdate, cdate_last
  character(len=8)   :: cy4m2d2, cy4m2d2_last
  character(len=2)   :: cmonth, cens2, chour, chour_last
  character(len=3)   :: cens3
  double precision   :: time_cur
  logical            :: file_exists, compress, found
  integer            :: monlen(12)
  integer            :: iargs, year, month, day, hour, iy4m2d2
  integer            :: idy, idy_beg, idy_end, eid, vid, tid, tstep
  integer            :: iyr_last, imo_last, idy_last
  integer            :: i, j, vid_u, vid_v, vid_p, vid_t, vid_rh, vid_sw, vid_lw, vid_prcp
  integer            :: incid, incid_last, oncid(eds+1)
  monlen  = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
!
! input variable names 
!
  ivnames = (/'TMP_2maboveground',  'PRES_surface',     'RH_2maboveground',    'DSWRF_surface', &
              'DLWRF_surface',      'APCP_surface',     'UGRD_10maboveground', 'VGRD_10maboveground'/)
!
! output variable names 
!
  ovnames = (/'temperature',        'surface_pressure', 'specific_humidity',   'solar_radiation', &
              'longwave_radiation', 'precipitation',    'wind_speed'/)
!
! base directory
!
  dir = '/scratch2/NCEPDEV/land/Zhichang.Guo/GEFS/'
!
! output internalally compressed netcdf file or not
!
  compress = .true.
!
! command line input (year and month)
!
  iargs = iargc()
  call getarg(0, exefile)
  call getarg(1, cyear)
  call getarg(2, cmonth)
  read(cyear,'(i4)') year
  read(cmonth,'(i2)') month
  idy_beg = 1
  if(year == 2020 .and. month == 9) then
      idy_beg = 23
      hour = 3
  endif
  iy4m2d2 = year*10000 + month*100 + idy_beg
  if(iy4m2d2 < 20200923 .or. year > 2022 .or. month < 1 .or. month > 12) then
     print*, 'something is wrong with date'
     stop
  endif
!
! Modify this for the leap year 
!
  idy_end = monlen(month) 
  write(cyear,'(i4)') year
  write(cmonth,'(i2.2)') month
!
! find the id for P/T/RH/U/V
!
  vid_p    = -1
  vid_t    = -1
  vid_rh   = -1
  vid_u    = -1
  vid_v    = -1
  vid_sw   = -1
  vid_lw   = -1
  vid_prcp = -1
  do vid = 1, vds
    if(trim(ivnames(vid)) == 'PRES_surface')        vid_p    = vid
    if(trim(ivnames(vid)) == 'TMP_2maboveground')   vid_t    = vid
    if(trim(ivnames(vid)) == 'RH_2maboveground')    vid_rh   = vid
    if(trim(ivnames(vid)) == 'UGRD_10maboveground') vid_u    = vid
    if(trim(ivnames(vid)) == 'VGRD_10maboveground') vid_v    = vid
    if(trim(ivnames(vid)) == 'DSWRF_surface')       vid_sw   = vid
    if(trim(ivnames(vid)) == 'DLWRF_surface')       vid_lw   = vid
    if(trim(ivnames(vid)) == 'APCP_surface')        vid_prcp = vid
  enddo
  if(vid_p < 0 .or. vid_t < 0  .or. vid_rh < 0 .or. vid_u < 0 .or. &
     vid_v < 0 .or. vid_sw < 0 .or. vid_lw < 0 .or. vid_prcp < 0) then
     print*, 'Variables P/T/RH/U/V/SW/LW/PRCP ids are not found'
     stop 10
  endif
!
! defind the output netcdf files
!
  do eid = 1, eds+1
    if(eid == eds+1) then
      if(compress) then
        ofname = trim(dir)//'/monthly'//'/gefs_forcing.'//trim(cyear)//'-'//trim(cmonth)//'.nc4'
      else
        ofname = trim(dir)//'/monthly'//'/gefs_forcing.'//trim(cyear)//'-'//trim(cmonth)//'.nc'
      endif
    else
      write(cens3,'(i3.3)') eid
      if(compress) then
        ofname = trim(dir)//'/monthly'//'/gefs.ens'//trim(cens3)//'.'//trim(cyear)//'-'//trim(cmonth)//'.nc4'
      else
        ofname = trim(dir)//'/monthly'//'/gefs.ens'//trim(cens3)//'.'//trim(cyear)//'-'//trim(cmonth)//'.nc'
      endif
    endif
    call define_output_file_one(ofname, lons, lats, vds-1, ovnames, oncid(eid), compress)
  enddo
!
! for each day of the month, read in variables, and do the following calculations:
!     1. calculate wind speed from U/V
!     2. calculate specific humidity from T/P/RH
!     3. calculate ensemble mean for each variables
!     4. calculate ensemble mean specific from ensemble mean T/P/RH
!
  tstep = 0
  do idy = idy_beg, idy_end
    write(cy4m2d2,'(i4,2i2.2)') year, month, idy
    write(cdate,'(i4,a1,i2.2,a1,i2.2)') year, '-', month, '-', idy
    iyr_last = year
    imo_last = month
    idy_last = idy - 1
    if(idy_last < 1) then
      imo_last = imo_last - 1
      if(imo_last < 1) then
        iyr_last = iyr_last - 1
        imo_last = 12
      endif
      idy_last = monlen(imo_last)
    endif
    print*, year, month, idy, iyr_last, imo_last, idy_last
    write(cy4m2d2_last,'(i4,2i2.2)') iyr_last, imo_last, idy_last
    write(cdate_last,'(i4,a1,i2.2,a1,i2.2)') iyr_last, '-', imo_last, '-', idy_last
    do tid = 1, tds
      write(chour,'(i2.2)') (tid-1)*3
      avg(:,:,:) = 0.0
      cnt(:,:,:) = 0.0
      do eid = 1, eds
        write(cens2,'(i2.2)') eid
        write(cens3,'(i3.3)') eid
        ifname  = trim(dir)//'/nc4/'//trim(cy4m2d2)//'/gefs.ens'//trim(cens2)//'.'//trim(cdate)//'_'//trim(chour)//'Z.nc'
        inquire(file=ifname, exist=file_exists)
        if(file_exists) then
          call open_source_file(ifname, incid)
          if(eid == 1) then
            tstep = tstep + 1
            call read_time(incid, time_cur)
          endif
          call output_time(oncid(eid), tstep, time_cur)
!
!         read in variables
!
          call read_source_variable(incid, lons, lats, vds, v3d(:,:,:), ivnames)
!
!         calculate wind speed and relative humidity, then output variables
!
          do vid = 1, vds-1
            if(vid == vid_rh) then
              call cal_rh2sh(lons, lats, v3d(:,:,vid_t), v3d(:,:,vid_p), v3d(:,:,vid_rh), q)
              call output_variable(oncid(eid), lons, lats, tstep, q(:,:), trim(ovnames(vid)))
            else
              if(vid == vid_u) then
                do i = 1, lons
                  do j = 1, lats
                    v3d(i,j,vid) = sqrt(v3d(i,j,vid_u)*v3d(i,j,vid_u) + v3d(i,j,vid_v)*v3d(i,j,vid_v))
                  enddo
                enddo
              else if(vid == vid_prcp .or. vid == vid_lw .or. vid == vid_sw) then
! tid = 1/3/5/7 are 6-hourly forecasts (00Z, 06Z, 12Z, 18Z)
! tid = 2/4/6/8 are 3-hourly forecasts (03Z, 09Z, 15Z, 21Z)
! deal with 6-hourly average of short-wave, long-wave radiation and precipitation
! second 3-hourly average = 2 * 6-hourly average - first 3-hourly average
                if(mod(tid-1,2) == 0) then
                  if(tid == 1) then
                    write(chour_last,'(i2.2)') (tid-2)*3 + 24
                    ifname_last = trim(dir)//'/nc4/'//trim(cy4m2d2_last)//'/gefs.ens'//trim(cens2)// & 
                                  '.'//trim(cdate_last)//'_'//trim(chour_last)//'Z.nc'
                  else
                    write(chour_last,'(i2.2)') (tid-2)*3
                    ifname_last = trim(dir)//'/nc4/'//trim(cy4m2d2)//'/gefs.ens'//trim(cens2)//'.'// &
                                  trim(cdate)//'_'//trim(chour_last)//'Z.nc'
                  endif
                  inquire(file=ifname_last, exist=found)
                  if(.not. found) then
                    print*, "last file "//trim(ifname_last)//' not found', cy4m2d2, year, month, idy
                  endif
                  call open_source_file(ifname_last, incid_last)
                  call read_source_variable_one(incid_last, lons, lats, v2d(:,:), ivnames(vid))
                  call close_file(incid_last)
                  do i = 1, lons
                    do j = 1, lats
                      v3d(i,j,vid) = 2.0*v3d(i,j,vid) - v2d(i,j)
                    enddo
                  enddo
                endif
              endif
              call output_variable(oncid(eid), lons, lats, tstep, v3d(:,:,vid), trim(ovnames(vid)))
            endif
!
!           accumulate each variable for calculating ensemble mean
!
            do i = 1, lons
              do j = 1, lats
                avg(i,j,vid) = avg(i,j,vid) + v3d(i,j,vid)
                cnt(i,j,vid) = cnt(i,j,vid) + 1.0
              enddo
            enddo
          enddo
          call close_file(incid)
        else
          write(*,'(A110)') 'Warning: '//trim(ifname)//' not found'
        endif
      enddo
!
!     calculate ensemble mean for each variable
!
      do vid = 1, vds-1
        do i = 1, lons
          do j = 1, lats
            if(cnt(i,j,vid) > 0.5) then
              avg(i,j,vid) = avg(i,j,vid)/cnt(i,j,vid)
            else
              avg(i,j,vid) = missing
            endif
          enddo
        enddo
      enddo
!
!     calculate specific humidity ensemble mean from ensemble mean of T/P/RH
!
      call cal_rh2sh(lons, lats, avg(:,:,vid_t), avg(:,:,vid_p), avg(:,:,vid_rh), q)
!
!     output ensemble mean
!
      do vid = 1, vds-1
        if(vid == vid_rh) then
          call output_variable(oncid(eds+1), lons, lats, tstep, q(:,:),       trim(ovnames(vid)))
        else
          call output_variable(oncid(eds+1), lons, lats, tstep, avg(:,:,vid), trim(ovnames(vid)))
        endif
      enddo
    enddo   ! each time step
  enddo     ! each day
  do eid = 1, eds+1
    call close_file(oncid(eid))
  enddo
  print*, 'Program ended normally!'
end program get_gefs_monthly
!--------------------------------------------
subroutine define_output_file_one(filename, lons, lats, vds, vnames, ncid, compress)
  use netcdf
  implicit none
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: lons, lats, vds
  character(len=*), intent(in)  :: vnames(vds)
  integer,          intent(out) :: ncid
  logical,          intent(in)  :: compress
  character(len=80)             :: lname
  character(len=3)              :: cens3
  integer                       :: dim_id_i, dim_id_j, dim_id_t  ! netcdf dimension identifiers
  integer                       :: status, vid, varid
  if(compress) then
    status = nf90_create(filename, NF90_NETCDF4, ncid)
  else
    status = nf90_create(filename, NF90_CLOBBER, ncid)
  endif
    if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

  status = nf90_def_dim(ncid, "lon"   , lons , dim_id_i)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "lat"   , lats , dim_id_j)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "time"   , nf90_unlimited , dim_id_t)
    if (status /= nf90_noerr) call handle_err(status)

! Define variables in the file.

  status = nf90_def_var(ncid, 'time', NF90_DOUBLE, (/dim_id_t/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, 'units', 'seconds since 1970-01-01 00:00:00')
    if (status /= nf90_noerr) call handle_err(status)
  do vid = 1, vds
    call get_long_name(vnames(vid), lname)
    status = nf90_def_var(ncid, trim(vnames(vid)), NF90_FLOAT, (/dim_id_i, dim_id_j, dim_id_t/), varid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, 'long_name', trim(lname))
      if (status /= nf90_noerr) call handle_err(status)

    if(compress) then
      status = nf90_def_var_deflate(ncid, varid, 1, 1, 5)
        if (status /= nf90_noerr) call handle_err(status)
    endif
  enddo
  status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)
end subroutine define_output_file_one
!--------------------------------------------
subroutine open_source_file(source_filename, ncid)
  use netcdf
  implicit none
  character(len=*), intent(in)  :: source_filename
  integer,          intent(out) :: ncid
  integer                       :: status
  status = nf90_open(source_filename, NF90_NOwrite, ncid)
    if (status /= nf90_noerr) call handle_err(status)
end subroutine open_source_file
!--------------------------------------------
subroutine output_time(ncid, tid, tstep)
  use netcdf
  implicit none
  integer,          intent(in) :: ncid, tid
  double precision, intent(in) :: tstep
  integer                      :: status, varid
  status = nf90_inq_varid(ncid, 'time', varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid, tstep, start = (/tid/))
    if (status /= nf90_noerr) call handle_err(status)
end subroutine output_time
!--------------------------------------------
subroutine read_time(ncid, tstep)
  use netcdf
  implicit none
  integer,          intent(in)  :: ncid
  double precision, intent(out) :: tstep
  integer                       :: status, varid
  status = nf90_inq_varid(ncid, 'time', varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid, tstep, start = (/1/))
    if (status /= nf90_noerr) call handle_err(status)
end subroutine read_time
!--------------------------------------------
subroutine read_source_variable(ncid, lons, lats, vds, source_input, vnames)
  use netcdf
  implicit none
  integer,          intent(in)  :: ncid, lons, lats, vds
  character(len=*), intent(in)  :: vnames(vds)
  real,             intent(out) :: source_input(lons,lats,vds)
  integer                       :: status, varid, vid
  do vid = 1, vds
    status = nf90_inq_varid(ncid, vnames(vid), varid)
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var(ncid, varid, source_input(:,:,vid), start = (/1,1,1/), count = (/lons,lats,1/))
      if (status /= nf90_noerr) call handle_err(status)
  enddo
end subroutine read_source_variable
!--------------------------------------------
subroutine read_source_variable_one(ncid, lons, lats, source_input, vname)
  use netcdf
  implicit none
  integer,          intent(in)  :: ncid, lons, lats
  character(len=*), intent(in)  :: vname
  real,             intent(out) :: source_input(lons,lats)
  integer                       :: status, varid, vid
  status = nf90_inq_varid(ncid, vname, varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid, source_input(:,:), start = (/1,1,1/), count = (/lons,lats,1/))
    if (status /= nf90_noerr) call handle_err(status)
end subroutine read_source_variable_one
!
!----------- calculate specific humidity from relative humidity -----------
!
subroutine cal_rh2sh(nx, ny, t, p, rh, q)
  implicit none
  real,    parameter   :: rd  = 287.058 ! dry air constant J/(K*kg)
  real,    parameter   :: rv  = 461.5   ! J/(K*kg)
  real,    parameter   :: t0  = 273.15  ! K
  real,    parameter   :: es0 = 611.0   !
  integer, intent(in)  :: nx, ny        ! dimensions
  real,    intent(in)  :: t(nx,ny)      ! temperature in K (or 273.15 + degree)
  real,    intent(in)  :: p(nx,ny)      ! pressure (Pa)
  real,    intent(in)  :: rh(nx,ny)     ! input: relative humidity (%)
  real,    intent(out) :: q(nx,ny)      ! output: specific humidity (kg/kg)
  real                 :: es            ! saturation pressure (pa) over ice and water 
                                        ! respectively, Bohren & Albrecht
                                        ! 2000, pp 197-200
  real                 :: ws, w
  integer              :: i, j
  do i = 1, nx
    do j = 1, ny
      if(t(i,j) < t0) then
        es = es0*exp(6293./t0 - 6293./t(i,j) - 0.555*log(t(i,j)/t0))
      else
        es = es0*exp(6808./t0 - 6808./t(i,j) - 5.09*log(t(i,j)/t0))
      endif
      ws = rd/rv*es/(p(i,j)-es)
      w = ws * rh(i,j) * 0.01
      q(i,j) = w/(1.+w)
    enddo
  enddo
end subroutine cal_rh2sh
!--------------------------------------------
subroutine output_variable(ncid, lons, lats, tid, var, vname)
  use netcdf
  implicit none
  integer,          intent(in) :: ncid, lons, lats, tid
  character(len=*), intent(in) :: vname
  real,             intent(in) :: var(lons,lats)
  integer                      :: status, varid
  status = nf90_inq_varid(ncid, vname, varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , var, start = (/1,1,tid/), count = (/lons,lats,1/))
    if (status /= nf90_noerr) call handle_err(status)
end subroutine output_variable
!--------------------------------------------
subroutine close_file(ncid)
  use netcdf
  implicit none
  integer, intent(in) :: ncid
  integer             :: status
  status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
end subroutine close_file
!--------------------------------------------
subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine handle_err
!--------------------------------------------
subroutine get_long_name(vname, key)
  implicit none
  character(len=*), intent(in)  :: vname
  character(len=*), intent(out) :: key
  if(trim(vname) == "temperature")        key = "near-surface air temperature [K]"
  if(trim(vname) == "surface_pressure")   key = "surface pressure [Pa]"
  if(trim(vname) == "specific_humidity")  key = "surface specific humidity [kg/kg]"
  if(trim(vname) == "solar_radiation")    key = "surface downward solar radiation [W/m^2]"
  if(trim(vname) == "longwave_radiation") key = "surface downward long-wave radiation [W/m^2]"
  if(trim(vname) == "wind_speed")         key = "surface wind speed [m/s]"
  if(trim(vname) == "precipitation")      key = "precipitation [mm/s]"
end subroutine get_long_name
