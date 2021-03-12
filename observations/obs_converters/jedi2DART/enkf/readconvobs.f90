module readconvobs
!$$$  module documentation block
!
! module: readconvobs...follows gsi/enkf/readconvobs.f90 and adpated for JEDI observation output
!
! abstract: read data from diag_conv* files (containing prepbufr data) written out
!  by GSI forward operator code.
!
! Public Subroutines:
!  get_num_convobs: determine the number of observations to read.
!  get_convobs_data: read the data.
!   
! Public Variables: None
!
!$$$
use kinds, only: r_kind,i_kind,r_single, r_double
use constants, only: one,zero, fv
use params, only: num_ufo_procs, obs_platforms, num_obs_platforms, process_sondes_Tv
use netcdf_mod, only: open_netcdf, close_netcdf, get_netcdf_dims, get_netcdf_var_1d
use mpisetup, only: stop2
use jedi_time_routines_mod, only : read_jedi_time
use genqsat_mod, only : genqsat1

implicit none

private

! public subroutines
public :: get_num_convobs, get_convobs_data

! observation platforms/types to read from netcdf files are a namelist variable read in
!character(len=50), dimension(num_obs_platforms):: obs_platforms = (/'aircraft', 'gnssro', 'satwind', 'sfc', 'sondes' /)

! when processing eastward_wind, also process northward_wind
integer(i_kind), parameter :: num_met_variables = 7
character(len=100), dimension(num_met_variables):: met_variables = (/'air_temperature', 'eastward_wind', &
                                                                  'specific_humidity', 'refractivity', 'surface_pressure', &
                                                                  'bending_angle', 'virtual_temperature' /)
!integer(i_kind), parameter :: num_met_variables = 7
!character(len=100), dimension(num_met_variables):: met_variables = (/'air_temperature', 'eastward_wind', 'northward_wind', &
!                                                                  'specific_humidity', 'refractivity', 'surface_pressure', &
!                                                                  'bending_angle' /)
contains

subroutine get_num_convobs(obspath,datestring,num_obs_tot,id)
    character (len=500), intent(in) :: obspath
    character (len=10), intent(in) :: datestring
    character(len=10), intent(in) :: id
    integer(i_kind), intent(out) :: num_obs_tot

    character(len=500) obsfile
    character(len=4) :: pe_name
    character(len=100) :: met_variable, this_variable
    character(len=50) :: ob_platform
    integer(i_kind) :: i, itype, ipe, j
    integer(i_kind) :: nobs_curr, num_obs_totdiag
    integer(i_kind),dimension(num_met_variables,2) :: nobs
    integer(i_kind) :: ncfileid
    real(r_kind) :: errorlimit,errorlimit2,error,pres,obmax
    real(r_kind) :: errorlimit2_obs,errorlimit2_bnd
    logical :: fexist, skip

    integer(i_kind), allocatable, dimension (:) :: effectiveQC, effectiveQC_V
    real(r_single), allocatable, dimension (:) :: pressure, altitude
    real(r_single), allocatable, dimension (:) :: finalObsError, GPS_Type
    real(r_single), allocatable, dimension (:) :: observation, observationV

    ! If ob error > errorlimit or < errorlimit2, skip it.
    errorlimit = 1._r_kind/sqrt(1.e9_r_kind)
    errorlimit2_obs = 1._r_kind/sqrt(1.e-6_r_kind)
    errorlimit2_bnd = 1.e3_r_kind*errorlimit2_obs
    num_obs_tot = 0
    num_obs_totdiag = 0
    nobs = 0

    obtypeloop: do itype=1, num_obs_platforms

       ob_platform = obs_platforms(itype)

       peloop: do ipe=0,num_ufo_procs-1 ! PEs start at 0, hence the -1

          write(pe_name,'(i4.4)') ipe

          ! Need to incorporate $id into filename
          if ( num_ufo_procs == 1 ) then
            !obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_all_'//trim(adjustl(id))//'.nc4'
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//trim(adjustl(id))//'_all.nc4'
          else
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
             inquire(file=obsfile,exist=fexist)
             if (.not.fexist) obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//trim(adjustl(id))//'_'//pe_name//'.nc4'
          endif

          ! Make sure file exists
          inquire(file=obsfile,exist=fexist)
          if (.not.fexist) cycle peloop

          ! File exists.  Open it, figure out how many obs this file has.
          call open_netcdf(obsfile,ncfileid)
          call get_netcdf_dims(ncfileid,'nlocs',nobs_curr)
          if ( nobs_curr .le. 0 ) then
             call close_netcdf(obsfile,ncfileid)
             cycle peloop ! If the file has no obs go to next file
          endif
          !write(*,*)' nlocs = ',nobs_curr

          allocate(pressure(nobs_curr), effectiveQC(nobs_curr))
          allocate(finalObsError(nobs_curr), observation(nobs_curr))
          allocate(observationV(nobs_curr),effectiveQC_V(nobs_curr))
          allocate(altitude(nobs_curr))

          if ( trim(adjustl(ob_platform)) .eq. 'gnssro') then
            !call get_netcdf_var_1d(ncfileid,'altitude@MetaData',nobs_curr,altitude)
             pressure = 500.
          else
             call get_netcdf_var_1d(ncfileid,'air_pressure@MetaData',nobs_curr,pressure) ! pressure in Pa
             pressure = pressure * 0.01 ! get into hPa
          endif
  
          do j = 1,num_met_variables

             met_variable = met_variables(j)

             call check_met_variables_for_this_platform(ob_platform,met_variable,skip) ! output is boolean $skip
             if ( skip ) cycle

             this_variable = trim(adjustl(met_variable))//'@ObsValue'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,observation)
             
             ! if we're processing U wind, also get V wind
             if ( trim(adjustl(met_variable)).eq.'eastward_wind' ) then
                call get_netcdf_var_1d(ncfileid,'northward_wind@ObsValue',nobs_curr,observationV)
                call get_netcdf_var_1d(ncfileid,'northward_wind@EffectiveQC',nobs_curr,effectiveQC_V)
             endif
             
             this_variable = trim(adjustl(met_variable))//'@EffectiveQC' ! this may be all that's needed
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,effectiveQC)

             this_variable = trim(adjustl(met_variable))//'@EffectiveError'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,finalObsError)
             
             ! convert Tv to T
             if ( trim(adjustl(met_variable)).eq.'virtual_temperature' ) then
               call convert_tv_to_temp('obs',ncfileid,nobs_curr,observation,effectiveQC) ! overwrite observation,effectiveQC
             endif

             do i = 1,nobs_curr
                nobs(j,1) = nobs(j,1) + 1  ! number of read observations
                num_obs_totdiag = num_obs_totdiag + 1 ! total number of obs read across all types/platforms...not necessarily total number used
                obmax = abs(observation(i))
                if ( trim(adjustl(met_variable)).eq.'eastward_wind' ) then
                   nobs(j,1) = nobs(j,1) + 1 ! account for v-wind also
                   num_obs_totdiag = num_obs_totdiag + 1 
                   obmax = max(abs(observation(i)), abs(observationV(i)))
                   effectiveQC(i) = max( effectiveQC(i), effectiveQC_V(i))
                endif

                if (effectiveQC(i) > 0 ) cycle
                if ( trim(adjustl(met_variable)).eq.'surface_pressure' ) then
                   pres = observation(i) * 0.01 ! Pa --> hPa
                else
                   pres = pressure(i)
                endif
                if ( pres < 0.001_r_kind .or. pres > 1200._r_kind .or. &
                         abs(obmax) > 1.e9_r_kind ) cycle
                ! If we're this far, the ob is good and can be assimilated.
                nobs(j,2) = nobs(j,2) + 1  ! number of "kept" observations for this met_variable
                num_obs_tot = num_obs_tot + 1 ! total number of "kept" observations across all types/platforms
                if ( trim(adjustl(met_variable)).eq.'eastward_wind' ) then
                   nobs(j,2) = nobs(j,2) + 1  ! account for v-wind also
                   num_obs_tot = num_obs_tot + 1 
                endif
             enddo

          enddo ! loop over "j" (number of obs for this specific platform and meteorological type)

          ! clean-up
          call close_netcdf(obsfile,ncfileid)
          deallocate(pressure, effectiveQC, finalObsError, observation, observationV, altitude, effectiveQC_V)

       enddo peloop ! end loop over "ipe"
    enddo obtypeloop ! end loop over "itype"

    write(6,*),num_obs_totdiag, ' conventional obs read from files'
    write(6,*),num_obs_tot,' conventional obs actually kept from files'
    write(6,*)'columns below obtype,nread, nkeep'
    do j = 1, num_met_variables
       write(6,100) met_variables(j), nobs(j,1), nobs(j,2)
    enddo
100       format(2x,a20,2x,i9,2x,i9)

end subroutine get_num_convobs

subroutine get_convobs_data(obspath, datestring, nobs_max, h_x_ensmean, h_xnobc, x_obs, x_err, &
           x_lon, x_lat, x_press, x_time, x_code, x_errorig, x_type, id, id2)

    character*500, intent(in) :: obspath
    character*10, intent(in) :: datestring
    integer(i_kind), intent(in) :: nobs_max
    character(len=10), intent(in) :: id,id2

    real(r_single), dimension(nobs_max) :: h_x_ensmean,h_xnobc,x_obs,x_err,x_lon,&
                               x_lat,x_press,x_time,x_errorig
    integer(i_kind), dimension(nobs_max) :: x_code
    character(len=20), dimension(nobs_max) ::  x_type

    integer(i_kind) :: i, itype, ipe, j, nob
    real(r_kind) :: errorlimit,errorlimit2,error,pres
    real(r_kind) :: errorlimit2_obs,errorlimit2_bnd, obmax
    real(r_kind) :: tolerance = 1.e-5
    real(r_kind) :: ges_tv(1,1), ges_prsl(1,1)
    real(r_double) :: qsat(1,1)
    logical :: twofiles, fexist, skip

    character(len=500) obsfile, obsfile2
    character(len=4) :: pe_name
    character(len=100) :: met_variable, this_variable
    character(len=50) :: ob_platform
    integer(i_kind) :: nobs_curr, nobs_curr2, num_obs_totdiag
    integer(i_kind) :: ncfileid, ncfileid2

    integer(i_kind), allocatable, dimension (:) :: effectiveQC, effectiveQC_V
    integer(i_kind), allocatable, dimension (:) :: Observation_Type, Observation_Type2 ! prepbufr code
    real(r_single), allocatable, dimension (:) :: altitude, altitude2
    real(r_single), allocatable, dimension (:) :: initialObsError, finalObsError
    real(r_single), allocatable, dimension (:) :: observation, observationV
    real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_unadjustedV, Obs_Minus_Forecast_unadjustedV2
!   real(r_single), allocatable, dimension (:) :: Forecast_Saturation_Spec_Hum
    real(r_single), allocatable, dimension (:) :: aux_Data, aux_Data2
    real(r_single), allocatable, dimension (:) :: latitude, longitude, pressure, time
    real(r_single), allocatable, dimension (:) :: latitude2, longitude2, pressure2, time2
    real(r_single), allocatable, dimension (:) :: GPS_Type
    real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_adjusted
    real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_unadjusted
    real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_adjusted2
    real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_unadjusted2
    character(len=20), allocatable, dimension(:) :: jedi_time_string, jedi_time_string2

    ! Error limit is made consistent with screenobs routine
    errorlimit = 1._r_kind/sqrt(1.e9_r_kind)
    errorlimit2_obs = 1._r_kind/sqrt(1.e-6_r_kind)
    errorlimit2_bnd = 1.e3_r_kind*errorlimit2_obs

    twofiles = id2 /= id
    
    nob = 0
    
    obtypeloop: do itype=1, num_obs_platforms

       ob_platform = obs_platforms(itype)

       peloop: do ipe=0,num_ufo_procs-1 ! PEs start at 0, hence the -1

          write(pe_name,'(i4.4)') ipe

          ! Need to incorporate $id into filename
          if ( num_ufo_procs == 1 ) then
            !obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_all_'//trim(adjustl(id))//'.nc4'
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//trim(adjustl(id))//'_all.nc4'
          else
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
             inquire(file=obsfile,exist=fexist)
             if (.not.fexist) obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//trim(adjustl(id))//'_'//pe_name//'.nc4'
          endif

          ! Make sure file exists
          inquire(file=obsfile,exist=fexist)
          if (.not.fexist) cycle peloop

          ! File exists.  Open it, figure out how many obs this file has.
          call open_netcdf(obsfile,ncfileid)
          call get_netcdf_dims(ncfileid,'nlocs',nobs_curr)
          if ( nobs_curr .le. 0 ) then
             call close_netcdf(obsfile,ncfileid)
             cycle peloop ! If the file has no obs go to next file
          endif
          
          if (twofiles) then
             ! Need to incorporate $id2 into filename
             if ( num_ufo_procs == 1 ) then
               !obsfile2 = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_all_'//trim(adjustl(id2))//'.nc4'
                obsfile2 = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//trim(adjustl(id2))//'_all.nc4'
             else
                obsfile2 = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//pe_name//'_'//trim(adjustl(id2))//'.nc4'
                inquire(file=obsfile2,exist=fexist)
                if (.not.fexist) obsfile2 = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//trim(adjustl(id2))//'_'//pe_name//'.nc4'
             endif

             ! Make sure file exists
             inquire(file=obsfile2,exist=fexist)
             if (.not.fexist) then
                write (6,*) ' ens mean file here but ...'
                write (6,*) ' missing member file ',trim(adjustl(obsfile2))
                call close_netcdf(obsfile,ncfileid) ! close ensemble mean file
                call stop2(73)
             endif

             ! Second file exists.  Open it, figure out how many obs this file has.
             call open_netcdf(obsfile2,ncfileid2)
             call get_netcdf_dims(ncfileid2,'nlocs',nobs_curr2)
             if ( nobs_curr2 .le. 0 .or. (nobs_curr .ne. nobs_curr2) ) then
                write (6,*) ' mismatching obs numbers'
                write (6,*) ' nobs, nobs2 = ',nobs_curr,nobs_curr2
                call close_netcdf(obsfile,ncfileid)
                call close_netcdf(obsfile2,ncfileid2)
                call stop2(74)
             endif
          endif

          allocate(latitude(nobs_curr), longitude(nobs_curr), pressure(nobs_curr), time(nobs_curr), &
                 effectiveQC(nobs_curr), initialObsError(nobs_curr), finalObsError(nobs_curr), &
                 Observation_Type(nobs_curr), observation(nobs_curr), altitude(nobs_curr), &
                 Obs_Minus_Forecast_adjusted(nobs_curr), Obs_Minus_Forecast_unadjusted(nobs_curr), &
                 jedi_time_string(nobs_curr),effectiveQC_V(nobs_curr))
          allocate(aux_Data(nobs_curr),aux_Data2(nobs_curr))
          allocate(Observation_Type2(nobs_curr),Obs_Minus_Forecast_unadjusted2(nobs_curr))
          allocate(Obs_Minus_Forecast_adjusted2(nobs_curr))
          allocate(observationV(nobs_curr),Obs_Minus_Forecast_unadjustedV(nobs_curr),Obs_Minus_Forecast_unadjustedV2(nobs_curr))

          call get_netcdf_var_1d(ncfileid,'latitude@MetaData',nobs_curr,latitude)
          call get_netcdf_var_1d(ncfileid,'longitude@MetaData',nobs_curr,longitude)
          call get_netcdf_var_1d(ncfileid,'datetime@MetaData',nobs_curr,jedi_time_string)
          ! get difference between analysis time and observation times (hrs)
          call read_jedi_time(datestring,nobs_curr,jedi_time_string,time) ! output is time
          if ( trim(adjustl(ob_platform)) .eq. 'gnssro') then
             call get_netcdf_var_1d(ncfileid,'altitude@MetaData',nobs_curr,altitude)
             pressure = 500. ! not used if GNSSRO--use altitude for vertical level
          else
             call get_netcdf_var_1d(ncfileid,'air_pressure@MetaData',nobs_curr,pressure) ! pressure in Pa
             pressure = pressure * 0.01 ! get into hPa
          endif

          if (twofiles) then
             allocate(latitude2(nobs_curr), longitude2(nobs_curr), pressure2(nobs_curr), time2(nobs_curr))
             allocate(jedi_time_string2(nobs_curr),altitude2(nobs_curr))
             call get_netcdf_var_1d(ncfileid2,'latitude@MetaData',nobs_curr,latitude2)
             call get_netcdf_var_1d(ncfileid2,'longitude@MetaData',nobs_curr,longitude2)
             call get_netcdf_var_1d(ncfileid2,'datetime@MetaData',nobs_curr,jedi_time_string2)
             ! get difference between analysis time and observation times (hrs)
             call read_jedi_time(datestring,nobs_curr,jedi_time_string2,time2) ! output is time2
             if ( trim(adjustl(ob_platform)) .eq. 'gnssro') then
                call get_netcdf_var_1d(ncfileid2,'altitude@MetaData',nobs_curr,altitude2)
                pressure2 = 500. ! not used if GNSSRO--use altitude for vertical level
             else
                call get_netcdf_var_1d(ncfileid2,'air_pressure@MetaData',nobs_curr,pressure2) ! pressure in Pa
                pressure2 = pressure2 * 0.01 ! get into hPa
             endif
          endif

          do j = 1,num_met_variables

             met_variable = met_variables(j)

             call check_met_variables_for_this_platform(ob_platform,met_variable,skip) ! output is boolean $skip
             if ( skip ) cycle

             this_variable = trim(adjustl(met_variable))//'@ObsValue'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,observation)

             this_variable = trim(adjustl(met_variable))//'@EffectiveQC'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,effectiveQC)

             this_variable = trim(adjustl(met_variable))//'@EffectiveError'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,finalObsError)

             ! convert Tv to T
             if ( trim(adjustl(met_variable)).eq.'virtual_temperature' ) then
               call convert_tv_to_temp('obs',ncfileid,nobs_curr,observation,effectiveQC) ! overwrite observation,effectiveQC
             endif

             this_variable = trim(adjustl(met_variable))//'@ObsError'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,initialObsError)

             !this_variable = trim(adjustl(met_variable))//'@depbg' ! make sure this is y-h(x) and not h(x) -y 
             !call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,Obs_Minus_Forecast_adjusted) ! y-h(x)

             this_variable = trim(adjustl(met_variable))//'@hofx'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,Obs_Minus_Forecast_unadjusted) ! h(x)
             
             if ( trim(adjustl(met_variable)).eq.'virtual_temperature' ) then
               call convert_tv_to_temp('hx',ncfileid,nobs_curr,Obs_Minus_Forecast_unadjusted) ! overwrite Obs_Minus_Forecast_unadjusted
             endif
             Obs_Minus_Forecast_unadjusted = observation - Obs_Minus_Forecast_unadjusted ! y - h(x)

             if ( trim(adjustl(met_variable)).eq.'eastward_wind' ) then
                call get_netcdf_var_1d(ncfileid,'northward_wind@ObsValue',nobs_curr,observationV)
                call get_netcdf_var_1d(ncfileid,'northward_wind@hofx',nobs_curr,Obs_Minus_Forecast_unadjustedV)
                Obs_Minus_Forecast_unadjustedV = observationV - Obs_Minus_Forecast_unadjustedV ! y - h(x)
                call get_netcdf_var_1d(ncfileid,'northward_wind@EffectiveQC',nobs_curr,effectiveQC_V)
             endif
                
             if ( trim(adjustl(ob_platform)) .eq. 'gnssro') then
                this_variable = 'occulting_sat_id@MetaData' ! GNSS satellite ID
             else
                this_variable = trim(adjustl(met_variable))//'@ObsType' ! prepbufr integer code
             endif
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,Observation_Type) ! get observation type
          
            !if ( trim(adjustl(met_variable)).eq.'specific_humidity' ) then
            !   call get_netcdf_var_1d(ncfileid,'air_temperature@ObsValue',nobs_curr,aux_Data)
            !endif

             if (twofiles) then
                this_variable = trim(adjustl(met_variable))//'@hofx'
                call get_netcdf_var_1d(ncfileid2,this_variable,nobs_curr,Obs_Minus_Forecast_unadjusted2) ! h(x)
                if ( trim(adjustl(met_variable)).eq.'virtual_temperature' ) then
                   call convert_tv_to_temp('hx',ncfileid2,nobs_curr,Obs_Minus_Forecast_unadjusted2) ! overwrite Obs_Minus_Forecast_unadjusted2
                endif
                Obs_Minus_Forecast_unadjusted2 = observation - Obs_Minus_Forecast_unadjusted2 ! y - h(x)

                if ( trim(adjustl(met_variable)).eq.'eastward_wind' ) then
                   call get_netcdf_var_1d(ncfileid2,'northward_wind@hofx',nobs_curr,Obs_Minus_Forecast_unadjustedV2)
                   Obs_Minus_Forecast_unadjustedV2 = observationV - Obs_Minus_Forecast_unadjustedV2 ! y - h(x)
                endif
                
                if ( trim(adjustl(ob_platform)) .eq. 'gnssro') then
                   this_variable = 'occulting_sat_id@MetaData'
                else
                   this_variable = trim(adjustl(met_variable))//'@ObsType'
                endif
                call get_netcdf_var_1d(ncfileid2,this_variable,nobs_curr,Observation_Type2)
               !if ( trim(adjustl(met_variable)).eq.'specific_humidity' ) then
               !   call get_netcdf_var_1d(ncfileid2,'air_temperature@ObsValue',nobs_curr,aux_Data2)
               !endif
             else
                Obs_Minus_Forecast_unadjusted2 = Obs_Minus_Forecast_unadjusted
                Observation_Type2 = Observation_Type
                if ( trim(adjustl(met_variable)).eq.'eastward_wind' ) Obs_Minus_Forecast_unadjustedV2 = Obs_Minus_Forecast_unadjustedV
               !aux_Data2 = aux_Data
             endif

             do i = 1,nobs_curr
                obmax = abs(observation(i))
                if ( trim(adjustl(met_variable)).eq.'eastward_wind' ) then
                   obmax = max(abs(observation(i)), abs(observationV(i)))
                   effectiveQC(i) = max( effectiveQC(i), effectiveQC_V(i))
                endif
                if (effectiveQC(i) > 0 ) cycle
                if ( trim(adjustl(met_variable)).eq.'surface_pressure' ) then
                   pres = observation(i) * 0.01 ! Pa --> hPa
                else
                   pres = pressure(i) ! already in hPa
                endif
                if ( pres < 0.001_r_kind .or. pres > 1200._r_kind .or. &
                         abs(obmax) > 1.e9_r_kind ) cycle

                ! checking bufr-type/lat/lon/time/pressure for consistency
                if (twofiles) then
                   if(Observation_Type(i) /= Observation_Type2(i) .or. abs(latitude(i)-latitude2(i)) .gt. tolerance .or. &
                     abs(longitude(i)-longitude2(i)) .gt. tolerance .or. abs(time(i)-time2(i)) .gt. tolerance .or. &
                     abs(pressure(i)-pressure2(i)) .gt. tolerance ) then
                     write (6,*) ' ob data inconsistency for i = ',i
                     write (6,*) 'obs_type,obs_type2 = ',Observation_Type(i),Observation_Type2(i)
                     write (6,*) 'lat,lat2 = ',latitude(i),latitude2(i)
                     write (6,*) 'long,long2 = ',longitude(i),longitude2(i)
                     write (6,*) 'time,time2 = ',time(i),time2(i)
                     write (6,*) 'pressure,pressure2 = ',pressure(i),pressure2(i)
                     write (6,*) 'id, id2 = ',id,id2
                     call close_netcdf(obsfile,ncfileid)
                     call close_netcdf(obsfile2,ncfileid2)
                     call stop2(-98)
                   end if
                   if ( trim(adjustl(ob_platform)) .eq. 'gnssro') then
                      if(altitude(i) /= altitude2(i)) then
                         write (6,*) ' ob data height inconsistency for GNSS i = ',i
                         write (6,*) 'altitude1, altitude2, = ',altitude(i),altitude2(i)
                         call close_netcdf(obsfile,ncfileid)
                         call close_netcdf(obsfile2,ncfileid2)
                         call stop2(-98)
                      endif
                   endif
                end if
             
                ! If we're this far, the ob is good and can be assimilated.
                nob = nob + 1

                ! prepBUFR integer code
                ! winds are > 200, T,Q,Ps are < 200
                x_code(nob)  = Observation_Type(i)
                if ( trim(adjustl(met_variable)).eq.'eastward_wind' .or. &
                     trim(adjustl(met_variable)).eq.'northward_wind' ) then
                   if ( x_code(nob) < 200 ) x_code(nob) = x_code(nob) + 100
                else if ( trim(adjustl(met_variable)).eq.'air_temperature' .or. &
                          trim(adjustl(met_variable)).eq.'specific_humidity' .or. &
                          trim(adjustl(met_variable)).eq.'virtual_temperature' .or. &
                          trim(adjustl(met_variable)).eq.'surface_pressure' ) then
                   if ( x_code(nob) > 200 ) x_code(nob) = x_code(nob) - 100
                endif

                ! observation location and time
                x_lat(nob)   = latitude(i)
                x_lon(nob)   = longitude(i)
                if ( trim(adjustl(ob_platform)) .eq. 'gnssro') then
                   x_press(nob) = altitude(i) ! should be in meters
                else
                   x_press(nob) = pres ! already in hPa
                endif
                x_time(nob)  = time(i) ! time relative to analysis time (hrs)

                ! observation errors
                ! orginally standard deviation, so get into variance

                !if (errororig > 1.e-5_r_kind) then
                !   x_errorig(nob) = (one/errororig)**2
                !else
                !   x_errorig(nob) = 1.e10_r_kind
                !endif
               !x_errorig(nob) = (initialObsError(i))**2
               !x_err(nob)   = (finalObsError(i))**2 !(one/error)**2

                ! special handling of gps error
               !if (obtype == 'gps' .and. x_errorig(nob) .gt. 1.e9) x_errorig(nob)=x_err(nob)

                ! observation and obs error; if pressure, get everything into hPa
                if ( trim(adjustl(met_variable)).eq.'surface_pressure' ) then
                   x_obs(nob) = observation(i) * 0.01 ! Pa --> hPa
                   x_errorig(nob) = (initialObsError(i)*0.01)**2
                   x_err(nob)   = (finalObsError(i)*0.01)**2 !(one/error)**2
                  !Obs_Minus_Forecast_adjusted(i) = Obs_Minus_Forecast_adjusted(i) * 0.01
                   Obs_Minus_Forecast_unadjusted(i) = Obs_Minus_Forecast_unadjusted(i) * 0.01
                   Obs_Minus_Forecast_unadjusted2(i) = Obs_Minus_Forecast_unadjusted2(i) * 0.01
                else
                   x_obs(nob)   = observation(i)
                   x_errorig(nob) = (initialObsError(i))**2
                   x_err(nob)   = (finalObsError(i))**2 !(one/error)**2
                endif

                ! hx and hxnobc
               !h_x_ensmean(nob) = x_obs(nob) - Obs_Minus_Forecast_adjusted(i)  ! probably need this eventually
                h_x_ensmean(nob) = x_obs(nob) - Obs_Minus_Forecast_unadjusted(i)

                ! ???
                !if (obtype == '  q' .or. obtype == 'spd' .or. obtype == ' dw' .or. &
                 ! obtype == ' pw') then
                 !  hx_mean_nobc(nob) = hx_mean(nob)
                !endif

                ! observation type
                if ( trim(adjustl(met_variable)).eq.'virtual_temperature' ) then
                   x_type(nob)  = 'air_temperature'
                else
                   x_type(nob)  = trim(adjustl(met_variable)) !obtype character
                endif
                if (x_type(nob) == ' uv')  x_type(nob) = '  u'
                if (x_type(nob) == 'tcp')  x_type(nob) = ' ps'
                if (x_type(nob) == ' rw')  x_type(nob) = '  u'

                ! get Hx for the member
               !if (nanal <= nanals) then
                h_xnobc(nob) = x_obs(nob) - Obs_Minus_Forecast_unadjusted2(i)
                
                ! if U-wind, also output V-wind
                if ( trim(adjustl(met_variable)).eq.'eastward_wind') then
                   nob = nob + 1
                   x_code(nob)  = x_code(nob-1)
                   x_lat(nob)   = x_lat(nob-1)
                   x_lon(nob)   = x_lon(nob-1)
                   x_press(nob) = x_press(nob-1)
                   x_time(nob)  = x_time(nob-1)
                   x_obs(nob)   = observationV(i)
                   x_errorig(nob) = x_errorig(nob-1)
                   x_err(nob)   = x_err(nob-1)
                   h_x_ensmean(nob) = x_obs(nob) - Obs_Minus_Forecast_unadjustedV(i)
                   x_type(nob)  = 'northward_wind' !obtype character
                   h_xnobc(nob) = x_obs(nob) - Obs_Minus_Forecast_unadjustedV2(i)
                endif
                
                ! Transform specific humidity to RH
               !if ( trim(adjustl(met_variable)).eq.'specific_humidity' ) then
               !   ges_prsl(1,1) = x_press(nob) ! needs to be hPa
               !   ges_tv(1,1) = aux_Data(i)*(one+fv*x_obs(nob))
               !   ! call genqsat1(sph,qsat,ges_prsl,ges_tv,ice,npts,nlevs)
               !   call genqsat1(x_obs(nob),qsat,ges_prsl,ges_tv,.true.,1,1)
               !   x_obs(nob)   = observation(i)/qsat(1,1)
               !   x_errorig(nob) = (initialObsError(i)/qsat(1,1))**2
               !   x_err(nob)   = (finalObsError(i)/qsat(1,1))**2 !(one/error)**2
               !   h_x_ensmean(nob) = h_x_ensmean(nob)/qsat(1,1)
               !   h_xnobc(nob) = h_xnobc(nob)/qsat(1,1)
               !endif
               
                  !if (obtype == '  q' .or. obtype == 'spd' .or. obtype == ' dw' .or. obtype == ' pw') then
                      !hx(nob) = Observation(i) - Obs_Minus_Forecast_adjusted2(i)
                  !endif
               !endif

             enddo ! loop over "i" (number of obs)

          enddo ! loop over "j" (number of obs for this specific platform and meteorological type)

          ! clean-up
          call close_netcdf(obsfile,ncfileid)
          if (twofiles) call close_netcdf(obsfile2,ncfileid2)

          deallocate(latitude, longitude, pressure, time, altitude)
          deallocate(effectiveQC, initialObsError, finalObsError,effectiveQC_V)
          deallocate(Observation_Type, observation)
          deallocate(Obs_Minus_Forecast_adjusted, Obs_Minus_Forecast_unadjusted)
          deallocate(jedi_time_string)
          deallocate(aux_Data,aux_Data2)
          deallocate(Observation_Type2, Obs_Minus_Forecast_adjusted2, Obs_Minus_Forecast_unadjusted2)
          deallocate(observationV,Obs_Minus_Forecast_unadjustedV,Obs_Minus_Forecast_unadjustedV2)
          if (twofiles) deallocate(latitude2, longitude2, pressure2, time2, jedi_time_string2,altitude2)

       enddo peloop ! end loop over "ipe"
    enddo obtypeloop ! end loop over "itype"

    if (nob .ne. nobs_max) then
       print *,'number of obs not what expected in get_convobs_data',nob,nobs_max
       call stop2(94)
    end if

end subroutine get_convobs_data

subroutine convert_tv_to_temp(what_for,ncid,nobs_curr,tv,effectiveQC)
   character(len=*), intent(in) :: what_for
   integer(i_kind), intent(in) :: ncid, nobs_curr
   real(r_single), dimension(nobs_curr),  intent(inout) :: tv ! input is virtual temperature. output is sensible temperature
   integer(i_kind), dimension(nobs_curr), intent(inout), optional :: effectiveQC
   
   integer(i_kind) :: i
   real(r_single), allocatable, dimension (:) :: q, q_qc, q_err

   allocate(q(nobs_curr), q_qc(nobs_curr), q_err(nobs_curr))

   ! get specific humidity (kg/kg). We can use this routine to convert
   !  observation values and H(x) values
   if ( trim(adjustl(what_for)) == 'obs' ) then
      call get_netcdf_var_1d(ncid,'specific_humidity@ObsValue',nobs_curr,q)
   else if ( trim(adjustl(what_for)) == 'hx' ) then
      call get_netcdf_var_1d(ncid,'specific_humidity@hofx',nobs_curr,q)
   endif

   ! convert Tv to T. RHS is Tv, LHS is T.  Deal with QC next.
   tv = tv / (1.0 + 0.61 * q)  ! formula in GSI prepbufr processing

   ! Now deal with QC.  If it's a Tv ob, we need to assume the Q ob is good, but we'll do some gross error checking
   if ( present(effectiveQC)) then
      call get_netcdf_var_1d(ncid,'specific_humidity@EffectiveError',nobs_curr,q_err) ! not used yet
      call get_netcdf_var_1d(ncid,'specific_humidity@EffectiveQC',nobs_curr,q_qc)
      do i = 1,nobs_curr
         if (effectiveQC(i) > 0 ) cycle ! effectiveQC(i) is initially from Tv--if already bad, nothing to do.
      !!!if ( q_qc(i) > 0 ) effectiveQC(i) = 100 ! If QC for Q is bad, reset effective QC
         if ( q(i) < 0.0 .or. q(i) > 0.05 ) effectiveQC(i) = 100 ! Quick gross-error check on Q 
         if ( tv(i) < 100.0 .or. tv(i) > 350.0 ) effectiveQC(i) = 100 ! Quick gross-error check on Tv
      enddo
   endif

   deallocate(q, q_qc, q_err)
end subroutine convert_tv_to_temp

subroutine check_met_variables_for_this_platform(platform,met_variable,skip)
   character(len=*), intent(in) :: platform, met_variable
   logical, intent(out) :: skip

   character(len=100), dimension(num_met_variables):: vars_to_keep
   integer(i_kind) :: i

   vars_to_keep(:) = '' ! initialize to blank

   if ( trim(adjustl(platform)) .eq. 'sondes' ) then
     if ( process_sondes_Tv ) then
        vars_to_keep(1:5) = (/ 'air_temperature', 'eastward_wind', 'northward_wind', 'specific_humidity', 'virtual_temperature' /)
     else
        vars_to_keep(1:4) = (/ 'air_temperature', 'eastward_wind', 'northward_wind', 'specific_humidity' /)
     endif
   else if ( trim(adjustl(platform)) .eq. 'aircraft' ) then
      vars_to_keep(1:4) = (/ 'air_temperature', 'eastward_wind', 'northward_wind', 'specific_humidity' /)
   else if ( trim(adjustl(platform)) .eq. 'gnssro'   ) then
     !vars_to_keep(1:2) = (/ 'refractivity', 'bending_angle' /)
      vars_to_keep(1) = 'refractivity'
   else if ( trim(adjustl(platform)) .eq. 'satwind' ) then
      vars_to_keep(1:2) = (/ 'eastward_wind', 'northward_wind' /)
   else if ( trim(adjustl(platform)) .eq. 'sfc'   ) then
      vars_to_keep(1) = 'surface_pressure'
   endif

   ! make sure everything is trimmed for following if/then
   do i = 1,num_met_variables
      vars_to_keep(i) = trim(adjustl(vars_to_keep(i)))
   enddo

   if( any( trim(adjustl(met_variable)).eq.vars_to_keep) ) then
      skip = .false.  ! don't skip this met_variable
   else
      skip = .true.   ! skip this met_variable
   endif

   return

end subroutine check_met_variables_for_this_platform

end module readconvobs
