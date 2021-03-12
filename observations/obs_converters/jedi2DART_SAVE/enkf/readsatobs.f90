module readsatobs
!$$$  module documentation block
!
! module: readsatobs                   read data from satellite radiance
!                                      diag* files.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract: read data from satellite radiance diag* files written out
!  by GSI forward operator code.
!
! Public Subroutines:
!  get_num_satobs: determine the number of observations to read.
!  get_satobs_data: read the data.
!   
! program history log:
!   2009-02-23  Initial version.
!
! attributes:
!   language: f95
!
!$$$
  
use kinds, only: r_kind,i_kind,r_single
use params, only: nsats_rad, nsatmax_rad, dsis, sattypes_rad, num_ufo_procs
use netcdf_mod, only: open_netcdf, close_netcdf, get_netcdf_dims, get_netcdf_var_1d, get_netcdf_var_2d
use mpisetup, only: stop2
use jedi_time_routines_mod, only : read_jedi_time
use constants, only : init_constants_derived, deg2rad
use netcdf

implicit none

private
public :: get_satobs_data, get_num_satobs

contains

subroutine get_num_satobs(obspath,datestring,num_obs_tot,id)

    implicit none

    character (len=500), intent(in) :: obspath
    character(len=10), intent(in) :: id, datestring
    integer(i_kind), intent(out) :: num_obs_tot

    character(len=500) obsfile, geofile, ydiagfile
    character(len=4) :: pe_name, ichan
    character(len=100) :: this_variable
    integer(i_kind) :: nsat, ipe, i, j, idx, ncfileid, ncfileid_geo, ncfileid_ydiag
    integer(i_kind) :: nobs_curr, nchans, nlevs, this_channel, num_obs_totdiag
    logical :: fexist, skip
    real(r_kind) :: errorlimit,errorlimit2

    integer(i_kind), allocatable, dimension(:) :: effectiveQC, channels
    integer(i_kind), allocatable, dimension(:) :: nobs
    real(r_single), allocatable, dimension(:) :: pressure, finalObsError, observation
    real(r_single), allocatable, dimension(:,:) :: pressure2d
    integer(i_kind), allocatable, dimension(:,:) :: model_level_at_peak_of_weightingfunction

!  make consistent with screenobs
    errorlimit=1._r_kind/sqrt(1.e9_r_kind)
    errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)
    
    allocate(nobs(nsats_rad))
    nobs = 0
    num_obs_tot = 0
    num_obs_totdiag = 0

    do nsat=1,nsats_rad

       peloop: do ipe=0,num_ufo_procs-1 ! PEs start at 0, hence the -1

          write(pe_name,'(i4.4)') ipe

          ! Need to incorporate $id into filename
          if ( num_ufo_procs == 1 ) then
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_all_'//trim(adjustl(id))//'.nc4'
          else
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
             geofile = trim(adjustl(obspath))//'/geovals_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
             ydiagfile = trim(adjustl(obspath))//'/ydiag_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
          endif
          !write(*,fmt='(a)')'obsfile = ',obsfile

          ! Make sure files exist
          inquire(file=obsfile,exist=fexist) ; if (.not.fexist) cycle peloop
          inquire(file=geofile,exist=fexist) ; if (.not.fexist) cycle peloop
          inquire(file=ydiagfile,exist=fexist) ; if (.not.fexist) cycle peloop

          ! File exists.  Open it, figure out how many obs this file has.
          call open_netcdf(obsfile,ncfileid)
          call get_netcdf_dims(ncfileid,'nlocs',nobs_curr)
          if ( nobs_curr .le. 0 ) then
             call close_netcdf(obsfile,ncfileid)
             cycle peloop ! If the file has no obs go to next file
          endif
          !write(*,*)' nlocs = ',nobs_curr

          allocate(effectiveQC(nobs_curr),pressure(nobs_curr))
          allocate(finalObsError(nobs_curr), observation(nobs_curr))
          allocate( model_level_at_peak_of_weightingfunction(1,nobs_curr))

          ! Get channels and number of channels
          call get_netcdf_dims(ncfileid,'nvars',nchans)
          !write(*,*)' there are  ',nchans,' channels'
          allocate(channels(nchans))
          !call get_netcdf_var_1d(ncfileid,'sensor_channel@VarMetaData',nchans,channels)
          call get_channels_from_var_names(ncfileid,nchans,channels)

          ! get pressure profile at each location; needed for pressure at peak of weighting function
          ! this opens geofile, closes geofile, and allocates pressure2d
          call get_pressure_profiles(geofile,nobs_curr,pressure2d)

          call open_netcdf(ydiagfile,ncfileid_ydiag) ! open ydiag file, don't do anything yet with it

          do j = 1,nchans

             this_channel = channels(j)
             write(ichan,'(i4)') this_channel

             ! Figure out if we want obs from this channel
             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@EffectiveQC'
             call check_channels(ncfileid,this_variable,skip)
             if ( skip ) cycle
             if ( ipe == 0 ) write(*,*)'keeping channel ',ichan,' for ',trim(adjustl(sattypes_rad(nsat)))
                
             ! If we're this far, we want obs from the channel
             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@ObsValue'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,observation)

             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@EffectiveQC' ! this may be all that's needed
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,effectiveQC)

             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@EffectiveError'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,finalObsError)

             this_variable = 'pressure_level_at_peak_of_weightingfunction_'//trim(adjustl(ichan)) ! 'pressure_level_at_peak_of_weightingfunction_7'
             call get_netcdf_var_2d(ncfileid_ydiag,this_variable,1,nobs_curr,model_level_at_peak_of_weightingfunction) ! Model level at peak of weighting function

             do i = 1,nobs_curr
                num_obs_totdiag = num_obs_totdiag + 1 ! total number of obs read...not necessarily total number used
                
                idx = model_level_at_peak_of_weightingfunction(1,i) ! index is top down
                pressure(i) = pressure2d(idx,i) ! pressure at peak of weighting function (hPa)

                if (effectiveQC(i) > 0 ) cycle
                if ( pressure(i) < 0.001_r_kind .or. pressure(i) > 1200._r_kind .or. &
                         abs(observation(i)) > 1.e9_r_kind ) cycle
                nobs(nsat) = nobs(nsat) + 1  ! number of "kept" observations for this satellite
                num_obs_tot = num_obs_tot + 1 ! number of total radiance obs across all satellites/channels
             enddo
          
          enddo ! loop over "j" (number of channels)

          ! clean-up
          call close_netcdf(obsfile,ncfileid)
          call close_netcdf(ydiagfile,ncfileid_ydiag)
          deallocate(pressure, effectiveQC, finalObsError, observation, channels)
          deallocate(pressure2d, model_level_at_peak_of_weightingfunction)

       enddo peloop ! end loop over "ipe"
       
    enddo ! end loop over satellites ("nsat")

    write(6,*),num_obs_totdiag, ' radiance obs read from all files'
    write(6,*),num_obs_tot,' radiance obs actually kept from all files'
    do nsat = 1,nsats_rad
       write(6,100) nsat,trim(adjustl(sattypes_rad(nsat))),nobs(nsat)
    enddo
100    format(2x,i3,2x,a20,2x,'num_obs_tot= ',i9)
    
    deallocate(nobs)

end subroutine get_num_satobs

subroutine get_satobs_data(obspath, datestring, nobs_max, h_x, h_xnobc, x_obs, x_err, &
           x_lon, x_lat, x_press, x_time, x_channum, x_errorig, x_type, x_biaspred, x_indx,id,id2)

! use radinfo, only: iuse_rad,nusis,jpch_rad,nuchan,npred,adp_anglebc,emiss_bc
  use radinfo, only: npred, adp_anglebc,emiss_bc, angord

  implicit none

  character*500, intent(in) :: obspath
  character(len=10), intent(in) ::  datestring
  integer(i_kind), intent(in) :: nobs_max
  character(len=10), intent(in) :: id,id2

  real(r_single), dimension(nobs_max) :: h_x,h_xnobc,x_obs,x_err,x_lon,&
                               x_lat,x_press,x_time,x_errorig
  real(r_single), dimension(npred+1,nobs_max) :: x_biaspred
  integer(i_kind), dimension(nobs_max) ::  x_channum,x_indx
  character(len=20), dimension(nobs_max) ::  x_type

  character(len=500) :: obsfile,obsfile2, geofile, ydiagfile
  character(len=20)  :: sat_type
  character(len=4)   :: pe_name, ichan
  character(len=100) :: this_variable
  integer(i_kind) :: nob, nsat, ipe, i, j, idx, nn, ii, nlevs
  integer(i_kind) :: ncfileid, ncfileid2, ncfileid_geo, ncfileid_ydiag
  integer(i_kind) :: nobs_curr, nobs_curr2, nchans, nchans2, this_channel
  logical         :: twofiles,fexist,skip
  real(r_kind)    :: errorlimit,errorlimit2
  real(r_kind) :: tolerance = 1.e-5

  integer(i_kind), allocatable, dimension (:) :: gsiUseFlag, effectiveQC, channels
  real(r_single), allocatable, dimension (:) :: initialObsError, finalObsError
  real(r_single), allocatable, dimension (:) :: observation, obsBias
  real(r_single), allocatable, dimension (:) :: latitude, longitude, pressure, time
  real(r_single), allocatable, dimension (:) :: latitude2, longitude2, time2
  real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_adjusted
  real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_unadjusted
  real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_adjusted2
  real(r_single), allocatable, dimension (:) :: Obs_Minus_Forecast_unadjusted2
  character(len=20), allocatable, dimension(:) :: jedi_time_string, jedi_time_string2

  real(r_single), dimension(:), allocatable :: BC_Fixed_Scan_Position, BCPred_Constant, BCPred_Scan_Angle
  real(r_single), dimension(:), allocatable :: BCPred_Cloud_Liquid_Water, BCPred_Lapse_Rate_Squared, BCPred_Lapse_Rate
  real(r_single), dimension(:), allocatable :: BCPred_Cosine_Latitude_times_Node, BCPred_Sine_Latitude
  real(r_single), dimension(:), allocatable :: BCPred_Emissivity, scanAngle
  real(r_single), allocatable, dimension (:,:) :: BCPred_angord
  real(r_single), allocatable, dimension(:,:) :: pressure2d
  integer(i_kind), allocatable, dimension(:,:) :: model_level_at_peak_of_weightingfunction

  ! make consistent with screenobs
    errorlimit=1._r_kind/sqrt(1.e9_r_kind)
    errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)

    twofiles = id /= id2

    nob = 0

    call init_constants_derived ! set deg2rad

    do nsat=1,nsats_rad

       peloop: do ipe=0,num_ufo_procs-1 ! PEs start at 0, hence the -1

          write(pe_name,'(i4.4)') ipe

          ! Need to incorporate $id into filename
          if ( num_ufo_procs == 1 ) then
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_all_'//trim(adjustl(id))//'.nc4'
          else
             obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
             geofile = trim(adjustl(obspath))//'/geovals_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
             ydiagfile = trim(adjustl(obspath))//'/ydiag_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_'//pe_name//'_'//trim(adjustl(id))//'.nc4'
          endif
          !write(*,fmt='(a)')'obsfile = ',obsfile

          ! Make sure file exists
          inquire(file=obsfile,exist=fexist) ; if (.not.fexist) cycle peloop
          inquire(file=geofile,exist=fexist) ; if (.not.fexist) cycle peloop
          inquire(file=ydiagfile,exist=fexist) ; if (.not.fexist) cycle peloop

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
                obsfile2 = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_all_'//trim(adjustl(id2))//'.nc4'
             else
                obsfile2 = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(sattypes_rad(nsat)))//'_'//pe_name//'_'//trim(adjustl(id2))//'.nc4'
             endif

             ! Make sure files exists
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

          ! allocate lots of variables
          allocate(latitude(nobs_curr), longitude(nobs_curr), pressure(nobs_curr), time(nobs_curr), &
                 gsiUseFlag(nobs_curr), effectiveQC(nobs_curr), initialObsError(nobs_curr), finalObsError(nobs_curr), &
                 observation(nobs_curr), jedi_time_string(nobs_curr), obsBias(nobs_curr), &
                 Obs_Minus_Forecast_adjusted(nobs_curr), Obs_Minus_Forecast_unadjusted(nobs_curr))
          allocate(Obs_Minus_Forecast_adjusted2(nobs_curr), Obs_Minus_Forecast_unadjusted2(nobs_curr))
          allocate( model_level_at_peak_of_weightingfunction(1,nobs_curr))

          ! bias correction variables
          allocate ( BCPred_Constant(nobs_curr), BCPred_Scan_Angle(nobs_curr), BCPred_Cloud_Liquid_Water(nobs_curr))
          allocate ( BCPred_Lapse_Rate_Squared(nobs_curr), BCPred_Lapse_Rate(nobs_curr), BCPred_Cosine_Latitude_times_Node(nobs_curr))
          allocate ( BCPred_Sine_Latitude(nobs_curr), BCPred_Emissivity(nobs_curr))
          allocate ( BC_Fixed_Scan_Position(nobs_curr)) 
          allocate ( scanAngle(nobs_curr), BCPred_angord(angord, nobs_curr))

          ! Get lat/lon/time/pressure from ensemble mean
          call get_netcdf_var_1d(ncfileid,'latitude@MetaData',nobs_curr,latitude)
          call get_netcdf_var_1d(ncfileid,'longitude@MetaData',nobs_curr,longitude)
          call get_netcdf_var_1d(ncfileid,'datetime@MetaData',nobs_curr,jedi_time_string)
          ! get difference between analysis time and observation times (hrs)
          call read_jedi_time(datestring,nobs_curr,jedi_time_string,time) ! output is time

          ! Get channels
          call get_netcdf_dims(ncfileid,'nvars',nchans)
          !write(*,*)' there are  ',nchans,' channels'
          allocate(channels(nchans))
         !call get_netcdf_var_1d(ncfileid,'sensor_channel@VarMetaData',nchans,channels)
          call get_channels_from_var_names(ncfileid,nchans,channels)

          ! get pressure profile at each location; needed for pressure at peak of weighting function
          ! this opens geofile, closes geofile, and allocates pressure2d
          call get_pressure_profiles(geofile,nobs_curr,pressure2d)

          if (twofiles) then
             allocate(latitude2(nobs_curr), longitude2(nobs_curr), time2(nobs_curr))
             allocate(jedi_time_string2(nobs_curr))
             call get_netcdf_var_1d(ncfileid2,'latitude@MetaData',nobs_curr,latitude2)
             call get_netcdf_var_1d(ncfileid2,'longitude@MetaData',nobs_curr,longitude2)
             call get_netcdf_var_1d(ncfileid2,'datetime@MetaData',nobs_curr,jedi_time_string2)
             ! get difference between analysis time and observation times (hrs)
             call read_jedi_time(datestring,nobs_curr,jedi_time_string2,time2) ! output is time2

             call get_netcdf_dims(ncfileid2,'nvars',nchans2)
             if ( nchans /= nchans2 ) then
                write (6,*) 'mismatching number of channels'
                write (6,*) 'nchannels,nchannels2 = ',nchans,nchans2
                call close_netcdf(obsfile,ncfileid)
                call close_netcdf(obsfile2,ncfileid2)
                call stop2(-98)
             endif
          endif

          call open_netcdf(ydiagfile,ncfileid_ydiag) ! open ydiag file, don't do anything yet with it

          do j = 1,nchans

             this_channel = channels(j)
             write(ichan,'(i4)') this_channel

             ! Figure out if we want obs from this channel
             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@EffectiveQC'
             call check_channels(ncfileid,this_variable,skip)
             if ( skip ) cycle

             ! If we're this far, we want obs from the channel
             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@ObsValue'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,observation)

             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@EffectiveQC'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,effectiveQC)

             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@EffectiveError'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,finalObsError)

             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@ObsError'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,initialObsError)

             this_variable = 'pressure_level_at_peak_of_weightingfunction_'//trim(adjustl(ichan)) ! 'pressure_level_at_peak_of_weightingfunction_7'
             call get_netcdf_var_2d(ncfileid_ydiag,this_variable,1,nobs_curr,model_level_at_peak_of_weightingfunction) ! Model level at peak of weighting function

            !this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@depbg' ! make sure this is y-h(x) and not h(x) -y
            !call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,Obs_Minus_Forecast_adjusted) ! y-h(x)

             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@ObsBias'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,obsBias) ! amount of total bias correction

             this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@hofx'
             call get_netcdf_var_1d(ncfileid,this_variable,nobs_curr,Obs_Minus_Forecast_unadjusted) ! h(x)
             Obs_Minus_Forecast_unadjusted = observation - Obs_Minus_Forecast_unadjusted ! y - h(x)

             if (twofiles) then
                this_variable = 'brightness_temperature_'//trim(adjustl(ichan))//'@hofx'
                call get_netcdf_var_1d(ncfileid2,this_variable,nobs_curr,Obs_Minus_Forecast_unadjusted2) ! h(x)
                Obs_Minus_Forecast_unadjusted2 = observation - Obs_Minus_Forecast_unadjusted2 ! y - h(x)
             else
                Obs_Minus_Forecast_unadjusted2 = Obs_Minus_Forecast_unadjusted ! y - h(x)
             endif
             
            ! bias correction variables from the ensemble mean
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@constantPredictor',nobs_curr,BCPred_Constant)
           !call get_netcdf_var_1d(ncfileid,'scan_angle_bias_correction_term@TestReference',nobs_curr,BCPred_Scan_Angle)
            BCPred_Scan_Angle = 0.0
           !call get_netcdf_var_1d(ncfileid,'cloud_liquid_water_bias_correction_term@TestReference',nobs_curr,BCPred_Cloud_Liquid_Water)
            BCPred_Cloud_Liquid_Water = 0.0
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@lapse_rate_order_2Predictor',nobs_curr,BCPred_Lapse_Rate_Squared)
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@lapse_ratePredictor',nobs_curr,BCPred_Lapse_Rate)
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@cosine_of_latitude_times_orbit_nodePredictor',nobs_curr,BCPred_Cosine_Latitude_times_Node)
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@sine_of_latitudePredictor',nobs_curr,BCPred_Sine_Latitude)
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@emissivityPredictor',nobs_curr,BCPred_Emissivity)
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@scan_angle_order_4Predictor',nobs_curr,BCPred_angord(1,:))
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@scan_angle_order_3Predictor',nobs_curr,BCPred_angord(2,:))
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@scan_angle_order_2Predictor',nobs_curr,BCPred_angord(3,:))
            call get_netcdf_var_1d(ncfileid,'brightness_temperature_'//trim(adjustl(ichan))//'@scan_anglePredictor',nobs_curr,BCPred_angord(4,:))

            ! Deal with scan angle. Not sure 'sensor_view_angle@MetaData' is correct.
           !call get_netcdf_var_1d(ncfileid,'sensor_view_angle@MetaData',nobs_curr,scanAngle)
           !scanAngle = scanAngle * deg2rad

            ! This variable is the actual bias correction [sum(beta*predictor)] for the scan angle
            !    This will need to come in through UFO where UFO has access to satbias_in to get beta
            BC_Fixed_Scan_Position(:) = 1.0 ! temporary hack
            !call get_netcdf_var_1d(ncfileid,'BC_Fixed_Scan_Position',nobs_curr,BC_Fixed_Scan_Position)

             do i = 1,nobs_curr

                idx = model_level_at_peak_of_weightingfunction(1,i) ! index is top down
                pressure(i) = pressure2d(idx,i) ! pressure at peak of weighting function (hPa)
                
               !if ( .not. twofiles )write(*,*)'qc,emiss = ',effectiveQC(i), BCPred_Emissivity(i)
                if (effectiveQC(i) > 0 ) cycle
                if ( pressure(i) < 0.001_r_kind .or. pressure(i) > 1200._r_kind .or. &
                         abs(observation(i)) > 1.e9_r_kind ) cycle

                ! checking lat/lon/time/pressure for consistency
                if (twofiles) then
                   if( abs(latitude(i)-latitude2(i)) .gt. tolerance .or. &
                     abs(longitude(i)-longitude2(i)) .gt. tolerance .or. &
                     abs(time(i)-time2(i)) .gt. tolerance ) then
                     write (6,*) ' ob data inconsistency '
                     write (6,*) 'lat,lat2 = ',latitude(i),latitude2(i)
                     write (6,*) 'long,long2 = ',longitude(i),longitude2(i)
                     write (6,*) 'time,time2 = ',time(i),time2(i)
                     call close_netcdf(obsfile,ncfileid)
                     call close_netcdf(obsfile2,ncfileid2)
                     call stop2(-98)
                   end if
                end if

                ! If we're this far, the ob is good and can be assimilated.
                nob = nob + 1 

                write(sat_type,'(a20)') trim(adjustl(sattypes_rad(nsat)))
                x_type(nob)= sat_type
                x_channum(nob) = this_channel
                x_indx(nob) = this_channel ! originally was a mapping to entries in satinfo file ! indxsat
                x_lon(nob) = longitude(i)
                x_lat(nob) = latitude(i)
                x_time(nob) = time(i)
                x_press(nob) = pressure(i) ! hPa
                x_obs(nob) = observation(i)

                ! bias corrected Hx based on ensemble mean
                h_x(nob) = (x_obs(nob) - obsBias(i)) - Obs_Minus_Forecast_unadjusted(i)

                ! un-bias corrected Hx for the member
                h_xnobc(nob) = x_obs(nob) - Obs_Minus_Forecast_unadjusted2(i)

                ! orginally standard deviation, so get into variance
                x_errorig(nob) = (initialObsError(i))**2
                x_err(nob) = (finalObsError(i))**2

        !! DTK:  **NOTE**
       !!       The bifix term will need to be expanded if/when the GSI/GDAS goes to using
       !!       a higher polynomial version of the angle dependent bias correction (if
       !!         and when it is moved into part of the varbc)
       !!         x_biaspred(1,nob) = data_chan1(n)%bifix! fixed angle dependent bias

                !x_biaspred(1,nob) = data_chan1(n)%bifix(1) ! fixed angle dependent bias
                !x_biaspred(2,nob) = data_chan1(n)%bicons ! constant bias correction
                !x_biaspred(3,nob) = data_chan1(n)%biang ! scan angle bias correction
                !x_biaspred(4,nob) = data_chan1(n)%biclw ! CLW bias correction
                !x_biaspred(5,nob) = data_chan1(n)%bilap2 ! square lapse rate bias corr
                !x_biaspred(6,nob) = data_chan1(n)%bilap ! lapse rate bias correction
                !if (npred == 7) then
                !   x_biaspred(7,nob) = data_chan1(n)%bicos ! node*cos(lat) bias correction for SSMIS
                !   x_biaspred(8,nob) = data_chan1(n)%bisin ! sin(lat) bias correction for SSMIS                    
                !endif
                !if (emiss_bc) x_biaspred(9,nob) = data_chan1(n)%biemis

                !if (adp_anglebc) then
                !   x_biaspred( 1,nob)  = data_chan1(n)%bifix(5) ! fixed angle dependent bias correction
                !   x_biaspred(npred-2,nob)  = data_chan1(n)%bifix(1) ! 4th order scan angle (predictor)
                !   x_biaspred(npred-1,nob)  = data_chan1(n)%bifix(2) ! 3rd order scan angle (predictor)
                !   x_biaspred(npred,nob)  = data_chan1(n)%bifix(3) ! 2nd order scan angle (predictor)
                !   x_biaspred(npred+1,nob)    = data_chan1(n)%bifix(4) ! 1st order scan angle (predictor)
                !endif

                !x_biaspred(:,nob) = 0. ! if no bias-correction update in EnKF, no need for predictors

! from radinfo: radiance bias correction terms are as follows:
!  pred(1,:)  = global offset
!  pred(2,:)  = zenith angle predictor, is not used and set to zero now
!  pred(3,:)  = cloud liquid water predictor for clear-sky microwave radiance assimilation
!  pred(4,:)  = square of temperature laps rate predictor
!  pred(5,:)  = temperature laps rate predictor
!  pred(6,:)  = cosinusoidal predictor for SSMI/S ascending/descending bias
!  pred(7,:)  = sinusoidal predictor for SSMI/S
!  pred(8,:)  = emissivity sensitivity predictor for land/sea differences
!  pred(9,:)  = fourth order polynomial of angle bias correction
!  pred(10,:) = third order polynomial of angle bias correction
!  pred(11,:) = second order polynomial of angle bias correction
!  pred(12,:) = first order polynomial of angle bias correction
                !if (lupd_satbiasc) then ! bias predictors only used if lupd_satbiasc=T
            
                    x_biaspred(1,nob) = BCPred_Constant(i) !data_chan1(n)%bicons ! constant bias correction
                    x_biaspred(2,nob) = BCPred_Scan_Angle(i) !data_chan1(n)%biang ! scan angle bias correction
                    x_biaspred(3,nob) = BCPred_Cloud_Liquid_Water(i) !data_chan1(n)%biclw ! CLW bias correction
                    x_biaspred(4,nob) = BCPred_Lapse_Rate_Squared(i) !data_chan1(n)%bilap2 ! square lapse rate bias corr
                    x_biaspred(5,nob) = BCPred_Lapse_Rate(i) !data_chan1(n)%bilap ! lapse rate bias correction
                    x_biaspred(6,nob) = BCPred_Cosine_Latitude_times_Node(i) !data_chan1(n)%bicos ! node*cos(lat) bias correction for SSMIS
                    x_biaspred(7,nob) = BCPred_Sine_Latitude(i) !data_chan1(n)%bisin ! sin(lat) bias correction for SSMIS                    
                    if (emiss_bc) then
                       x_biaspred(8,nob) = BCPred_Emissivity(i) !data_chan1(n)%biemis
                       nn = 9
                    else
                       nn = 8
                    endif

                    ! Can simplify to get rid of variable BCPred_angord
                   !do ii = 1,angord
                   !   BCPred_angord(angord-ii+1,i) = scanAngle(i)**ii
                   !enddo
                    if (adp_anglebc) then
                       x_biaspred(nn  ,nob)  = BCPred_angord(1,i) !data_chan1(n)%bifix(1) ! 4th order scan angle (predictor)
                       x_biaspred(nn+1,nob)  = BCPred_angord(2,i) ! data_chan1(n)%bifix(2) ! 3rd order scan angle (predictor)
                       x_biaspred(nn+2,nob)  = BCPred_angord(3,i) !data_chan1(n)%bifix(3) ! 2nd order scan angle (predictor)
                       x_biaspred(nn+3,nob)  = BCPred_angord(4,i) !data_chan1(n)%bifix(4) ! 1st order scan angle (predictor)
                    endif
                    x_biaspred(npred+1,nob) = BC_Fixed_Scan_Position(i) !data_chan1(n)%bifix(5) ! CSS bug fix. fixed angle dependent bias
                !else
                !   x_biaspred(:,nob) =zero ! lupd_satbiasc=F, don't need bias predictors
                !endif

             enddo  ! loop over "i" (number of obs)

          enddo ! loop over "j" (number of channels)
          
          ! clean-up
          call close_netcdf(obsfile,ncfileid)
          call close_netcdf(ydiagfile,ncfileid_ydiag)
          if (twofiles) call close_netcdf(obsfile2,ncfileid2)

          deallocate(latitude, longitude, pressure, time)
          deallocate(gsiUseFlag, effectiveQC, initialObsError, finalObsError)
          deallocate(observation,jedi_time_string, obsBias)
          deallocate(Obs_Minus_Forecast_adjusted, Obs_Minus_Forecast_unadjusted)
          deallocate(Obs_Minus_Forecast_adjusted2, Obs_Minus_Forecast_unadjusted2)
          deallocate(channels)
          deallocate(BCPred_Constant,BCPred_Scan_Angle,BCPred_Cloud_Liquid_Water,BCPred_Lapse_Rate_Squared,BCPred_Lapse_Rate)
          deallocate(BCPred_Cosine_Latitude_times_Node,BCPred_Sine_Latitude,BCPred_Emissivity,BC_Fixed_Scan_Position,BCPred_angord)
          deallocate(scanAngle)
          deallocate(pressure2d, model_level_at_peak_of_weightingfunction)
          if(twofiles) deallocate(latitude2, longitude2, time2, jedi_time_string2)

       enddo peloop ! end loop over "ipe"
    enddo ! end loop over satellites ("nsat")

    if (nob .ne. nobs_max) then
       print *,'number of obs not what expected in get_satobs_data',nob,nobs_max
       call stop2(92)
    end if

end subroutine get_satobs_data

subroutine get_channels_from_var_names(ncfileid,nchans,channels)
   integer(i_kind), intent(in) :: ncfileid, nchans
   integer(i_kind), intent(inout) :: channels(nchans)

   character(len=25) :: variable_names(nchans)
   integer(i_kind) :: j, ipos

   call get_netcdf_var_1d(ncfileid,'variable_names@VarMetaData',nchans,variable_names)
   !!!!!!write(*,*) SCAN("FORTRAN", "AO", .TRUE.)  ! 6, found 'A'
   do j = 1,nchans
      ipos = scan(trim(adjustl(variable_names(j))), "_", .true.) ! get location of last instance of "_"
      read(variable_names(j)(ipos+1:), '(i)') channels(j)
   enddo
end subroutine get_channels_from_var_names

subroutine check_channels(ncfileid,this_variable,skip)
   integer, intent(in)           :: ncfileid
   character(len=*), intent(in)  :: this_variable
   logical, intent(out)          :: skip

   integer(i_kind) :: ncvarid, ncstatus

   ! if variable not there, the channel was not selected for assimilation
   ncstatus = nf90_inq_varid(ncfileid,trim(adjustl(this_variable)),ncvarid)

   if ( ncstatus /= 0 ) then
      skip = .true. ! skip this channel
   else
      skip = .false. ! keep this channel
   endif
             
end subroutine check_channels

subroutine get_pressure_profiles(fname,nobs_curr,pressure2d)
   character(len=*), intent(in) :: fname
   integer(i_kind), intent(in)  :: nobs_curr
   real(r_single), allocatable, dimension(:,:), intent(inout) :: pressure2d

   integer(i_kind) :: ncid, nlevs

   call open_netcdf(fname,ncid) ! open file
   call get_netcdf_dims(ncid,'air_pressure_nval',nlevs)
   allocate(pressure2d(nlevs,nobs_curr))
   call get_netcdf_var_2d(ncid,'air_pressure',nlevs,nobs_curr,pressure2d) ! Pa; 3D model pressure at each ob location
   pressure2d = pressure2d * 0.01 ! get into hPa
   call close_netcdf(fname,ncid)
end subroutine get_pressure_profiles

end module readsatobs
