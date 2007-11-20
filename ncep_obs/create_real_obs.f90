! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_real_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use types_mod,        only : r8, deg2rad, PI
use obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                             static_init_obs_sequence, destroy_obs_sequence 
use     real_obs_mod, only : real_obs_sequence
use    utilities_mod, only : initialize_utilities, register_module, do_output, &
                             error_handler, timestamp, E_ERR, E_MSG, logfileunit, &
                             find_namelist_in_file, check_namelist_read

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq

character(len = 129) :: output_name
character(len = 8 ) :: obsdate
integer :: iunit, io, ii, day1, kkk, kbeg, kend

character(len = 2) :: obstime(4), hour1
data obstime/'06','12','18','24'/

real(r8) :: bin_beg(5), bin_end(5)
data bin_beg/ 3.01_r8,  9.01_r8, 15.01_r8, 21.01_r8, 3.01_r8/
data bin_end/ 9.00_r8, 15.00_r8, 21.00_r8, 27.00_r8, 27.00_r8/

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
integer :: year = 2003, month =1, day =1, tot_days = 31
integer :: max_num = 800000, select_obs = 0
character(len = 129) :: ObsBase = 'temp_obs.'
logical :: ADPUPA = .false., AIRCAR = .false., AIRCFT = .false., SATWND = .false., &
           SATEMP = .false., SFCSHP = .false., ADPSFC = .false.

logical :: obs_U  = .false., obs_V  = .false., obs_T  = .false. , &
           obs_PS = .false., obs_QV = .false., daily_file = .true.

real(r8) :: lon1 = 0.0_r8,    &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound 
            lat1 = -89.0_r8,  &   !  lower latitude bound
            lat2 = 89.0_r8        !  upper latitude bound

namelist /ncepobs_nml/year, month, day, tot_days, max_num, select_obs, ObsBase, &
        ADPUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND, &
        obs_U, obs_V, obs_T, obs_PS, obs_QV, daily_file, lon1, lon2, lat1, lat2

! ----------------------------------------------------------------------
! Select observation types using NCEP categories (when select_obs \= 0).
!
!  ADPUPA: upper-air reports (mostly radiosonde plus few dropsonde, PIBAL)
!  AIRCFT: Conv. (AIREP, PIREP) and ASDAR aircraft reports
!  AIRCAR: ACARS sircraft reports
!  SATEMP: ATOVS retrived temperature
!  SFCSHP: SURFACE MARINE reports
!  ADPSFC: SURFACE LAND SYNOPTIC STATION reports
!  SATWND: Satellite derived wind reports
! ----------------------------------------------------------------------

!  Select variables of U, V, T, QV, PS using the logicals:
!  obs_U   obs_V   obs_PS   obs_T   obs_QV  

call initialize_utilities('create_real_obs')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()

! Read the namelist entry
call find_namelist_in_file("input.nml", "ncepobs_nml", iunit)
read(iunit, nml = ncepobs_nml, iostat = io)
call check_namelist_read(iunit, io, "ncepobs_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'create_real_obs','ncepobs_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=ncepobs_nml)
if (do_output()) write(     *     , nml=ncepobs_nml)

lon1 = min(max(lon1,0.0_r8),360.0_r8)
lon2 = min(max(lon2,0.0_r8),360.0_r8)
if ( lon1 > lon2 ) lon2 = lon2 + 360.0_r8

! Loop through the days interested.

do ii = 1, tot_days
 
  day1 = day + ii - 1

  ! define observation filename
  write(obsdate, '(i4.4,i2.2,i2.2)') year, month, day1

  ! set the obs sequence of the day (daily or 6 hourly)
  if(daily_file) then
    kbeg = 5
    kend = 5
    output_name = 'obs_seq'//obsdate
  else
    kbeg = 1
    kend = 4
  endif

  do kkk = kbeg, kend

    if (daily_file) then
      hour1 = ''
    else
      hour1 = obstime(kkk)
    end if

    seq = real_obs_sequence(year, month, day1, hour1, max_num, select_obs, &
         ObsBase, ADPUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND, &
         obs_U, obs_V, obs_T, obs_PS, obs_QV, bin_beg(kkk), bin_end(kkk), & 
         lon1, lon2, lat1, lat2)

    ! output the daily sequence to a file
    if(.not. daily_file) output_name = 'obs_seq'//obsdate//obstime(kkk)
    call write_obs_seq(seq, output_name)

    ! release the memory of the seq.
    call destroy_obs_sequence(seq)

  enddo

enddo

call timestamp(source,revision,revdate,'end') ! close the log file.

end program create_real_obs
