module params

use kinds, only : i_kind, r_single, r_kind

implicit none

! All variables and subroutines are public

! namelist variables, required by GSI/enkf
!  nsats_rad: the total number of satellite data types to read.
!  sattypes_rad:  strings describing the satellite data type (which form part
!   of the diag* filename).
!  sattypes_oz :  strings describing the ozone satellite data type (which form
!   part of the diag* filename).
integer(i_kind), parameter :: nsatmax_rad = 200
integer(i_kind), parameter :: nsatmax_oz  = 100
integer(i_kind), parameter :: ntypes_to_compute_FO_max = 200
integer(i_kind), parameter :: num_obs_platforms_max = 200
integer(i_kind)    :: nsats_rad,nsats_oz, num_obs_platforms

character(len=20), dimension(nsatmax_rad) :: sattypes_rad = ' ' 
character(len=20), dimension(nsatmax_oz)  ::sattypes_oz   = ' ' 
character(len=10)  :: datestring                        = '0000000000' ! if 0000000000 will not be used.
character(len=500) :: datapath                          = ' ' ! mandatory, path to data directory (include trailing slash)
character(len = 129) :: write_FO_for_these_obs_types(ntypes_to_compute_FO_max) = ' ' 
character(len = 129) :: exclude_these_obs_types(ntypes_to_compute_FO_max)      = ' '
! observation platforms/types to read from netcdf files
character(len=129), dimension(num_obs_platforms_max):: obs_platforms = ' ' !!  (/'aircraft', 'gnssro', 'satwind', 'sfc', 'sondes' /)

! namelist variables for DART application
character(len=129) :: obs_seq_out_filename      = 'obs_seq.out'
logical            :: convert_conv              = .true. 
logical            :: convert_sat               = .false.
logical            :: write_prior_copies        = .false.
integer(i_kind)    :: ens_size                  = 1
integer(i_kind)    :: output_option             = 1
logical            :: lie_about_ob_times        = .false.
logical            :: recenter_about_mean_prior = .true.
logical            :: output_bias_preds         = .true.
integer(i_kind)    :: num_ufo_procs = 1
integer(i_kind)    :: npred         = 7        ! number of radiance biases predictors
logical            :: adp_anglebc   = .true.   ! logical to turn off or on the variational radiance angle bias correction
integer(i_kind)    :: angord        = 4        ! order for angle bias correction
logical            :: emiss_bc      = .true.   ! logical to turn off or on the emissivity predictor
logical            :: process_sondes_Tv = .false.   ! process Virtual Temperature radiosonde obs?

! observation arrays filled by EnKF routines
real(r_single),    allocatable, dimension(:)   :: obsprd_prior
real(r_single),    allocatable, dimension(:)   :: ensmean_obnobc
real(r_single),    allocatable, dimension(:)   :: ensmean_ob
real(r_single),    allocatable, dimension(:)   :: ob
real(r_single),    allocatable, dimension(:)   :: oberrvar
real(r_single),    allocatable, dimension(:)   :: oberrvar_orig
real(r_single),    allocatable, dimension(:)   :: obloclon
real(r_single),    allocatable, dimension(:)   :: obloclat
real(r_single),    allocatable, dimension(:)   :: obpress
real(r_single),    allocatable, dimension(:)   :: obtime
real(r_single),    allocatable, dimension(:,:) :: biaspreds
real(r_single),    allocatable, dimension(:,:) :: anal_ob
integer(i_kind),   allocatable, dimension(:)   :: stattype
integer(i_kind),   allocatable, dimension(:)   :: indxsat
character(len=20), allocatable, dimension(:)   :: obtype
real(r_single),    allocatable, dimension(:,:) :: anal_ob_chunk ! For mpi_scatterv implementation of parallel writes

contains

subroutine obsmod_cleanup()
! deallocate arrays...these were allocated in mpi_getobs
   if (allocated(obsprd_prior))   deallocate(obsprd_prior)
   if (allocated(ensmean_ob))     deallocate(ensmean_ob)
   if (allocated(ensmean_obnobc)) deallocate(ensmean_obnobc)
   if (allocated(ob))             deallocate(ob)
   if (allocated(oberrvar))       deallocate(oberrvar)
   if (allocated(oberrvar_orig))  deallocate(oberrvar_orig)
   if (allocated(obloclon))       deallocate(obloclon)
   if (allocated(obloclat))       deallocate(obloclat)
   if (allocated(obpress))        deallocate(obpress)
   if (allocated(obtime))         deallocate(obtime)
   if (allocated(biaspreds))      deallocate(biaspreds)
   if (allocated(anal_ob))        deallocate(anal_ob)
   if (allocated(stattype))       deallocate(stattype)
   if (allocated(indxsat))        deallocate(indxsat)
   if (allocated(obtype))         deallocate(obtype)
end subroutine obsmod_cleanup

end module params
