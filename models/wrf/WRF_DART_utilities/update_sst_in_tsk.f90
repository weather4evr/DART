! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program update_sst_in_tsk
! This routine is to update TSK over the ocean with SST data from an external source.
! So-Young Ha (MMM/NCAR) 3-15-2011

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/wrf/WRF_DART_utilities/update_sst_in_tsk.f90 $
! $Id: update_sst_in_tsk.f90 4620 2011-01-04 17:07:38Z thoar $
! $Revision: 4620 $
! $Date: 2011-01-04 10:07:38 -0700 (Tue, 04 Jan 2011) $

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_MSG, E_ERR,       &
                          open_file, close_file, nc_check, get_next_filename, &
                          find_namelist_in_file, check_namelist_read,         &
                          do_nml_file, do_nml_term, nmlfileunit,              &
                          initialize_utilities, finalize_utilities
!use        model_mod, only : get_wrf_date   ! CSS commented out. Being dependent on model_mod means being dependent on 
                                             !  many other modules that have changed, making this code harder to compile.Took subrtou
                                             !  Simple took get_wrf_date from model_mod and pasted it below
use time_manager_mod, only : time_type, read_time, set_date, print_time
use netcdf

implicit none

! version controlled file description for error handling, do not edit
! split into separate lines; getting too long for the absoft compiler
character(len=128), parameter :: source = &
 "$URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/wrf/WRF_DART_utilities/update_sst_in_tsk.f90 $"
character(len=128), parameter :: revision = "$Revision: 4620 $"
character(len=128), parameter :: revdate  = "$Date: 2011-01-04 10:07:38 -0700 (Tue, 04 Jan 2011) $"

! variables used to read the netcdf info
integer, parameter :: maxd = 7
integer :: i, j, ndims, odims, ncrc, mdims, k
integer :: ncinid, ncoutid        ! netcdf id for file
integer :: invarid, outvarid, maskvarid
integer :: tinid, toutid
integer :: dimid(maxd), dimlen(maxd), odimid(maxd), odimlen(maxd), mdimid(maxd)
integer :: tdimid(2)
character(128) :: dimname(maxd), odimname(maxd)
integer ::  ninDimensions,  ninVariables,  ninAttributes,  inunlimitedDimID
integer :: noutDimensions, noutVariables, noutAttributes, outunlimitedDimID

! arrays for all possible dimensions
!real(r8), pointer ::   oned(:)       
!real(r8), pointer ::   twod(:,:)
!real(r8), pointer :: threed(:,:,:)
!real(r8), pointer ::  fourd(:,:,:,:)
!real(r8), pointer ::  fived(:,:,:,:,:)
!real(r8), pointer ::   sixd(:,:,:,:,:,:)
!real(r8), pointer :: sevend(:,:,:,:,:,:,:)

real(r8), allocatable ::   xin(:,:,:)
real(r8), allocatable ::  xout(:,:,:)
real(r8), allocatable :: xmask(:,:,:)

logical, save :: module_initialized = .false.

character(len=80)  :: varname
character(len=19)  :: timestring1, timestring2
character(len=128) :: msgstring, msgstring2, tmpstring
integer :: iunit, io
integer :: ivtype, tdims,idims(2)
integer :: in_yr, in_mo, in_dy, in_hr, in_mn, in_sc
integer :: out_yr, out_mo, out_dy, out_hr, out_mn, out_sc
type(time_type)   :: dart_time_in, dart_time_out

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
character(len=128) ::  update_in_file_name = 'wrfinput_d01_old'
character(len=128) :: update_out_file_name = 'wrfinput_d01_new'
character(len=128) ::  input_field = 'SST'
character(len=128) :: output_field = 'TSK'
character(len=128) ::   mask_field = 'LANDMASK'
logical :: date_check = .true.		! Check if the dates for input_file (SST) and output_file (TSK) are same.
logical :: debug = .false.		! or .true.

namelist /update_sst_in_tsk_nml/  update_in_file_name, update_out_file_name, &
                                  input_field, output_field, mask_field, date_check, debug

! main code here
 
! flow:
!   initialization
call initialize_utilities('update_sst_in_tsk')
call initialize_module()

! Read the namelist entry
call find_namelist_in_file("input.nml", "update_sst_in_tsk_nml", iunit)
read(iunit, nml = update_sst_in_tsk_nml, iostat = io)
call check_namelist_read(iunit, io, "update_sst_in_tsk_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=update_sst_in_tsk_nml)
if (do_nml_term()) write(     *     , nml=update_sst_in_tsk_nml)

if (debug) then
   call error_handler(E_MSG, 'update_sst_in_tsk', ' debug on')
endif

! whether to fail or just warn if a field is not found

call error_handler(E_MSG, 'update_sst_in_tsk', ' reading file: '//trim(update_in_file_name))
call error_handler(E_MSG, 'update_sst_in_tsk', ' overwriting file: '//trim(update_out_file_name))

!   do they exist?  can they be opened?
!   update_in_file_name & update_out_file_name are netcdf
!   fieldlist is ascii, one wrf fieldname per line

! open the files
call nc_check(nf90_open( update_in_file_name, NF90_NOWRITE,    ncinid), 'nf90_open',  'update_in_file_name')
call nc_check(nf90_open(update_out_file_name, NF90_WRITE,     ncoutid), 'nf90_open', 'update_out_file_name')

if (debug) then
   call nc_check(nf90_inquire( ncinid,  ninDimensions,  ninVariables, &
                  ninAttributes,  inunlimitedDimID), 'nf90_inquire',  'update_in_file_name')
   call nc_check(nf90_inquire(ncoutid, noutDimensions, noutVariables, &
                 noutAttributes, outunlimitedDimID), 'nf90_inquire', 'update_out_file_name')

   write(msgstring, *) ' update_in_file_name ndim, nvar, nattr:', ninDimensions, &
                       ninVariables, ninAttributes
   call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
   write(msgstring, *) 'update_out_file_name ndim, nvar, nattr:', noutDimensions, &
                       noutVariables, noutAttributes
   call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
endif

! Check dates between input_file and output_file
if ( date_check ) then
   call nc_check( nf90_inq_varid( ncinid, "Times",  tinid), 'update_sst_in_tsk', &
                 'inq_varid Times in input_file' )
   call nc_check( nf90_inq_varid(ncoutid, "Times", toutid), 'update_sst_in_tsk', &
                 'inq_varid Times in output_file' )
   call nc_check( nf90_inquire_variable( ncinid,  tinid, varname, xtype=ivtype, &
                  ndims=tdims, dimids=tdimid), 'update_sst_in_tsk', &
                 'inquire_variable Times in input_file' )
   do i=1,tdims
      call nc_check( nf90_inquire_dimension(ncinid, tdimid(i), &
                     len=idims(i)),'update_sst_in_tsk','inquire_dimensions Times' )
      if(debug) write(*,*) ' dimension ',i,idims(i)
   enddo

   call nc_check( nf90_get_var( ncinid,  tinid, timestring1, &
                  start = (/ 1, idims(2) /)), 'update_sst_in_tsk', &
                 'get_var Times in input_file' ) 
   call nc_check( nf90_get_var( ncinid, toutid, timestring2, &
                  start = (/ 1, idims(2) /)), 'update_sst_in_tsk', &
                 'get_var Times in output_file' ) 

   call get_wrf_date(timestring1, in_yr, in_mo, in_dy, in_hr, in_mn, in_sc)
   !dart_time_in  = set_date(in_yr, in_mo, in_dy, in_hr, in_mn, in_sc)
   call get_wrf_date(timestring2, out_yr, out_mo, out_dy, out_hr, out_mn, out_sc)
   !dart_time_out = set_date(out_yr, out_mo, out_dy, out_hr, out_mn, out_sc)
   if((in_yr.ne.out_yr).or.(in_mo.ne.out_mo).or.(in_dy.ne.out_dy)) then
       write(msgstring,*)   'input_time:',in_yr,in_mo,in_dy,&
                          ' output_time:',out_yr,out_mo,out_dy
       call error_handler(E_ERR,'update_sst_in_tsk: Time is different.',msgstring)
   endif
endif

!    input file to get data from
!    list of netcdf fields to copy over
!    output file to be updated in place

   ! inquire in input for fieldname
   ! inquire in output
   ncrc = nf90_inq_varid( ncinid, trim(input_field),  invarid)
   if (ncrc /= NF90_NOERR) then
      tmpstring = ' not found in input file '//trim(update_in_file_name)
      msgstring = 'variable '//trim(input_field)//trim(tmpstring)
      call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, source, revision, revdate) 
   endif
   ncrc = nf90_inq_varid(ncoutid, trim(mask_field),maskvarid)
   if (ncrc /= NF90_NOERR) then
      tmpstring = ' not found in output file '//trim(update_out_file_name)
      msgstring = 'variable '//trim(mask_field)//trim(tmpstring)
      call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, source, revision, revdate) 
   endif
   ncrc = nf90_inq_varid(ncoutid, trim(output_field), outvarid)
   if (ncrc /= NF90_NOERR) then
      tmpstring = ' exists in output file '//trim(update_out_file_name)
      msgstring = 'variable '//trim(output_field)//trim(tmpstring)
      msgstring2 = 'but was not found in output file '//trim(update_out_file_name)
      call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, &
                         source, revision, revdate, text2=msgstring2) 
   endif

   if (debug) then
      write(msgstring, *) ' invarid: ', trim(input_field)//' ',  invarid
      call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
      write(msgstring, *) 'maskvarid: ', trim(mask_field)//' ', maskvarid
      call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
      write(msgstring, *) 'outvarid: ', trim(output_field)//' ', outvarid
      call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
   endif

   ! get dimensions and make sure they match

   call nc_check(nf90_inquire_variable( ncinid,  invarid, ndims=ndims,  dimids=dimid), &
                 'nf90_inquire_variable',  'update_in_file_name/'//trim(input_field))
   call nc_check(nf90_inquire_variable(ncoutid,maskvarid, ndims=mdims, dimids=mdimid), &
                 'nf90_inquire_variable', 'update_out_file_name/'//trim(mask_field))
   call nc_check(nf90_inquire_variable(ncoutid, outvarid, ndims=odims, dimids=odimid), &
                 'nf90_inquire_variable', 'update_out_file_name/'//trim(output_field))

   if (ndims /= odims) then
      write(msgstring, *) 'variable ', trim(input_field), ' and ', trim(output_field), &
         ' have different numbers of dimensions in the two files'
      call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
      write(msgstring, *) 'input dimension size ', ndims, ' does not match output ', odims
      call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, source, revision, revdate)
   endif
   if (mdims /= odims) then
      write(msgstring, *) 'variable ', trim(mask_field), ' and ', trim(output_field), &
         ' have different numbers of dimensions in the two files'
      call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
      write(msgstring, *) 'input dimension size ', mdims, ' does not match output ', odims
      call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, source, revision, revdate)
   endif
   
   do j=1,ndims
      call nc_check( nf90_inquire_dimension( ncinid,  dimid(j),  dimname(j),  dimlen(j)), &
                   'nf90_inquire_dimension',  'update_in_file_name/'//trim( dimname(j)) )
      if (debug) then
         write(msgstring, '(2A,I5,A,I8,2A)') trim(update_in_file_name), ' dim: ', j, ' len: ', dimlen(j), ' name: ', trim(dimname(j))
         call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
      endif
      call nc_check( nf90_inquire_dimension(ncoutid, odimid(j), odimname(j), odimlen(j)), &
                   'nf90_inquire_dimension', 'update_out_file_name/'//trim(odimname(j)) )
      if (debug) then
         write(msgstring, '(2A,I5,A,I8,2A)') trim(update_out_file_name), ' dim: ', j, ' len: ', odimlen(j), ' name: ', trim(odimname(j))
         call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
      endif
      
      if (dimlen(j) /= odimlen(j)) then
         write(msgstring, *) 'variable ', trim(input_field), ' and ', trim(output_field), &
                             ' have different dimensions in the two files'
         call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
         write(msgstring, *) 'input dim length ', dimlen(j), ' does not match output ', odimlen(j)
         call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, source, revision, revdate)
      endif

      ! only possible if the unlimited dim is declared but hasn't been written to
      if (dimlen(j) == 0) then
         write(msgstring, *) trim(input_field), 'will be skipped because it is empty in input file'
         call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
      endif

   enddo


   select case(ndims)
      case (1)
         write(tmpstring, '(2A,1I5,2A)') trim(input_field), '(', dimlen(1),   ') to ', trim(output_field)
      case (2)
         write(tmpstring, '(2A,2I5,2A)') trim(input_field), '(', dimlen(1:2), ') to ', trim(output_field)
      case (3)
         write(tmpstring, '(2A,3I5,2A)') trim(input_field), '(', dimlen(1:3), ') to ', trim(output_field)
      case (4)
         write(tmpstring, '(2A,4I5,2A)') trim(input_field), '(', dimlen(1:4), ') to ', trim(output_field)
      case (5)
         write(tmpstring, '(2A,5I5,2A)') trim(input_field), '(', dimlen(1:5), ') to ', trim(output_field)
      case (6)
         write(tmpstring, '(2A,6I5,2A)') trim(input_field), '(', dimlen(1:6), ') to ', trim(output_field)
      case (7)
         write(tmpstring, '(2A,7I5,2A)') trim(input_field), '(', dimlen(1:7), ') to ', trim(output_field)
      case default
         ! "can't happen"
         write(msgstring, *) 'array dimension is illegal value: ', ndims
         call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, source, revision, revdate)
   end select

   ! announce what we're about to do
   write(msgstring, *) 'copying ', trim(tmpstring)
   call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)

   ! allocate right dim array
   ! read/write and then deallocate
   ! FIXME: we want to allow only the 3-dimensional case for the time being. (SYHA)

   select case(ndims)
      case (3)	
         allocate(  xin(dimlen(1),dimlen(2),dimlen(3)))
         allocate( xout(dimlen(1),dimlen(2),dimlen(3)))
         allocate(xmask(dimlen(1),dimlen(2),dimlen(3)))
         call nc_check(nf90_get_var( ncinid,  invarid,  xin), 'nf90_get_var',  'update_in_file_name')
         call nc_check(nf90_get_var(ncoutid,maskvarid,xmask), 'nf90_get_var', 'update_out_file_name')
         call nc_check(nf90_get_var(ncoutid, outvarid, xout), 'nf90_get_var', 'update_out_file_name')

         !call fillup_tsk( dimlen(1), dimlen(2), dimlen(3), xin, xmask, xout )
         do k = 1, dimlen(3)
          do j = 1, dimlen(2)
           do i = 1, dimlen(1)

              if( ( xmask(i,j,k) < 0.5_r8 ) .and. &
                  ( xin(i,j,k) > 170.0_r8 ) .and. ( xin(i,j,k) < 400.0_r8 ) ) then
                    xout(i,j,k) = xin(i,j,k)
              end if

           end do
          end do
         end do

         call nc_check(nf90_put_var(ncoutid, outvarid, xout), 'nf90_put_var', 'update_out_file_name')
         deallocate(xin)
         deallocate(xout)
         deallocate(xmask)
      case default
         ! "really can't happen"
         write(msgstring, *) 'array dimension is illegal value: ', ndims
         call error_handler(E_ERR, 'update_sst_in_tsk', msgstring, source, revision, revdate)
   end select

!  close up
call nc_check(nf90_close( ncinid), 'nf90_close',  'update_in_file_name')
call nc_check(nf90_close(ncoutid), 'nf90_close', 'update_out_file_name')

if (debug) then
   write(msgstring, *) 'closing files',  trim(update_in_file_name), ' and ', trim(update_out_file_name)
   call error_handler(E_MSG, 'update_sst_in_tsk', msgstring)
endif

call finalize_utilities('update_sst_in_tsk')

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call register_module(source, revision, revdate)
  module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine fillup_tsk( nx, ny, nz, xin, xmask, xout )
    
    integer :: i, j, k
    integer :: nx, ny, nz
    real,   intent(in) ::    xin(nx,ny,nz)               
    real,   intent(in) ::  xmask(nx,ny,nz)               
    real, intent(inout) ::  xout(nx,ny,nz)               
    
    do k = 1, nz
     do j = 1, ny
      do i = 1, nx

         if( ( xmask(i,j,k) < 0.5_r8 ) .and. &
             ( xin(i,j,k) > 170.0_r8 ) .and. ( xin(i,j,k) < 400.0_r8 ) ) then
               xout(i,j,k) = xin(i,j,k)
         end if    

      end do
     end do
    end do

end subroutine fillup_tsk

!----------------------------------------------------------------------
! Returns integers taken from tstring
! It is assumed that the tstring char array is as YYYY-MM-DD_hh:mm:ss

subroutine get_wrf_date (tstring, year, month, day, hour, minute, second)

integer,           intent(out) :: year, month, day, hour, minute, second
character(len=19), intent(in)  :: tstring

read(tstring( 1: 4),'(i4)') year
read(tstring( 6: 7),'(i2)') month
read(tstring( 9:10),'(i2)') day
read(tstring(12:13),'(i2)') hour
read(tstring(15:16),'(i2)') minute
read(tstring(18:19),'(i2)') second

return

end subroutine get_wrf_date


!----------------------------------------------------------------------
end program
