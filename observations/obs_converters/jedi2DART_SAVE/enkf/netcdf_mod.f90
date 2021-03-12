!ifort  -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include test.f90 -o a.out
module netcdf_mod

use netcdf
use mpisetup, only: stop2

implicit none

private

! public subroutines
public :: open_netcdf, close_netcdf, get_netcdf_dims, get_netcdf_var_1d, get_netcdf_var_2d

! variables visible to this module nly
integer :: ncstatus

interface get_netcdf_var_1d
   module procedure get_netcdf_var_1d_real
   module procedure get_netcdf_var_1d_integer
   module procedure get_netcdf_var_1d_char
end interface
interface get_netcdf_var_2d
   module procedure get_netcdf_var_2d_real
   module procedure get_netcdf_var_2d_integer
end interface

contains

subroutine open_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer, intent(out) :: ncfileid

   ncstatus = nf90_open(path=trim(adjustl(fname)),mode=nf90_nowrite,ncid=ncfileid)  ! open file
   if ( ncstatus .eq. 0 ) then
!     write(*,fmt='(a)') 'opened '//trim(adjustl(fname))//' for reading'
!     write(*,fmt='(a,i8)') 'fileid = ',ncfileid
   else
      write(*,fmt='(a)') 'error reading '//trim(adjustl(fname))
      call stop2(31) ! stop
   endif

   return
end subroutine open_netcdf

subroutine close_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ncfileid
   ncstatus = nf90_close(ncfileid) ! close file
   if ( ncstatus .ne. 0 ) then
      write(*,fmt='(a)') 'error closing '//trim(adjustl(fname))
      call stop2(32) ! stop
   endif
end subroutine close_netcdf

subroutine get_netcdf_dims(fileid,variable,output)
   integer, intent(in) :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(out) :: output

   integer :: ncdimid, ierr

   ierr = 0
   ncstatus = nf90_inq_dimid(fileid,trim(adjustl(variable)),ncdimid) ; ierr = ierr + ncstatus
   ncstatus = nf90_inquire_dimension(fileid,ncdimid,len=output)      ; ierr = ierr + ncstatus
   if ( ierr /= 0 ) then
      write(0,*) 'Error reading dimension for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(33) ! stop
   endif
!  write(*,fmt='(a,i8)')variable//' = ', output ! print out dimensions for the variable

   return

end subroutine get_netcdf_dims

subroutine get_netcdf_var_1d_real(fileid,variable,dim1,output)

   integer, intent(in) :: fileid, dim1
   character(len=*), intent(in) :: variable
   real, intent(inout), dimension(dim1)  :: output

   integer :: ncvarid, ierr

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(34) ! stop
   endif

   return

end subroutine get_netcdf_var_1d_real

subroutine get_netcdf_var_1d_integer(fileid,variable,dim1,output)

   integer, intent(in) :: fileid, dim1
   character(len=*), intent(in) :: variable
   integer, intent(inout), dimension(dim1)  :: output

   integer :: ncvarid, ierr

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(35) ! stop
   endif

   return

end subroutine get_netcdf_var_1d_integer

subroutine get_netcdf_var_1d_char(fileid,variable,dim1,output)

   integer, intent(in) :: fileid, dim1
   character(len=*), intent(in) :: variable
   character(len=*), intent(inout), dimension(dim1)  :: output

   integer :: ncvarid, ierr

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(36) ! stop
   endif

   return

end subroutine get_netcdf_var_1d_char

subroutine get_netcdf_var_2d_real(fileid,variable,dim1,dim2,output)

   integer, intent(in) :: fileid, dim1, dim2
   character(len=*), intent(in) :: variable
   real, intent(inout), dimension(dim1,dim2)  :: output

   integer :: ncvarid, ierr

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(37) ! stop
   endif

   return

end subroutine get_netcdf_var_2d_real

subroutine get_netcdf_var_2d_integer(fileid,variable,dim1,dim2,output)

   integer, intent(in) :: fileid, dim1, dim2
   character(len=*), intent(in) :: variable
   integer, intent(inout), dimension(dim1,dim2)  :: output

   integer :: ncvarid, ierr

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(38) ! stop
   endif

   return

end subroutine get_netcdf_var_2d_integer

end module netcdf_mod
