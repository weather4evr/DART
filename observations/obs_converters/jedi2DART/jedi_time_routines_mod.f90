module jedi_time_routines_mod

implicit none

private

! public subroutines
public :: read_jedi_time

contains

subroutine read_jedi_time(datestring,nobs,jedi_time_string,diff_hrs)
   character(len=10), intent(in) :: datestring
   integer, intent(in) :: nobs
   character(len=20), dimension(nobs), intent(in) :: jedi_time_string
   real, dimension(nobs), intent(out) :: diff_hrs

   !   IDATE(1)       YEAR OF CENTURY
   !   IDATE(2)       MONTH OF YEAR
   !   IDATE(3)       DAY OF MONTH
   !   IDATE(4)       HOUR OF DAY
   !   IDATE(5)       MINUTE OF HOUR
   integer :: idate_ref(5), idate(5)
   integer :: i, nmin, nmin_ref

   ! get the reference time (analysis time)
   read(datestring, '(i4,3i2)') idate_ref(1), idate_ref(2), idate_ref(3), idate_ref(4)
   idate_ref(5) = 0

   !2018-04-14T21:05:59Z
   do i = 1,nobs
      read(jedi_time_string(i)(1:4),*)   idate(1)
      read(jedi_time_string(i)(6:7),*)   idate(2)
      read(jedi_time_string(i)(9:10),*)  idate(3)
      read(jedi_time_string(i)(12:13),*) idate(4)
      read(jedi_time_string(i)(15:16),*) idate(5)

      ! get number of minutes from a historical reference date
      call W3FS21(idate_ref, nmin_ref) ! analysis 
      call W3FS21(idate, nmin) ! observaton

      diff_hrs(i) = real(nmin - nmin_ref)/60.
   enddo

end subroutine read_jedi_time

! Subroutine taken from NCEP w3lib
SUBROUTINE W3FS21(IDATE, NMIN)
   integer, intent(in) :: IDATE(5)
   integer, intent(inout) :: NMIN

   integer :: iyear, ijdn, ndays, JDN78
   DATA  JDN78 / 2443510 /

   NMIN  = 0

   IYEAR = IDATE(1)
!
!  COMPUTE JULIAN DAY NUMBER FROM YEAR, MONTH, DAY
!
   IJDN  = IW3JDN(IYEAR,IDATE(2),IDATE(3))
!
!  SUBTRACT JULIAN DAY NUMBER OF JAN 1,1978 TO GET THE
!  NUMBER OF DAYS BETWEEN DATES
!
   NDAYS = IJDN - JDN78

!  NUMBER OF MINUTES
   NMIN = NDAYS * 1440 + IDATE(4) * 60 + IDATE(5)

   RETURN
END SUBROUTINE W3FS21

! Function taken from NCEP w3lib
FUNCTION IW3JDN(IYEAR,MONTH,IDAY)
   integer, intent(in) :: IYEAR, MONTH, IDAY
   integer :: IW3JDN
   IW3JDN  =    IDAY - 32075 &
              + 1461 * (IYEAR + 4800 + (MONTH - 14) / 12) / 4 &
              + 367 * (MONTH - 2 - (MONTH -14) / 12 * 12) / 12 &
              - 3 * ((IYEAR + 4900 + (MONTH - 14) / 12) / 100) / 4 
   RETURN
END FUNCTION

end module jedi_time_routines_mod
