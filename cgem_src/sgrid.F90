subroutine sgrid(TC_8,dT_out)

use date_time
use cgem_vars

implicit none

integer(8), intent (out) :: TC_8  !Integer seconds
integer, intent (out) :: dT_out  !pass back dT

!Allocates nea, km(nvrt), nospA/Z, ff
call cgem_vars_allocate
call cgem_init
! Compute starting time of run in seconds since Model_dim::iYr0:
TC_8 = TOTAL_SECONDS( iYr0, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )

call rad_init(TC_8)

dT_out = dT

!nvrt = km 
!idry_e = 0
!!idry_e(:) - wet/dry flag
!allocate(idry_e(nea))


end subroutine sgrid
