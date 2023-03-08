subroutine update_grid(TC_8,istep)

use date_time
use grid_vars
use cgem_vars

implicit none

integer(8), intent (in) :: TC_8  !Integer seconds
integer, intent(in) :: istep

!Allocates nea, km(nvrt), nospA/Z, ff
!Update outside vars - rad, T, S, wind
call getSolar( iYr0, TC_8, lon, lat, Rad)


#ifdef NCFILE
if ( mod( istep, iout ) .eq. 0 ) then
 istep_out = istep_out + 1
 call Model_Output_CGEM( istep_out )
endif
#endif

end subroutine update_grid 
