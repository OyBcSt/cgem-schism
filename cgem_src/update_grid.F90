subroutine update_grid(TC_8,istep)

use date_time
use grid
!use cgem

implicit none

integer(8), intent (in) :: TC_8  !Integer seconds
integer, intent(in) :: istep

!Allocates nea, km(nvrt), nospA/Z, ff
!Update outside vars - rad, T, S, wind
call getSolar( iYrS, TC_8, lon, lat, Rad)

return
end subroutine update_grid 
