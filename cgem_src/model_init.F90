subroutine model_init(TC_8,dT_dum,nsteps)

use date_time
use grid
use cgem

implicit none

integer(8), intent (out) :: TC_8  !Integer seconds
integer, intent (out) :: dT_dum, nsteps  !pass back dT

!Initialize grid
call grid_setup

!Initialize cgem
call cgem_setup

! Compute starting time of run in seconds since Model_dim::iYrS:
TC_8 = TOTAL_SECONDS( iYrS, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )
nsteps = nstep !nstep is from grid_vars
call DailyRad_Init(TC_8)
dT_dum = dT

return
end subroutine model_init 
