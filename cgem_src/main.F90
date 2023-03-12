program main

use cgem, only:dT,nstep
use cgem_misc

implicit none

integer         :: TC_8  ! Current time in seconds since Model_dim::iYrS.
integer         :: istep ! Current time step
integer         :: ivar,ivark  !Which variable to print

!Initialize grid
call grid_setup

!Initialize cgem
call cgem_setup

! Compute starting time of run in seconds since Model_dim::iYrS:
TC_8 = TOTAL_SECONDS( iYrS, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )

call DailyRad_Init(TC_8)

call Check_InputFile

!Need to define vars before figuring out which to print
call Command_Line_Args(ivar,ivark)

#ifdef DEBUG
write(6,*) "In main, before cgem_step, TC_8:",TC_8
#endif

do istep=1,nstep
call cgem_step(TC_8,istep,ivar,ivark)
TC_8 = TC_8 + dT
call grid_update(TC_8)
enddo

#ifdef DEBUG
write(6,*) "In main, after cgem_step, TC_8,istep:",TC_8,istep
#endif

end program main
