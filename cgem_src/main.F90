program main

use grid, only:dT,nstep
use cgem_misc

implicit none

integer         :: TC_8  ! Current time in seconds since Model_dim::iYrS.
integer         :: istep ! Current time step
integer         :: ivar,ivark  !Which variable to print
integer         :: k,skip

!Initialize grid
call grid_setup

!Initialize cgem
call cgem_setup(22)

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

ff(:,:) = ff(:,:) + ff_new(:,:)*dTd

do k=1,km

  ff(k,iA(:)) = AMAX1(ff(k,iA(:)),1.)
  if(which_quota.eq.1) then
    ff(k,iQn(:)) = AMAX1(ff(k,iQn(:)),QminN(:))
  else
    ff(k,iQn(:)) = AMIN1(AMAX1(ff(k,iQn(:)),QminN(:)),QmaxN(:))
  endif
  if(which_quota.eq.1) then
    ff(k,iQp(:)) = AMAX1(ff(k,iQp(:)),QminP(:))
  else
    ff(k,iQp(:)) = AMIN1(AMAX1(ff(k,iQp(:)),QminP(:)),QmaxP(:))
  endif
  ff(k,iZ(:)) = AMAX1(ff(k,iZ(:)),1.)

  skip = 3*nospA+2*nospZ+1

  ff(k,skip:nf) = AMAX1(ff(k,skip:nf),0.)


enddo

TC_8 = TC_8 + dT
call grid_update(TC_8)
enddo

#ifdef DEBUG
write(6,*) "In main, after cgem_step, TC_8,istep:",TC_8,istep
#endif

end program main
