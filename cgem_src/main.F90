program main

implicit none

integer(kind=8) :: TC_8  ! Current time in seconds since Model_dim::iYr0.
integer  :: istep     ! Current time step
integer :: dT !Timestep
integer :: id !element, loop over nea

#ifdef DEBUG
write(6,*) "In main, before sgrid"
#endif

call sgrid(TC_8,dT)

!!======================================================================     
!    Subroutine CGEM( TC_8, istep )

istep = 1
#ifdef DEBUG
write(6,*) "In main, before cgem, TC_8,istep:",TC_8,istep
#endif

!id will be from loop over nea
id = 1
call cgem(TC_8,istep,id)
TC_8 = TC_8 + dT

#ifdef DEBUG
write(6,*) "In main, after cgem, TC_8,istep:",TC_8,istep
#endif


end program main
