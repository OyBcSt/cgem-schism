program main

implicit none

integer(kind=8) :: TC_8  ! Current time in seconds since Model_dim::iYr0.
integer  :: istep     ! Current time step
integer :: dT !Timestep
integer :: id !element, loop over nea
integer ivar !Which variable to print

#ifdef DEBUG
write(6,*) "In main, before sgrid"
#endif

!Initializes everything
call sgrid(TC_8,dT)

!Need to define vars before figuring out which to print
call Command_Line_Args(ivar)


!!======================================================================     
!    Subroutine CGEM( TC_8, istep )

#ifdef DEBUG
write(6,*) "In main, before cgem, TC_8,istep:",TC_8,istep
#endif

!id will be from loop over nea
id = 1

do istep=1,100
call cgem(TC_8,istep,id,ivar)
TC_8 = TC_8 + dT
enddo

#ifdef DEBUG
write(6,*) "In main, after cgem, TC_8,istep:",TC_8,istep
#endif


end program main
