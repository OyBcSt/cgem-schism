program main

implicit none

integer(kind=8) :: TC_8  ! Current time in seconds since Model_dim::iYrS.
integer  :: istep,nstep     ! Current time step, total timesteps
integer :: dT !Timestep
integer ivar !Which variable to print

#ifdef DEBUG
write(6,*) "In main, before model_init"
#endif


!Initializes everything
call model_init(TC_8,dT,nstep)

call Check_InputFile

!call check_grid(TC_8,dT)

!Need to define vars before figuring out which to print
call Command_Line_Args(ivar)

!!======================================================================     
!    Subroutine CGEM( TC_8, istep )

#ifdef DEBUG
write(6,*) "In main, before run_cgem, TC_8:",TC_8
#endif

do istep=1,nstep
call run_cgem(TC_8,istep,ivar)
TC_8 = TC_8 + dT
call update_grid(TC_8,istep)
enddo

#ifdef DEBUG
write(6,*) "In main, after cgem, TC_8,istep:",TC_8,istep
#endif

#ifdef NCFILE
call Model_Finalize_CGEM
#endif

end program main
