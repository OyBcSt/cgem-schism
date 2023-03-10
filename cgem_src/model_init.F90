subroutine model_init(TC_8,dT_dum,nsteps)

use date_time
use grid
use cgem

#ifdef NCFILE
use OUTPUT_NETCDF_CGEM
use OUTPUT
#endif


implicit none

integer(8), intent (out) :: TC_8  !Integer seconds
integer, intent (out) :: dT_dum, nsteps  !pass back dT
character(100) :: BASE_NETCDF_OUTPUT_FILE_NAME

!Initialize grid
call grid_setup

!Initialize cgem
call cgem_setup

! Compute starting time of run in seconds since Model_dim::iYr0:
TC_8 = TOTAL_SECONDS( iYr0, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )
nsteps = nstep !nstep is from grid_vars
call rad_init(TC_8)
dT_dum = dT

#ifdef NCFILE
call OUTPUT_NETCDF_CGEM_allocate
BASE_NETCDF_OUTPUT_FILE_NAME = "cgemForSchism"
call Init_Output_CGEM(BASE_NETCDF_OUTPUT_FILE_NAME)
#endif

return
end subroutine model_init 
