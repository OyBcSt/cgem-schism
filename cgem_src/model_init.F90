subroutine model_init(TC_8,dT_dum,nsteps)

use date_time
use grid_vars
use cgem_vars
use OUTPUT_NETCDF_CGEM
use OUTPUT

implicit none

integer(8), intent (out) :: TC_8  !Integer seconds
integer, intent (out) :: dT_dum, nsteps  !pass back dT
character(100) :: BASE_NETCDF_OUTPUT_FILE_NAME
!Allocates nea, km(nvrt), nospA/Z, ff

call grid_vars_allocate
call grid_init
call cgem_vars_allocate
call cgem_init
! Compute starting time of run in seconds since Model_dim::iYr0:
TC_8 = TOTAL_SECONDS( iYr0, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )
nsteps = nstep !nstep is from grid_vars
call rad_init(TC_8)
dT_dum = dT

call OUTPUT_NETCDF_CGEM_allocate
BASE_NETCDF_OUTPUT_FILE_NAME = "cgemForSchism"
call Init_Output_CGEM(BASE_NETCDF_OUTPUT_FILE_NAME)

!nvrt = km 
!idry_e = 0
!!idry_e(:) - wet/dry flag
!allocate(idry_e(nea))


end subroutine model_init 
