module grid_vars

!CGEM STATE VARIABLES
use date_time
use, intrinsic :: iso_fortran_env, only: stderr => error_unit
!use schism_glbl, only: rkind

implicit none

SAVE

!Grid parameters
integer, parameter ::  nea = 1
integer, parameter :: km = 5
integer, parameter :: iYr0 = 2002
integer(8) :: START_SECONDS,END_SECONDS

!--Run Specifics---------------
integer iYrS,iMonS,iDayS,iHrS,iMinS,iSecS
integer iYrE,iMonE,iDayE,iHrE,iMinE,iSecE
integer dT,nstep

!Constants
real, parameter :: SDay = 86400.0  ! # of sec in 24 hr day
integer :: StepsPerDay                   ! Time steps per day
real :: dTd                           ! Timestep in days 
integer :: dT_out !Output timestep
integer iout
!-output
integer istep_out

!This should all be from SCHISM
!Solar radiation
real :: Rad
!Lat/lon of elements
real :: lat,lon
!Depth stuff
real, allocatable :: d(:) ! depth cell bottom to surface
real, allocatable :: d_sfc(:) ! depth cell center to surface 
real, allocatable :: dz(:)   !Thickness of cell
!'tracers'
real, allocatable :: S(:),T(:)

contains

subroutine grid_vars_allocate()

integer i
integer ierr
integer :: counter = 0


#ifdef DEBUG
write(6,*) "Begin cgem_vars_allocate" 
#endif

!Depth stuff
allocate(d(km),stat=ierr) ! depth cell bottom to surface  
if(ierr.ne.0) write(6,*) "error in allocating:d"
allocate(d_sfc(km),stat=ierr) ! depth cell center to surface  
if(ierr.ne.0) write(6,*) "error in allocating:d_sfc"
allocate(dz(km),stat=ierr) !Thickness of cell
if(ierr.ne.0) write(6,*) "error in allocating:dz"


!SCHISM 'tracers'
allocate(S(km),stat=ierr) ! Salinity (psu) 
if(ierr.ne.0) write(6,*) "error in allocating:S"
allocate(T(km),stat=ierr) ! Temperature (C) 
if(ierr.ne.0) write(6,*) "error in allocating:T"

#ifdef DEBUG
write(6,*) "End grid_vars_allocate"
#endif


return
end subroutine grid_vars_allocate

subroutine grid_init 

integer                          :: istat,iunit
integer :: k
real S_init,T_init,depth_in,Rad_in
real lat_in,lon_in
character(len=1000) :: line
!http://degenerateconic.com/namelist-error-checking.html
namelist /grid/ iYrS,iMonS,iDayS,iHrS,iMinS,iSecS,iYrE,iMonE,iDayE,iHrE,iMinE,iSecE,&
 dT,dT_out,lon_in,lat_in,depth_in,Rad_in,S_init,T_init 

!Later, look at this python namelist thing:
!https://github.com/marshallward/f90nml
!install
!conda create --prefix ./env_f90 -c conda-forge f90nml
!(conda activate it)
!conda install -c conda-forge ascii_graph

#ifdef DEBUG
write(6,*) "Begin grid_init"
#endif

open(action='read',file='grid.nml',iostat=istat,newunit=iunit)

!namelist /grid/
read(nml=grid,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist grid: '//trim(line)
endif

close(iunit)


!Initialize Time stuff
dTd = dT/SDay         ! Timestep length in units of days
!StepsPerDay = SDay / dT ! Time steps in a day
StepsPerDay = 86400 / dT

! Compute starting time of run in seconds since Model_dim::iYr0:
START_SECONDS = &
TOTAL_SECONDS( iYr0, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )
#ifdef DEBUG
write(6,*) "In ReadInput"
write(6,*) "StartSeconds",START_SECONDS
#endif
END_SECONDS = &
TOTAL_SECONDS( iYr0, iYrE, iMonE, iDayE, iHrE, iMinE, iSecE )
#ifdef DEBUG
write(6,*) "In ReadInput"
write(6,*) "EndSeconds,dT,dT_out",START_SECONDS,dT,dT_out
#endif

nstep = ( END_SECONDS - START_SECONDS ) / dT !number of timesteps in a run
iout = dT_out/dT !output time-interval in timesteps
!nstep_sed = dT_sed / dT   ! number of steps in-between calls to sediment
!diagenesis

!Keep a counter for netCDF timesteps
istep_out = 0

#ifdef DEBUG
write(6,*) "nstep,iout",nstep,iout
#endif

!Set lat/lon
lat = lat_in
lon = lon_in
!Radiation
!Rad = Rad_in
call getSolar( iYr0, START_SECONDS, lon, lat, Rad)

#ifdef DEBUG
write(6,*) "lat,lon,Rad",lat,lon,Rad
#endif


!Define depths
d = depth_in          !Depth of the water column
dz = depth_in/km      !Thickness of a cell (they are the same for now)

!Distance from surface to bottom of cell
d(1) = dz(1)
!Distance from surface to center of cell
d_sfc(1) = dz(1)/2.   !First cell is half of total thickness of first cell
do k=2,km   !Okay, this is stupid and only works for this case...fix later
 d(k) = d(k-1) + dz(k)
 d_sfc(k) = d_sfc(k-1) + dz(k)
enddo

!"Tracers"
S = S_init
T = T_init

#ifdef DEBUG
write(6,*) "S,T",S(1),T(1)
write(6,*) "End grid_init"
#endif

end subroutine grid_init

end module grid_vars
