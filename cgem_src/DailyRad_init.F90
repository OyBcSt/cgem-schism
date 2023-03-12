subroutine DailyRad_init(TC_8)

   use grid, only: km,iYrS,d,dz,d_sfc,lat,lon,dT
   use cgem
   use cgem_light
   use date_time

   implicit none

!inputs
    integer, intent(in) ::  TC_8  ! Model time (seconds from iYrS)

!loops
    integer k,isp

! Variables needed for light routine and calc_Agrow
    real       :: aDailyRad_k(km)
    real, dimension(nospA,km) :: A_k    ! Phytoplankton number density (cells/m3)
    real, dimension(km) :: OM1A_k, OM1Z_k, OM1SPM_k, OM1BC_k !POC in g/m3
    real, dimension(km) :: CDOM_k    ! CDOM, ppb
    real, parameter :: C_cf  = 12.0E-3    ! C conversion factor (mmol-C/m3 to g-C/m3) 

  integer :: jul_day
  integer :: iYr, iMon, iDay, iHr, iMin, iSec
  real :: rhr ! decimal hour with fraction thereof
  real :: bottom_depth(km)
  real :: Zenith,  calc_solar_zenith !Solar Zenith Angle
  real :: SfcRad ! solar irradiance just below the surface
  ! The avg. clear-sky solar energy is 1200 W/m2. Convert it to photons/cm2/sec
  real, parameter :: solconst = 1200.00 * 2.77e14  ! Morel and Smith (1974) 
  real, parameter :: RADCONV = 1./6.0221413*1.e-19 ! Convert quanta/cm2/s to mol quanta/m2/s
                                             ! mol/m2/s = quanta/cm2/s * 1
                                             ! mol/Avogadro# * 10,000cm2/m2
                                             !          =  (1/6.022e23) * 1.0e4
                                             !          = (1/6.022)e-23 * 1.0e4
                                             !          = (1/6.0221413)e-19
  real :: Chla_tot_k(km)! Total Chl-a concentration (mg/m3) 
  real :: aRadMid(km) ! Holds desired output from Call_IOP_Par
  real :: aRadbot
!Conversions
!-- Integer Number of sec in 24 hr day
      integer, parameter :: iSDay = 86400


       ! Initialize previous day's irradiance for Chl:C calculation
                do k = 1, km 
                   do isp = 1, nospA
                      A_k(isp,k) = ff(k,isp) ! Phytoplankton in group isp, cells/m3
                   enddo
                   CDOM_k(k)  = ff(k,iCDOM)  ! CDOM is in ppb
                                 ! Convert mmol/m3 to g carbon/m3; CF_SPM is
                                 ! river specific
                                 ! and converts river OM to riverine SPM
                   OM1SPM_k(k) = (ff(k,iOM1_R) * C_cf) / CF_SPM
                   OM1Z_k(k)   = ff(k,iOM1_Z)  * C_cf   ! Convert mmol/m3 to g carbon/m3
                   OM1A_k(k)   = ff(k,iOM1_A)  * C_cf   ! Convert mmol/m3 to g carbon/m3
                   OM1BC_k(k)  = ff(k,iOM1_BC) * C_cf   ! Convert mmol/m3 to g carbon/m3
                enddo
!--------------
!call DailyRad_init(TC_8, lat, lon, d, dz, d_sfc, A_k, &
!                     & CDOM_k, OM1A_k, OM1Z_k, OM1SPM_k, OM1BC_k, aDailyRad_k)

  ! Use fixed C:Chla to estimate chlorophyll a concentration
  do k = 1, km
     Chla_tot_k(k) = 0.0
     do isp = 1, nospA
        Chla_tot_k(k) =  Chla_tot_k(k)  + A_k(isp,k) * Qc(isp) * 12. * (1./CChla(isp))
     enddo
     ! Init for later
     aRadSum(k) = 0.0
     bottom_depth(k) = d(k) ! Depth from Surface to bottom of cell
  enddo


  ! Convert time counter into date so we can get the
  ! day to calculate the correct solar path.

  CALL DATE_TIMESTAMP( iYrS, TC_8, iYr, iMon, iDay, iHr, iMin, iSec )

! Now calculate the Julian Day associated with model time TC_8
      jul_day = JDAY_IN_YEAR(iYr, iMon, iDay)

  do iSec = 1, iSDay, dT
     rhr = REAL(iSec) / REAL(3600)
     Zenith = calc_solar_zenith(lat,lon,rhr,jul_day)
     SfcRad = solconst * AMAX1( COS(Zenith), 0.0)    ! COS(Z)<= 0 means night

     if(SfcRad .gt. 0.) then
        Call Call_IOP_PAR(TC_8, SfcRad, CDOM_k, Chla_tot_k, &
        & OM1A_k, OM1Z_k, OM1SPM_k, OM1BC_k, d(km), dz, km, d_sfc, aRadMid,aRadbot)

        ! Add to running total
        aRadSum(:) = aRadSum(:) + aRadMid(:)
     endif
  enddo


  ! Copy result to return variable
  ! Need to convert from quanta/cm2/s to average mol quanta/m2/d
  aDailyRad_k(:) = aRadSum(:) * RADCONV * dT

!-------------
  aDailyRad = aDailyRad_k



!------------

return

END Subroutine DailyRad_init

