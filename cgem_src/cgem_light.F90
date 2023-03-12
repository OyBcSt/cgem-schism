!----------------------------------------------------------------------
      FUNCTION calc_solar_zenith(lat, lon, rhr, julianDay )    &      
      &        RESULT(sunangle) 
      
!-----------------------------------------------------------------------
!     Original MATLAB code Written by: Brad Penta/NRL
!
!     Translated into FORTRAN and 
!     Revised by                     : Barry Herchenroder/EMVL, April 2010
!                                                               June 2011
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Calculate the solar zenith angle (radians) given the latitude 
!  (lat, degrees) longitude (lon, degrees) of a point in the Gulf of 
!  Mexico (GOM) as well as the Julian day of the year (julianDay),
!  the hour in that day (ihr-integer)
!------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! Definitions
  !
  ! sunangle     -> angle of the sun in radians, computed from latitude,
  !                 longitude, Julian day, and time--This is the zenith 
  !
  ! Compute the sun zenith angle in radians for a given Julian day (julianDay),
  ! GMT (ihr;in hours), latitude (lat), and longitude (lon).
  !
  ! Sun angle of pi/2 radians (90 deg) = sun at or below horizon i.e. dark
  ! Sun angle of 0 radians (0 deg)     = sun directly overhead.
!------------------------------------------------------------------------------

  IMPLICIT NONE 
   
!-------------------------------------------------
! Declare variables in the interface
!------------------------------------------------- 
  REAL   , INTENT(IN)    :: lat          ! Latitude (deg); lat > 0 means
                                         !                 North of Equator
  REAL   , INTENT(IN)    :: lon          ! Longitude (deg E, 0 <= lon < 360); 
  REAL   , INTENT(IN)    :: rhr          ! decimal Hour of Julian Day  
  INTEGER, INTENT(IN)    :: julianDay    ! Julian Day GMT
  
!-------------------------------------------------
! Local variables
!------------------------------------------------- 
  REAL   :: costmp
  REAL   :: sdec
  REAL   :: rad
  REAL   :: rlat 
  REAL   :: rsdec   
  REAL   :: rthez 
  REAL   :: sunangle  ! returned as RESULT 
  REAL   :: tc  
  REAL   :: thez
  REAL   :: xha


#ifdef DEBUG
write(6,*) "In calc_solar_zenith"
#endif



!---------------------------------------------------
!  Compute solar declination angle
  rad = 180.0/ (4.0 * atan(1.0))             ! rad = radian to degrees 
                                             !       conversion factor
  thez = 360.0*((FLOAT(julianDay-1)) )/365. ! thez = theta zero orbital 
                                                 !        position (degs)
  rthez = thez/rad;  ! in radians

!  sdec = solar declination angle in degrees
   sdec = 0.396372                                                    &
   &    -22.91327 *cos(rthez)     + 4.02543*sin(rthez)                & 
   &    - 0.387205*cos(2.0*rthez) + 0.051967*sin(2.0*rthez)           &
   &    - 0.154527*cos(3.0*rthez) + 0.084798*sin(3.0*rthez)           
   rsdec = sdec/rad; ! in radians

#ifdef DEBUG
!write(6,*) "In calc_solar_zenith,rsdec,rad=",rsdec,rad
#endif

!  Time correction (tc) for solar hour angle
   tc = 0.004297                                                      &
   &  + 0.107029*cos(rthez)     - 1.837877*sin(rthez)                 &
   &  - 0.837378*cos(2.0*rthez) - 2.342824*sin(2.0*rthez)

!  xha = solar hour angle in degrees
   xha = (rhr-12.0 )*15.0 + lon + tc    
     
   IF (xha > 180.0) THEN
      xha = xha - 360.0
   ENDIF 
   IF (xha > -180.0) THEN
      xha = xha + 360.0
   ENDIF 

!  Sun zenith angle sunangle. 
! Cosine(sunangle) == costmp
! rsdec            == sun declination angle
! xha              == local hour angle  

   rlat = lat/rad
   costmp = sin(rlat)*sin(rsdec)                                      & 
   &      + cos(rlat)*cos(rsdec) * cos(xha/rad)

! The next block if is for failsafe purposes
   IF (abs(costmp) > 1.0) THEN
       IF (costmp > 0.0) THEN
           costmp = 1.0
       ELSE
           costmp = -1.0
       ENDIF 
   ENDIF 

   sunangle = acos(costmp)  ! solar zenith angle in radians. Note that
                            ! 2.0*atan(1.0)= Pi/2 radians= 90 deg

   IF (sunangle > (2.0*atan(1.0))) THEN
       sunangle = 2.0*atan(1.0);
   ENDIF 

#ifdef DEBUG
write(6,*) "End calc_solar_zenith"
#endif

 
   RETURN  
   END FUNCTION calc_solar_zenith   
   Subroutine Call_IOP_PAR(                                            &
		 & PARsurf     , sun_zenith,                           &
                 & CDOM_k      , totChl    ,                           &
                 & OM1_A       , OM1_Z     ,                           &
                 & OM1_R       , OM1_BC    ,                           &
                 & depth, dz, nz, d_sfc,                                   &
	         & PARdepth  )   

!-----------------------------------------------------------------------
! Code is based on Brad Penta's code
!---------------------------------------------------------------------
  USE cgem, only: km,aw490,astar490,astarOMA,astarOMZ,astarOMR,astarOMBC 

  IMPLICIT NONE 

!                 & d(nz), dz, nz, d_sfc,                          &

!----------------------------------------------------------------------------
  real   , intent(in)    :: PARsurf      ! Irradiance just below sea surface
  real   , intent(in)    :: sun_zenith   ! Angle of the sun
  real   , intent(in)    :: totChl(km)  ! total Chl-a (mg/m3)
  real   , intent(in)    :: CDOM_k(km)  ! CDOM (ppb) 
  real   , intent(in)    :: OM1_A(km)   ! Concentration of particulate
  real   , intent(in)    :: OM1_Z(km)   ! Concentration of particulate
                                         ! fecal pellets (g/m3)
  real   , intent(in)    :: OM1_R(km)   ! Concentration of particulate
                                         ! river generated SPM (g/m3)
  real   , intent(in)    :: OM1_BC(km)  ! Concentration of particulate
                                         ! initial and boundary condition generated SPM (g/m3)
  real   , intent(in)    :: depth        ! depth at bottom of cell k from surface
  real   , intent(in)    :: dz(km)      ! depth of cell
  real   , intent(in)    :: d_sfc(km)   ! depth at center of cell k from surface
  integer, intent(in)    :: nz          ! Total number of layers in water column
  real  , intent(out)    :: PARdepth(km)    ! PARdepth(k) is Par, the
                                             ! visible irradiance
                                             ! at the middle of layer k.
                                             ! (quanta/cm**2/sec)

!----------------------------------------------------------------
! Calculate absorption (490 nm) components: seawater, chl, SPM from rivers, CDOM,
! detritus (dead cells), fecal pellets ...
      real Chla_tot(km), CDOM_tot(km), OM1A_tot(km), OM1Z_tot(km), OM1R_tot(km), OM1BC_tot(km), CDOM(km)
      real Chla_mass(km), CDOM_mass(km), OM1A_mass(km), OM1Z_mass(km), OM1R_mass(km), OM1BC_mass(km)
      real a490_mid, aSw_mid, aChl490_mid, aCDOM490_mid, bbChl490_mid, bb490_mid
!      real a490_bot, aSw_bot, aChl490_bot, aCDOM490_bot, bbChl490_bot, bb490_bot
      real aOM1A490_mid, aOM1Z490_mid, aOM1R490_mid, aOM1BC490_mid
!      real aOM1A490_bot, aOM1Z490_bot, aOM1R490_bot, aOM1BC490_bot
      real cell_depth  !, bd_km1
      integer :: k  

#ifdef DEBUG
write(6,*) "Begin Call_IOP_Par"
write(6,*) "PARsurf,sun_zenith,totChl",PARsurf,sun_zenith,totChl
write(6,*) "nz,km",nz,km
#endif


! First, convert CDOM(ppb) into CDOM, a490 (m-1)
! Once the CDOM (QSE ppb) is in the model domain, we advect and mix using the same 
! routine as for other dissolved constituents. However, to use the CDOM in the light 
! attenuation models, we need to calculate a490 (Penta et al. 2008). 
! 1) convert CDOM(QSE ppb) back to a312: a312 = (CDOM(QSE ppb)-0.538)/2.933 
! 2) convert a312 to a490: a490 = a312*exp(-0.016*(490-312)), where here S = 0.016 
! (mean value from D'Sa and DiMarco (2008)
   do k=1,nz
      CDOM(k) = (CDOM_k(k) - 0.538)/2.933 !ppb to a312
      CDOM(k) = CDOM(k) * exp(-0.016*(490.-312.))
      CDOM(k) = AMAX1(CDOM(k),0.)
   enddo 

#ifdef DEBUG
write(6,*) "Call_IOP_Par: initialized CDOM"
#endif


!Initialize counters for Chla, CDOM, and detritus:
   Chla_tot = 0.
   CDOM_tot = 0.
   OM1A_tot = 0.
   OM1Z_tot = 0.
   OM1R_tot = 0.
   OM1BC_tot = 0.
!   bd_km1 = 0.

!Mass in each cell at layer k (area of volume part cancels out)
!The unit is mg[mmol] / m2
   do k=1,nz
!      cell_depth = bottom_depth(k) - bd_km1
!      bd_km1 = bottom_depth(k)
      cell_depth = dz(k)
      Chla_mass(k) = totChl(k)*cell_depth
      CDOM_mass(k) = CDOM(k)*cell_depth
      OM1A_mass(k) = OM1_A(k)*cell_depth
      OM1Z_mass(k) = OM1_Z(k)*cell_depth
      OM1R_mass(k) = OM1_R(k)*cell_depth
      OM1BC_mass(k) = OM1_BC(k)*cell_depth
   enddo


#ifdef DEBUG
write(6,*) "Call_IOP_Par:cell mass"
#endif

!Mass from surface to center of cell at layer k
!Is the sum of the mass of all previous k layers plus 
!half of the current k layer 
!Concentration is that divided by the distance
!from the surface to the center of cell at layer k
!(Division by d_sfc is in the next loop)  
      Chla_tot(1) = 0.5*Chla_mass(1)
      CDOM_tot(1) = 0.5*CDOM_mass(1)
      OM1A_tot(1) = 0.5*OM1A_mass(1)
      OM1Z_tot(1) = 0.5*OM1Z_mass(1)
      OM1R_tot(1) = 0.5*OM1R_mass(1)
      OM1BC_tot(1) = 0.5*OM1BC_mass(1)
   do k=2,nz
      Chla_tot(k)  = 0.5*Chla_mass(k) + SUM(Chla_mass(1:k-1))
      CDOM_tot(k)  = 0.5*CDOM_mass(k) + SUM(CDOM_mass(1:k-1))
      OM1A_tot(k)  = 0.5*OM1A_mass(k) + SUM(OM1A_mass(1:k-1))
      OM1Z_tot(k)  = 0.5*OM1Z_mass(k) + SUM(OM1Z_mass(1:k-1))
      OM1R_tot(k)  = 0.5*OM1R_mass(k) + SUM(OM1R_mass(1:k-1))
      OM1BC_tot(k) = 0.5*OM1BC_mass(k)+ SUM(OM1BC_mass(1:k-1))
   enddo

#ifdef DEBUG
write(6,*) "Call_IOP_Par:total mass"
#endif

#ifdef DEBUG
write(6,*) "Call_IOP_Par:depth,dz,d_sfc:",depth,dz,d_sfc
#endif


   do k=1,nz
!Calculate absorption coefficients:
      aSw_mid = aw490  !Sea water absorption at mid cell
      aChl490_mid = astar490 * Chla_tot(k) / d_sfc(k)        !Chla absorption at mid cell
      aCDOM490_mid = CDOM_tot(k) / d_sfc(k)        !CDOM absorption at mid cell
      aOM1A490_mid = astarOMA * OM1A_tot(k) / d_sfc(k) ! absorption at mid cell
      aOM1Z490_mid = astarOMZ * OM1Z_tot(k) / d_sfc(k) ! absorption at mid cell
      aOM1R490_mid = astarOMR * OM1R_tot(k) / d_sfc(k) ! absorption at mid cell
      aOM1BC490_mid = astarOMBC * OM1BC_tot(k) / d_sfc(k) ! absorption at mid cell
      a490_mid = aSw_mid + aChl490_mid + aCDOM490_mid + aOM1A490_mid + aOM1Z490_mid + aOM1R490_mid + aOM1BC490_mid
!Calculate backscattering coefficients:
      bbChl490_mid = 0.015 * (0.3*((Chla_tot(k) / d_sfc(k))**0.62)*(550./490.)) !Chla backscatter at mid cell
      bb490_mid = bbChl490_mid !Only Chla backscatters for now
! Calculate PAR at depth
      !Why would we check if a490_mid=0??
      !If it is zero, stop.
      !  if(a490_mid.le.0) then
      !     write(6,*) k,CDOM(k),CDOM_k(k)
      !     write(6,*) "a490_mid.le.0, =",a490_mid,aSw_mid,aChl490_mid,aCDOM490_mid
      !     !stop
      !  endif
      call IOP_PARattenuation(a490_mid, bb490_mid, PARsurf, sun_zenith, d_sfc(k), PARdepth(k)) 

#ifdef DEBUG
write(6,*) "Just called IOP_Paratt, PARdepth=:",k,PARdepth(k)
#endif

   enddo

! Calculate PAR at sea bottom
!      aSw_bot = aw490  !Sea water absorption at bottom of cell
!      aChl490_bot = astar490 * (Chla_tot(nz)+0.5*Chla_mass(nz)) / depth !Chla absorption at bottom
!      aCDOM490_bot = CDOM_tot(nz)+(0.5*CDOM_mass(nz)) / depth !CDOM absorption at bottom
!      aOM1A490_bot = astarOMA * (OM1A_tot(nz)+0.5*OM1A_mass(nz)) / depth    !A absorption at bottom
!      aOM1Z490_bot = astarOMZ * (OM1Z_tot(nz)+0.5*OM1Z_mass(nz)) / depth !FP absorption at bottom
!      aOM1R490_bot = astarOMR * (OM1R_tot(nz)+0.5*OM1R_mass(nz)) / depth !SPM absorption at bottom
!      aOM1BC490_bot = astarOMbc * (OM1BC_tot(nz)+0.5*OM1BC_mass(nz)) / depth !INIT/BC absorption at bottom
!      a490_bot = aSw_bot + aChl490_bot + aCDOM490_bot + aOM1A490_bot + aOM1Z490_bot + aOM1R490_bot + aOM1BC490_bot
!      bbChl490_bot = 0.015 * (0.3*(((Chla_tot(nz)+0.5*Chla_mass(nz)) / depth)**0.62)*(550./490.))
!!Chla backscatter at bottom
!      bb490_bot = bbChl490_bot !Only Chla backscatters for now
!      call IOP_PARattenuation(a490_bot, bb490_bot, PARsurf, sun_zenith, d_sfc(nz), PARbot)

   END Subroutine Call_IOP_PAR
!----------------------------------------------------------------------
!*****************************************************************************************
      subroutine IOP_PARattenuation(a490, bb490, PARsurf, sun_zenith, d_sfc, PARdepth)

! Penta et al., 2008 PAR penetration model: based on the IOP (inherent optical properties)
! absorption and backscatter at 490nm and the zenith angle of the sun; modified from Lee 
! et al., 2005 PAR penetration model: based on satellite absorption and backscater at 490nm

!    Lee, Z., K. Du, R. Arnone, S. Liew, and B. Penta, .Penetration of solar radiation
!         in the upper ocean: a numerical model for oceanic and coastal waters,. 
!         J. Geophys. Res. 110, C09019 doi:10.1029/2004JC002780 (2005).
!    Penta, B., Z. Lee, R.M. Kudela, S.L. Palacios, D.J. Gray, J.K. Jolliff, and 
!         I.G. Shulman, .An underwater light attenuation scheme for marine ecosystem 
!         models., Optics Express, 16, 16582-16591 (2008).

! Absorption (a490) is the total absorption at 490 nm calculated from all of the 
! constituents of the model that absorb light (Seawater, chlorophyll, CDOM, detritus, etc.
! Backscatter (bb490) is similar 

! Use the average value from the sea-surface to the depth of the calculation 
! (these calculations moved to a new subroutine(s) and the a490 and bb490 will be passed 
! into this subroutine  

! PAR (photosynthetically active radiation) just below the sea surface is used as a 
! starting value to be multiplied by the attenuation factor computed in this subroutine. 
! The starting value does not affect any of these calculations

! We generally use a factor of 0.4815 to represent the loss of PAR passing across the
! Air/Sea interface - this is included already in the data 

! Define coefficients for light attenuation model       
      real a490, alpha0, alpha1, alpha2, bb490, chi0, chi1, chi2, d_sfc
      real k1, k2, kpar, PARdepth, PARsurf, sun_zenith, zeta0, zeta1
      real zeta2 
     
! sun_zenith => Solar Zenith Angle in radians for time and location
! d_sfc => depth of center of cell from elevated sea surface        

! Set values of coefficients for light attenuation model 
! (Lee et al., 2005; Penta et al., 2008)       
      parameter(alpha0=0.090, alpha1=1.465, alpha2=-0.667, chi0=-0.057, chi1=0.482, & 
     &          chi2=4.221, zeta0=0.183, zeta1=0.702, zeta2=-2.567)
     
! Calculate the attenuation (Equation 2 Penta et al., 2008 errata)
      k1 = ((chi0 + chi1*(a490**0.5) + chi2*bb490) & 
     &     * (1 + alpha0 * sin(sun_zenith)))
     
      k2 = ((zeta0 + zeta1*a490 + zeta2*bb490) &
     &     * (alpha1 + alpha2 * cos(sun_zenith) ))
     
      kpar = k1 + (k2 / (1+(d_sfc))**0.5)
     
      PARdepth = PARsurf * exp(-(kpar*(d_sfc)))

      return
      end
!*****************************************************************************************
module cgem_light

use grid, only: km,d,dz,d_sfc,iYrS,lat,lon
use cgem, only: which_irradiance,Kw,Kchla,Kspm,Kcdom,aRadSum
use date_time

implicit none

contains

subroutine calc_PARdepth( TC_8,PARSurf,S_k,Chla_k,CDOM_k,OM1A_in,OM1Z_in,      &
 &                         OM1R_in,OM1BC_in,PARdepth_k,PARbot,aRadSum_k )

!---------------------------------------------
! Interface variables
!---------------------------------------------------------------------
    !Inputs
    integer(kind=8), intent(in)      :: TC_8              ! Model time (seconds from iYrS)
    real, dimension(km), intent(in)  :: S_k               ! Salinity (psu)
    real, dimension(km), intent(in)  :: Chla_k            ! Total amount of Chl-a in all the
                                                          !  phytoplankton species (mg/m3) per cell
    real, dimension(km), intent(in)  :: OM1A_in, OM1Z_in  !  
    real, dimension(km), intent(in)  :: OM1R_in, OM1BC_in ! POC in g/m3
    real, dimension(km), intent(in)  :: CDOM_k            ! CDOM, ppb
    real,                intent(in)  :: PARsurf           ! Irradiance just below the sea surface (quanta/cm2/s) 
    !Outputs
    real, dimension(km), intent(out) :: aRadSum_k      ! Update running total of current day's irradiance
    real, dimension(km), intent(out) :: PARdepth_k(km) ! Irradiance at center of layer k (quanta/cm2/s)
    real,                intent(out) :: PARbot         ! Irradiance at sea floor (quanta/cm2/s)
!---------------------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------
    real, dimension(km)  :: OM1A_k, OM1Z_k    !  
    real, dimension(km)  :: OM1SPM_k, OM1BC_k   ! POC in g/m3

    integer        ::  k, isp, isz ! Loop indicies, isp/isz is for phytoplankton/zooplankton species
    integer        ::  Is_Day            ! Switch for day/night for phytoplankton nutrient uptake only, Is_Day=0 means night
!------------------------------------ 
! Time variables  
    real, parameter :: one_d_365  = 1.0/365.0 ! Convert 1/yr to 1/day
    real, parameter :: OneD60     = 1.0/60.0  ! Convert 1/min to 1/sec
    real            :: HrTC          ! Decimal hour of day
    integer         :: iYrTC, iMonTC, iDayTC, iHrTC, iMinTC, iSecTC !Time variables
    integer         :: julianDay     ! Holds Julian Day
!-----------------------------------------------------------------------
! Variables needed for light routine and calc_Agrow
    real    :: SunZenithAtm       ! Solar beam zenith angle
    real    :: calc_solar_zenith  ! Function, calculates solar beam zenith angle
    real    :: Katt               ! Attenuation coefficient for Irradiance model 2 
    real    :: tmpexp             ! Intermediate calculation
    real    :: PARbotkm1          ! Irradiance at bottom of layer k-1 (quanta/cm2/s)
    real    :: PARtopk            ! Irradiance at top of layer k (quanta/cm2/s)
    real, parameter :: RADCONV = 1./6.0221413*1.e-19 ! Convert quanta/cm2/s to mol/m2/s:
                                               !  = quanta/cm2/s * 1 mol/Avogadro# * 10,000cm2/m2
                                               !  = (1/6.022e23) * 1.0e4 = (1./6.022)e-23 * 1.0e4
                                               !  = (1./6.0221413)*1.e-19
                                       !  phytoplankton species (mg/m3) per cell
    real, parameter :: C_cf  = 12.0E-3    ! C conversion factor (mmol-C/m3 to g-C/m3) 
!-----------------------------------------------------------------------
#ifdef DEBUG_PAR
write(6,*) "Begin calc_light:calc_PARdepth, TC_8",TC_8
#endif

!-----------------------------------------------------------------
!   Begin main ij loop for the biogeochemistry 
!   calculations at time-level istep
!-----------------------------------------------------------------
 !---------------------------------------------------------
 ! Calculate and convert variables needed for light routine
 !---------------------------------------------------------
      do k = 1, km

#ifdef DEBUG_PAR
write(6,*) "init loop, km",k
#endif
         !-------------------------------------------------
         ! -- Convert units for light model 
         !    C_cf == conversion factor (mmol-C/m3 to g-C/m3) 
         !-----
         ! Organic Matter from dead phytoplankton (mmol/m3) 
         !            converted to equivalent (g carbon/m3)
           OM1A_k(k)  = OM1A_in(k) * C_cf        
         !-----
         ! Organic Matter from fecal pellets      (mmol/m3)
         !            converted to equivalent (g carbon/m3)
           OM1Z_k(k) = OM1Z_in(k) * C_cf   
         !-----
         ! Organic Matter from rivers            
         !  Suspended Particulate Matter (SPM)    (mmol/m3) 
         !            converted to equivalent (g carbon/m3)
         ! There is 1.8% Organic Matter in SPM originating from the rivers.
           OM1SPM_k(k) = OM1R_in(k) * C_cf / 0.018
         !-----
         ! Organic Matter from boundary conditions(mmol/m3) 
         !            converted to equivalent (g carbon/m3)
           OM1BC_k(k)  = OM1BC_in(k) * C_cf 

      enddo ! End  "do k = 1, km" loop


!----------------------------------------------------------------------
! Execute the desired atmospheric light model.  To calculate PARsurf,
! the effect amount of downward spectrally integrated irradiance 
! just below the sea surface.  'Rad' is just above sea surface. 
!----------------------------------------------------------------------
#ifdef DEBUG_PAR
write(6,*) "In cgem_light:calc_PARdepth, calculate date_time is next, iYrS=",iYrS
#endif

 ! First calculate the Julian(GMT) model year (iYrTC), month (iMonTC), 
 ! day (iDayTC), hour (iHrTC), minute (iMinTC), and second (iSecTC) 
 ! associated with the midpoint of the present timestep istep TC_8

      CALL DATE_TIMESTAMP( iYrS, TC_8, &
                           iYrTC, iMonTC, iDayTC, iHrTC, iMinTC, iSecTC )

 ! Calc HrTC, the decimal hour of day
       HrTC = real(iHrTC,4) + OneD60*iMinTC + OneD60*iSecTC

 ! Now calculate the Julian Day associated with model time TC_8
      julianDay = JDAY_IN_YEAR(iYrTC, iMonTC, iDayTC)

 ! Now calculate SunZenithAtm, the solar beam zenith angle in radians
 ! for a given GMT Julian day, hour, longitude and latitude 
     SunZenithAtm = calc_solar_zenith(lat, lon,  HrTC, julianDay )

!----------------------------------------------------------------------------
! Execute the desired underwater light model to calculate the 1-D radiation
! arrays PARdepth_k, Esed and PAR_percent_k radiation arrays for
! vertical grid column (i).
!
! PARdepth_k(k) is the downward irradiance (photons/cm2/sec) at the middle
!                    of cell(k).
!
! PARbot is the downward irradiance (photons/cm2/sec) at the sea bottom
!
! PAR_percent_k(k)    is the % of incoming irradiance PARsurf that PARdepth_k(k)
!                 represents. PARsurf is the downward irradiance
!                 (photons/cm2/sec) just below the sea surface.
!-------------------------------------------------------------------------
         select case (which_irradiance)

                 !--------------------------------------------
         case (1)! Upgraded form of the underwater light model
                 ! developed by Brad Penta of NRL is used
                 !--------------------------------------------

#ifdef DEBUG_PAR
  write(6,*) "In calc_PARdepth, calling Brad's light model, km,PARsurf=",km,PARsurf
#endif

                call Call_IOP_PAR(                    &
                 & PARsurf, Chla_k, CDOM_k,           &
                 & OM1A_k, OM1Z_k, OM1SPM_k, OM1BC_k, &
                 & d(km), dz, km, d_sfc,              &
                 & PARdepth_k, PARbot)

#ifdef DEBUG_PAR
     write(6,*) "In cgem_light:calc_PARdepth."
     write(6,*) "PARdepth_k,PARbot=",PARdepth_k
     write(6,*) "PARbot=",PARbot
#endif

                 !-------------------------------------------------
         case (2)! Upgraded form of the original underwater light
                 ! model of Pete Eldridge is used. Now accounts for
                 ! light attenuation in each k layer rather than for the whole
                 ! mixed layer as in the Eldridge & Roelke(2010) code 
                 !-------------------------------------------------
                 PARbotkm1 = PARsurf             ! initialize at top of
                                                 ! column i., i.e.
                                                 ! at bottom of layer "zero".
                 do k = 1, km
                   !Calculate attenuation coefficient
                     Katt    = Kw                                                   &
                     &       + Kchla * Chla_k(k)                                    &
                     &       + Kspm  * (OM1SPM_k(k)+OM1A_k(k)+OM1Z_k(k)+OM1BC_k(k)) &
                     &       + Kcdom * CDOM_k(k)                                    &
                     &       + (((0.0022*S_k(k))-0.158)*S_k(k)+3.03)
                     PARtopk = PARbotkm1         ! irradiance at top of
                                                 ! layer k is same as
                                                 ! irradiance at bottom
                                                 ! of layer km1
                     tmpexp  = exp(-0.5*Katt*dz(k))
                     PARdepth_k(k) = PARtopk * tmpexp    ! irradiance at middle of layer k
                     PARbot  = PARdepth_k(k) * tmpexp    ! irradiance at bottom of layer k
                     PARbotkm1 = PARbot         ! reinitialize for next top layer
                 enddo 

         case default
             write(6,*) "Error in irradiance switch",which_irradiance
             stop
         end select

!---------------------End Underwater Light Model-----------------------------------
#ifdef DEBUG_PAR
write(6,*) "In cgem_light:calc_PARdepth, Finished Underwater Light Model"
#endif

! Update running total of current day's irradiance
         aRadSum_k(:) = aRadSum(:) + PARdepth_k(:)

#ifdef DEBUG_PAR
write(6,*) "Running total of current day's irradiance:"
write(6,*) "aRadSum_k=",aRadSum_k
#endif


   return
   end subroutine calc_PARdepth  
!---------------------------------------------------------------------- 

end module
