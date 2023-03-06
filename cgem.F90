!!======================================================================     
subroutine cgem(it)

!======================================================================
  use cgem_vars
  use schism_glbl, only: rkind

!SCHISM global variables
!  use schism_glbl, only : rkind,dt,ne,nea,npa,nvrt,bdy_frc,idry_e,kbe,ze,&
!      & tr_el,xlon_el,ylat_el,xlon,ylat,irange_tr,ntrs,ielg,iplg,elnode,srad,&
!      & su2,sv2,elside,iegl,eta2,i34,windx,windy,wsett,flx_sf,flx_bt,rnday

  implicit none 

! it= iteration
  integer,intent(in) :: it
  integer :: id,nz,k,isp

!From the hydro
    real(rkind), dimension(km) :: S, T  ! Salinity and Temperature(Celsius)
!-------------------------------------------------------------------------
! Phytoplankton parameters
! Phytoplankton uptake and growth
    real(rkind), dimension(nospA,km) :: A      ! Phytoplankton number density (cells/m3)
    real(rkind)    :: Agrow                       ! Phytoplankton growth (cells/m3/d)
    real(rkind), dimension(nospA,km) :: Agrow_k  ! Phytoplankton growth (cells/m3/d)
    real(rkind)    :: uA                 ! Specific growth rate (1/d)
    real(rkind)    :: uA_k(km,nospA)    ! Specific growth rate (1/d)
    real(rkind)    :: uN_k(km,nospA)    ! Nitrogen Limited growth rate (1/d)
    real(rkind)    :: uP_k(km,nospA)    ! Phosphorus limited growth rate (1/d)
    real(rkind)    :: uE_k(km,nospA)    ! Light limited growth rate (1/d)
    real(rkind)    :: uSi_k(km,nospA)   ! Silica limited growth rate (1/d)
    real(rkind)    :: f_Qn(nospA)        ! Quota model for N
    real(rkind)    :: f_Qp(nospA)        ! Quota model for P
    real(rkind)    :: Qn(nospA,km) ! Phytoplankton Nitrogen Quota (mmol-N/cell)
    real(rkind)    :: Qp(nospA,km) ! Phytoplankton Phosphorus Quota (mmol-P/cell)
    real(rkind)    :: vN    ! Phytoplankton uptake rate of Nitrogen (mmol-N/cell/d)
    real(rkind)    :: vP    ! Phytoplankton uptake rate of Phosphorus (mmol-P/cell/d)
    real(rkind)    :: vSi   ! Phytoplankton uptake rate of Silica (mmol-Si/cell/d)
    real(rkind)    :: AupN  ! Total Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
    real(rkind)    :: AupP  ! Total Phytoplankton uptake of Phosphorus (mmol-P/m3/d)
    real(rkind)    :: AupSi ! Total Phytoplankton uptake of Silica (mmol-Si/m3/d)
    integer :: RLN   ! Rate Limiting Nutrient of N, P, and Si
 ! Monod equations for phytoplankton
    real(rkind), dimension(nospA)    :: monodN  !Monod term in nitrogen uptake
    real(rkind), dimension(nospA)    :: monodP  !Monod term in phosphorus uptake
    real(rkind), dimension(nospA)    :: monodSi !Monod term in Si uptake
    real(rkind)    :: Ntotal   ! Total N (mmol/m3)
 ! Phytoplankton nutrient loss
    real(rkind), dimension(nospA)    :: Amort ! Dead phytoplankton (cells/m3/day)
    real(rkind), dimension(nospA)    :: AexudN_A    ! Phytoplankton group exudation (mmol-N/cell/d)
    real(rkind), dimension(nospA)    :: AexudP_A    ! Phytoplankton group exudation (mmol-P/cell/d) 
    real(rkind)    :: AexudN          ! Sum of Exudation of N from all phytoplankton groups (mmol-N/m3/d)
    real(rkind)    :: AexudP          ! Sum of Exudation of P from all phytoplankton groups (mmol-P/m3/d)
    real(rkind)    :: Aresp           ! Total respiration from a phytoplankton group (cells/m3/d)
    real(rkind)    :: Aresp_k(nospA,km) ! Total respiration from a phytoplankton group (cells/m3/d)
    real(rkind)    :: ArespC          ! Phytoplankton equivalent carbon loss from respiration (mmol-C/m3/d)
!---------------------------------------------------------
! Variables needed for light routine and calc_Agrow
    real(rkind)   :: SunZenithAtm       ! Solar beam zenith angle
    real(rkind)   :: calc_solar_zenith  ! Function, calculates solar beam zenith angle
    real(rkind)   :: Katt               ! Attenuation coefficient for Irradiance model 
    real(rkind)   :: tmpexp             ! Intermediate calculation
    real(rkind)   :: PARbotkm1          ! Irradiance at bottom of layer k-1 (quanta/cm2/s)
    real(rkind)   :: PARtopk            ! Irradiance at top of layer k (quanta/cm2/s)
    real(rkind)   :: PARsurf            ! Irradiance just below the sea surface (quanta/cm2/s) 
    real(rkind)   :: PARbot             ! Irradiance at sea floor (quanta/cm2/s)
    real(rkind)   :: PARdepth_k(km)    ! Irradiance at center of layer k (quanta/cm2/s)
    real(rkind)   :: PAR_percent_k(km) ! Percent irradiance at center of layer k (quanta/cm2/s)
    real(rkind)   :: aDailyRad_k(km), aRadSum_k(km)
    real(rkind), parameter :: RADCONV = 1./6.0221413*1.e-19 ! Convert quanta/cm2/s to mol/m2/s:
                                               !  = quanta/cm2/s * 1
                                               !  mol/Avogadro# * 10,000cm2/m2
                                               !  = (1/6.022e23) * 1.0e4 =
                                               !  (1./6.022)e-23 * 1.0e4
                                               !  = (1./6.0221413)*1.e-19
    real(rkind), dimension(km) :: Chla_tot  ! Total amount of Chl-a in all the
                                     !  phytoplankton species (mg/m3) per cell
    real(rkind), dimension(nospA,km) :: Chl_C_k     ! Chl:C
    real(rkind), dimension(km) :: OM1A_k, OM1Z_k, OM1SPM_k, OM1BC_k !POC in g/m3
    real(rkind), dimension(km) :: CDOM_k    ! CDOM, ppb
    real(rkind), dimension(km) :: N_k       ! Nitrogen, mmol/m3
    real(rkind), dimension(km) :: P_k       ! Phosphorus, mmol/m3
    real(rkind), dimension(km) :: Si_k      ! Silica, mmol/m3
    real(rkind), dimension(km) :: S_k, T_k  ! Salinity and Temperature(Celsius)
    real(rkind), parameter :: C_cf  = 12.0E-3    ! C conversion factor (mmol-C/m3 to g-C/m3) 

!------------------------------------------------------------------
! Zooplankton parameters
 !Zooplankton uptake and growth
    real(rkind), dimension(nospZ,km)   :: Z      ! Zooplankton number density (indv./m3)
    real(rkind), dimension(nospZ)       :: optNP    ! Optimal nutrient ratio for zooplankton
    real(rkind), dimension(nospZ)       :: Zgrow    ! Zooplankton growth (indv./m3/d)
    real(rkind), dimension(nospA,nospZ) :: Zgrazvol ! Grazing rate in units of biovolume (um3/m3/d)
    real(rkind), dimension(nospA,nospZ) :: ZgrazA   ! Zooplankton grazing of phytoplankton (cells/m3/d)
    real(rkind), dimension(nospA)       :: ZgrazA_tot ! Total zooplankton grazing of phytoplankton (cells/m3/d)
    real(rkind), dimension(nospZ)       :: ZgrazN   ! Zooplankton grazing uptake of Nitrogen (mmol-N/m3/d)
    real(rkind), dimension(nospZ)       :: ZgrazP   ! Zooplankton grazing uptake of Phosphorus (mmol-P/m3/d)
    real(rkind), dimension(nospZ)       :: ZgrazC   ! Zooplankton grazing uptake of Carbon (mmol-C/m3/d)
    real(rkind), dimension(nospZ)       :: ZinN     ! Zooplankton ingestion of Nitrogen (mmol-N/m3/d)
    real(rkind), dimension(nospZ)       :: ZinP     ! Zooplankton ingestion of Phosphorus (mmol-P/m3/d)
    real(rkind), dimension(nospZ)       :: ZinC     ! Zooplankton ingestion of Carbon (mmol-C/m3/d)
 !Monod equations for zooplankton ingestion of phytoplankton
    real(rkind), dimension(nospA,nospZ) :: monodZ   ! Monod term for zooplankton grazing
    real(rkind)                         :: Abiovol  ! Algae biovolume vector (um3/m3)
    real(rkind), dimension(nospA,nospZ) :: top_A    ! Monod numerator value for phytoplankton group
    real(rkind), dimension(nospA,nospZ) :: bottom_A ! Monod Denominator value for phytoplankton group
    real(rkind), dimension(nospZ)       :: bottom   ! Sum of Monod Denominator value for all phytoplankton groups
 !Zooplankton nutrient loss
    real(rkind), dimension(nospZ)       :: Zresp    ! Zooplankton respiration (individuals/m3/d)
    real(rkind)                         :: ZrespC   ! Carbon loss from zooplankton respiration (mmol-C/m3/day)
    real(rkind), dimension(nospZ)       :: ZunC     ! Unassimilated ingested Carbon (mmol-C/m3/d)
    real(rkind), dimension(nospZ)       :: ZunN     ! Unassimilated ingested Nitrogen (mmol-N/m3/d)
    real(rkind), dimension(nospZ)       :: ZunP     ! Unassimilated ingested Phosphorus (mmol-P/m3/d)
    real(rkind), dimension(nospZ)       :: ZunSi    ! Unassimilated ingested Silica (mmol-Si/m3/d)
    real(rkind), dimension(nospZ)       :: Zmort    ! Dead zooplankton (individuals/m3/d)
    real(rkind) :: ZmortC(nospZ), ZmortC_tot        ! Carbon released from dead zooplankton (mmol-C/m3/d)
    real(rkind) :: ZmortN(nospZ), ZmortN_tot        ! Nitrogen released from dead zooplankton (mmol-N/m3/d)
    real(rkind) :: ZmortP(nospZ), ZmortP_tot        ! Phosphorus released from dead zooplankton (mmol-P/m3/d)
    real(rkind) :: ZslopC(nospZ), ZslopC_tot        ! Carbon lost to sloppy feeding (mmol-C/m3/d)
    real(rkind) :: ZslopN(nospZ), ZslopN_tot        ! Nitrogen lost to sloppy feeding (mmol-N/m3/d)
    real(rkind) :: ZslopP(nospZ), ZslopP_tot        ! Phosphorus lost to sloppy feeding (mmol-P/m3/d)
    real(rkind), dimension(nospZ)       :: ZexN     ! Excretion from zooplankton (mmol-N/m3/d)
    real(rkind), dimension(nospZ)       :: ZexP     ! Excretion from zooplankton (mmol-P/m3/d)
    real(rkind), dimension(nospZ)       :: ZegC     ! Egestion from zooplankton (mmol-C/m3/d)
    real(rkind), dimension(nospZ)       :: ZegN     ! Egestion from zooplankton (mmol-N/m3/d)
    real(rkind), dimension(nospZ)       :: ZegP     ! Egestion from zooplankton (mmol-P/m3/d)
    real(rkind), dimension(nospZ)       :: ZegSi    ! Egestion from zooplankton (mmol-Si/m3/d)
    real(rkind) :: OM1_Ratio, OM2_Ratio             ! Separates sloppy feeding into OM1 and OM2
!---------------------------------------------------------------------------
! reaction and Nitrification subroutine variables
    real(rkind), dimension(km)    :: NO3, NH4, PO4, DIC, O2, CDOM, Si         ! Nutrient input to subroutines
    real(rkind), dimension(km)    :: OM1_A, OM1_Z, OM1_R, OM2_A, OM2_Z, OM2_R ! OM inputs to subroutines
    real(rkind), dimension(km)    :: OM1_BC, OM2_BC                           ! OM inputs to subroutines
    real(rkind), dimension(km)    :: ALK, Tr                                  ! Alkalinity and Tracer


 write(6,*) "nea",nea
  nz = km 

!Loop over elements
!nea - local number of elements in augmented subdomain (ne+neg)
  do id=1,nea
!    !Exit loop if the element is dry
!    if(idry_e(id)==1) cycle  !element becomes dry
  
!If it is not dry, calculate away
!Other than that, I don't know what this is doing    
!    call update_cosine_vars(i)
      do k = 1, nz
!      Variables from SCHISM code
         ! Temperature (celsius) and Salinity in columns
           T(k)     = 16. !T(i,j,k)
           S(k)     = 30. !S(i,j,k)

!      CGEM variables
       do isp = 1, nospA
               A(k,isp) = ff(id,k,isp) ! Phytoplankton in group isp, cells/m3
               Qn(k,isp) = ff(id,k,iQn(1)-1+isp)
               Qp(k,isp) = ff(id,k,iQp(1)-1+isp)
             enddo
         !Save Zooplanton to k array
             do isp = 1,nospZ
                Z(isp,k) = ff(id,k,iZ(isp)) ! Zooplankton in group isp, ind./m3
             enddo

           NO3(k)     = ff(id,k,iNO3)
           NH4(k)     = ff(id,k,iNH4)
           PO4(k)     = ff(id,k,iPO4)
           DIC(k)     = ff(id,k,iDIC)
           O2(k)      = ff(id,k,iO2)
           OM1_A(k)   = ff(id,k,iOM1_A)
           OM2_A(k)   = ff(id,k,iOM2_A)
           OM1_Z(k)   = ff(id,k,iOM1_Z)
           OM2_Z(k)   = ff(id,k,iOM2_Z)
           OM1_R(k)   = ff(id,k,iOM1_R)
           OM2_R(k)   = ff(id,k,iOM2_R)
           CDOM(k)    = ff(id,k,iCDOM)
           Si(k)      = ff(id,k,iSi)
           OM1_BC(k)  = ff(id,k,iOM1_BC)
           OM2_BC(k)  = ff(id,k,iOM2_BC)
           ALK(k)     = ff(id,k,iALK)
           TR(k)      = ff(id,k,iTR)
     enddo !end loop over layers

      !loop over layers again
      do k = 1, nz


     enddo !end loop over layers


!----------------------------------------------------------------
! Get chlorophyll-a quantity per layer
!----------------------------------------------------------------
!      select case (Which_chlaC)
!      case (1) ! Use fixed C:Chla 
    do k = 1, nz
       Chla_tot(k) = 0.0
       do isp = 1, nospA
          Chla_tot(k) =  Chla_tot(k) + A(isp,k) * Qc(isp) * 12. * (1./CChla(isp))
       enddo ! isp = 1, nospA
    enddo ! k = 1, nz




  enddo  !end loop over elements


end subroutine cgem 





!Global Variables are here:  Core/schism_glbl.F90
!rkind - default real(rkind) datatype
!dt
!ne - local number of resident elements
!nea - local number of elements in augmented subdomain (ne+neg)
!npa - local number of nodes in augmented subdomain (np+npg)
!nvrt - number of vertical layers
!bdy_frc(:,:,:) - body force at prism center Q_{i,k}
!idry_e(:) - wet/dry flag
!kbe(:) - element bottom vertical indices
!ze(:,:) - z-coord (local frame - vertical up)
!tr_el(ntracers,nvrt,nea2) - tracer concentration at prism center; used as temp. storage. tr_el(ntracers,nvrt,nea2) but last index usually 
! is valid up to nea only except for TVD
!xlon_el(:) - element center lon coordinates in _degrees_
!ylat_el(:) - element center lat coordinates in _degrees_
!xlon(:) - node lon coordinates in _radians_
!ylat(:) - node lat coordinates in _radians_
!irange_tr(2,natrm) - misc. variable shared between routines
!natrm - number of tracer models at the moment (including T,S)
!ntrs(natrm)
!ielg(:) - local-to-global element index table (augmented)
!iplg(:) - local-to-global node index table (augmented)
!elnode(:,:) - element-node tables
!srad(:) - shortwave radiation
!su2(:,:) - x-component of velocity at side centers & whole levels
!sv2(:,:) - y-component of velocity at side centers & whole levels
!elside(:,:) - element-side tables
!iegl(:) - global-to-local element index table
!eta2(:) - elevation at nodes at current timestep
!i34(:) - element type (3 or 4)
!windx(:) - x-component of wind vector
!windy(:) - y-component of wind vector
!wsett(ntracers,nvrt,nea) - settling velocity (positive downward) [m/s] for each tracer
!flx_sf(:,:) - surface boundary condition \kappa*dC/dz = flx_sf (at element center)
!flx_bt(:,:) - bottom boundary condition
!rnday


!Phytoplankton population growth (cells m−3day−1) was calculated for each of three phytoplankton func-tional types (PFT): 
!1) freshwater diatoms (equivalent spherical diameter(ESD) = 15μm), 
!2) saltwater diatoms (ESD = 15μm), and 
!3) dinoflagellates (ESD = 20μm).

!Zooplankton population dynamics and phytoplankton grazing were based on Roelke (2000) and include a
!1) macrozooplankton group (ESD = 250μm) and a 
!2) microzooplankton group (ESD = 30μm). 
