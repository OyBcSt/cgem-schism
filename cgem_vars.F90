module cgem_vars

!CGEM STATE VARIABLES
use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use schism_glbl, only: rkind

implicit none

SAVE

!Grid parameters
integer :: nospA
integer :: nospZ
integer :: km
integer :: nea

!---------------------------------------------------------      
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------  
      integer, dimension(:), allocatable :: iA(:)
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
      integer, dimension(:), allocatable :: iQn(:)
!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
      integer, dimension(:), allocatable :: iQp(:)
!--------------------------------------------------------------------
!-Z: Zooplankton number density (individuals/m3);
!--------------------------------------------------------------------
      integer, dimension(:), allocatable :: iZ(:)
!-------------------------------
!-NO3; Nitrate (mmol-N/m3)
!-------------------------------
      integer :: iNO3
!--------------------------------      
!-NH4; Ammonium (mmol-N/m3)
!--------------------------------
      integer :: iNH4
!-------------------------------------------        
!-PO4: Phosphate (mmol-P/m3)
!--------------------------------------
      integer :: iPO4 
!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!---------------------------------------------------------
      integer :: iDIC 
!-------------------------------------------        
!-O2: Molecular Oxygen (mmol-O2/m3)
!------------------------------
      integer :: iO2 
!-------------------------------------------------------------
!-OM1_A: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton
!-------------------------------------------------------------
      integer :: iOM1_A
!-----------------------------------------------------------------
!-OM2_A: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton 
!------------------------------------------------------------------
      integer :: iOM2_A
!-------------------------------------------------------------
!-OM1_Z:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets.
!-------------------------------------------------------------
      integer :: iOM1_Z
!-------------------------------------------------        
!-OM2_Z:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets.
!-----------------------------------------------
      integer :: iOM2_Z
!--------------------------------------------------------------------
!-OM1_R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!--------------------------------------------------------------------
      integer :: iOM1_R
!-------------------------------------------------      
!-OM2_R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!--------------------------------------------------------------------
      integer :: iOM2_R
!-------------------------------------------
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-------------------------------------------
      integer :: iCDOM
!---------------------------------------------
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-------------------------------------------
      integer :: iSi
!--------------------------------------------------------------------
!-OM1_BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!--------------------------------------------------------------------
      integer :: iOM1_BC
!-------------------------------------------------
!-OM2_BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!--------------------------------------------------------------------
      integer :: iOM2_BC
!-------------------------------------------
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!-------------------------------------------
      integer :: iALK
!Tracer
      integer :: iTR
!Total number of state variables
      integer :: nf

!State Variable Array
      real(rkind),allocatable :: ff(:,:,:) !state variable array

!----INPUT_VARS_CGEM
!--Switches in GEM---------
integer Which_fluxes(8)
integer Which_uptake
integer Which_quota
integer Which_irradiance
integer Which_chlaC
integer Which_photosynthesis
integer Which_growth
integer Which_temperature
!--Temperature
real, allocatable :: KTg1(:)
real, allocatable :: KTg2(:)
real, allocatable :: Tref(:)
real, allocatable :: Ea(:)
real, allocatable :: N(:)
!--Optics-----------------------
real Kw
real Kcdom
real Kspm
real Kchla
!--in module LIGHT_VARS
real astar490
real aw490
real astarOMA
real astarOMZ
real astarOMR
real astarOMBC
real PARfac
!---Phytoplankton 
real, allocatable :: ediblevector(:,:)
real, allocatable :: umax(:)
real, allocatable :: CChla(:)
real, allocatable :: alpha(:)
real, allocatable :: beta(:)
real, allocatable :: respg(:)
real, allocatable :: respb(:)
real, allocatable :: QminN(:)
real, allocatable :: QminP(:)
real, allocatable :: QmaxN(:)
real, allocatable :: QmaxP(:)
real, allocatable :: Kn(:)
real, allocatable :: Kp(:)
real, allocatable :: Ksi(:)
real, allocatable :: KQn(:)
real, allocatable :: KQp(:)
real, allocatable :: nfQs(:)
real, allocatable :: vmaxN(:)
real, allocatable :: vmaxP(:)
real, allocatable :: vmaxSi(:)
real, allocatable :: aN(:)
real, allocatable :: volcell(:)
real, allocatable :: Qc(:)
real, allocatable :: Athresh(:)
real, allocatable :: mA(:)
real, allocatable :: A_wt(:)
!---Zooplankton
real, allocatable :: Zeffic(:)
real, allocatable :: Zslop(:)
real, allocatable :: Zvolcell(:)
real, allocatable :: ZQc(:)
real, allocatable :: ZQn(:)
real, allocatable :: ZQp(:)
real, allocatable :: ZKa(:)
real, allocatable :: Zrespg(:)
real, allocatable :: Zrespb(:)
real, allocatable :: Zumax(:)
real, allocatable :: Zm(:)
!---Organic Matter              
real KG1
real KG2
real KG1_R
real KG2_R
real KG1_BC
real KG2_BC
real KNH4
real nitmax
real KO2
real KstarO2
real KNO3
real pCO2
real stoich_x1R
real stoich_y1R
real stoich_z1R
real stoich_x2R
real stoich_y2R
real stoich_z2R
real stoich_x1BC
real stoich_y1BC
real stoich_z1BC
real stoich_x2BC
real stoich_y2BC
real stoich_z2BC
real KGcdom
real CF_SPM
!----Other including Boundary Conditions------------
real KH_coeff
integer Which_Outer_BC
real wt_l, wt_o
real wt_pl, wt_po
real m_OM_init,m_OM_BC,m_OM_sh
!real Stoich_x1A_init,Stoich_y1A_init,Stoich_z1A_init
!real Stoich_x2A_init,Stoich_y2A_init,Stoich_z2A_init
!real Stoich_x1Z_init,Stoich_y1Z_init,Stoich_z1Z_init
!real Stoich_x2Z_init,Stoich_y2Z_init,Stoich_z2Z_init
real KG_bot

!Light curve parameters
real, allocatable :: alphad(:) ! Initial slope of photosynthesis-irradiance curve / Vmax
real, allocatable :: betad(:)  ! Photoinhibition constant / Vmax

!Diatom array
integer, allocatable :: is_diatom(:)



contains

subroutine cgem_vars_allocate()

integer i
integer ierr
integer :: counter = 0

nospA = 3
nospZ = 2
nea = 2
km = 8

!---------------------------------------------------------      
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------  
       allocate (iA(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iA"
       do i=1,nospA
          counter = counter+1
          iA(i) = counter
       enddo
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
       allocate (iQn(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iQn"
       do i=1,nospA
          counter = counter+1
          iQn(i) = counter
       enddo
!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
      allocate (iQp(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iQp"
       do i=1,nospA
          counter = counter+1
          iQp(i) = counter
       enddo
!--------------------------------------------------------------------
!-Z: Zooplankton number density (individuals/m3);
!--------------------------------------------------------------------
      allocate (iZ(nospZ),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iZ"
       do i=1,nospZ
          counter = counter+1
          iZ(i) = counter
       enddo
!-------------------------------
!-NO3; Nitrate (mmol-N/m3)
!-------------------------------
      iNO3 = counter+1
!--------------------------------      
!-NH4; Ammonium (mmol-N/m3)
!--------------------------------
      iNH4 = counter+2
!-------------------------------------------        
!-PO4: Phosphate (mmol-P/m3)
!--------------------------------------
      iPO4 = counter+3
!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!---------------------------------------------------------
      iDIC = counter+4
!-------------------------------------------        
!-O2: Molecular Oxygen (mmol-O2/m3)
!------------------------------
      iO2 = counter+5
!-------------------------------------------------------------
!-OM1_A: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton
!-------------------------------------------------------------
      iOM1_A = counter+6
!-----------------------------------------------------------------
!-OM2_A: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton 
!------------------------------------------------------------------
      iOM2_A = counter+7
!-------------------------------------------------------------
!-OM1_Z:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets.
!-------------------------------------------------------------
      iOM1_Z = counter+8
!-------------------------------------------------        
!-OM2_Z:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets.
!-----------------------------------------------
      iOM2_Z = counter+9
!--------------------------------------------------------------------
!-OM1_R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!--------------------------------------------------------------------
      iOM1_R = counter+10
!-------------------------------------------------      
!-OM2_R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!--------------------------------------------------------------------
      iOM2_R = counter+11
!-------------------------------------------
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-------------------------------------------
      iCDOM = counter+12
!---------------------------------------------
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-------------------------------------------
      iSi = counter+13
!--------------------------------------------------------------------
!-OM1_BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!--------------------------------------------------------------------
      iOM1_BC = counter+14
!-------------------------------------------------
!-OM2_BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!--------------------------------------------------------------------
      iOM2_BC = counter+15
!-------------------------------------------
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!-------------------------------------------
      iALK = counter+16
!Tracer
      iTR = counter+17

!How many state variables
      nf = iTR

      allocate(ff(nea,km,nf),stat=ierr)
      if(ierr.ne.0) write(6,*) "error in allocating:ff"


!----allocate INPUT_VARS_CGEM

!---Phytoplankton 
allocate( ediblevector(nospZ,nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ediblevector"

allocate( umax(nospA),stat=ierr  )
if(ierr.ne.0) write(6,*) "error in allocating:umax"

allocate( CChla(nospA),stat=ierr  )
if(ierr.ne.0) write(6,*) "error in allocating:CChla"

allocate( alpha(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:alpha"
allocate( beta(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:beta"
allocate( respg(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:respg"
allocate( respb(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:respb"
allocate( QminN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QminN"
allocate( QminP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QminP"
allocate( QmaxN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QmaxN"
allocate( QmaxP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QmaxP"
allocate( Kn(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Kn"
allocate( Kp(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Kp"
allocate( Ksi(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Ksi"
allocate( KQn(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:KQn"
allocate( KQp(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:KQp"
allocate( nfQs(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:nfQs"
allocate( vmaxN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxN"
allocate( vmaxP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxP"
allocate( vmaxSi(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxSi"
allocate( aN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:aN"
allocate( volcell(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:volcell"
allocate( Qc(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Qc"
allocate( Athresh(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Athresh"
allocate( mA(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:mA"
allocate( A_wt(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:A_wt"

!---Zooplankton
allocate( Zeffic(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zeffic"
allocate( Zslop(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zslop"
allocate( Zvolcell(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zvolcell"
allocate( ZQc(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQc"
allocate( ZQn(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQn"
allocate( ZQp(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQp"
allocate( ZKa(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZKa"
allocate( Zrespg(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zrespg"
allocate( Zrespb(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zrespb"
allocate( Zumax(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zumax"
allocate( Zm(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zm"

!Light curve parameters
allocate( alphad(nospA),stat=ierr ) ! Initial slope of photosynthesis-irradiance curve / Vmax
if(ierr.ne.0) write(6,*) "error in allocating:alphad"
allocate( betad(nospA),stat=ierr )  ! Photoinhibition constant / Vmax
if(ierr.ne.0) write(6,*) "error in allocating:betad"

!Diatom array
allocate( is_diatom(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:is_diatom"

!Temperature parameters for growth rates
allocate(Tref(nospA+nospZ),stat=ierr)                   !Tref(nospA+nospZ): Optimum temperature for growth(C)
if(ierr.ne.0) write(6,*) "error in allocating:Tref"
allocate(KTg1(nospA+nospZ),stat=ierr)                   !KTg1(nospA+nospZ): Effect of T below Topt(C^2)
if(ierr.ne.0) write(6,*) "error in allocating:KTg1"
allocate(KTg2(nospA+nospZ),stat=ierr)                   !KTg2(nospA+nospZ): Effect of T above Topt(C^2)
if(ierr.ne.0) write(6,*) "error in allocating:KTg2"
allocate(Ea(nospA+nospZ),stat=ierr)                     !Ea(nospA+nospZ): Slope of Arrhenius plot(eV)
if(ierr.ne.0) write(6,*) "error in allocating:Ea"


return
end subroutine cgem_vars_allocate

subroutine cgem_init 

integer                          :: fu, rc

real sinkA,sinkOM1_A,sinkOM2_A,sinkOM1_Z,sinkOM2_Z,sinkOM1_R,sinkOM2_R,sinkOM1_BC,sinkOM2_BC

namelist /switches/ Which_fluxes,Which_temperature,Which_uptake,Which_quota,Which_irradiance,Which_chlaC,Which_photosynthesis,Which_growth
namelist /optics/ Kw,Kcdom,Kspm,Kchla,astar490,aw490,astarOMA,astarOMZ,astarOMR,astarOMBC,PARfac
namelist /temperature/ Tref,KTg1,KTg2,Ea
namelist /phytoplankton/ ediblevector,umax,CChla,alpha,beta,respg,respb,&
 QminN,QminP,QmaxN,QmaxP,Kn,Kp,Ksi,KQn,KQp,nfQs,vmaxN,vmaxP,vmaxSi,aN,volcell,Qc,Athresh,sinkA,mA,A_wt
namelist /zooplankton/ Zeffic,Zslop,Zvolcell,ZQc,ZQn,ZQp,ZKa,Zrespg,Zrespb,Zumax,Zm
namelist /OM/ KG1,KG2,KG1_R,KG2_R,KG1_BC,KG2_BC,KNH4,nitmax,KO2,KstarO2,KNO3,pCO2,&
 stoich_x1R,stoich_y1R,stoich_x2R,stoich_y2R,stoich_x1BC,stoich_y1BC,stoich_x2BC,stoich_y2BC,&
 sinkOM1_A,sinkOM2_A,sinkOM1_Z,sinkOM2_Z,sinkOM1_R,sinkOM2_R,sinkOM1_BC,sinkOM2_BC,KGcdom,CF_SPM,KG_bot

!Later, look at this python namelist thing:
!https://github.com/marshallward/f90nml
!install
!conda create --prefix ./env_f90 -c conda-forge f90nml
!(conda activate it)
!conda install -c conda-forge ascii_graph
open(action='read',file='cgem.nml',iostat=rc,newunit=fu)

!namelist /switches/
read(nml=switches,iostat=rc,unit=fu)
if (rc /= 0) write (stderr, '("Error: invalid Namelist format, switches")')

!namelist /optics/
read(nml=optics,iostat=rc,unit=fu)
if (rc /= 0) write (stderr, '("Error: invalid Namelist format, optics")')

!namelist /temperature/
read(nml=temperature,iostat=rc,unit=fu)
if (rc /= 0) write (stderr, '("Error: invalid Namelist format, temperature")')
close(fu)

!namelist /phytoplankton/
read(nml=phytoplankton,iostat=rc,unit=fu)
if (rc /= 0) write (stderr, '("Error: invalid Namelist format, phytoplankton")')

!namelist /zooplankton/
read(nml=zooplankton,iostat=rc,unit=fu)
if (rc /= 0) write (stderr, '("Error: invalid Namelist format, zooplankton")')

!namelist /OM/
read(nml=OM,iostat=rc,unit=fu)
if (rc /= 0) write (stderr, '("Error: invalid Namelist format, OM")')


end subroutine cgem_init



end module cgem_vars

