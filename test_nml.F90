program test_nml
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
implicit none

real :: Kw,Kcdom,Kspm,Kchla,astar490,aw490,astarOMA,astarOMZ,astarOMR,astarOMBC,PARfac
integer                          :: fu, rc
namelist /optics/ Kw,Kcdom,Kspm,Kchla,astar490,aw490,astarOMA,astarOMZ,astarOMR,astarOMBC,PARfac

!Kw = 0.146            !Kw: AOP, light attenuation due to water
!Kcdom = 0.001         !Kcdom: AOP, light attenuation due to CDOM
!Kspm = 0.029          !Kspm: AOP, light attenuation due to SPM
!Kchla = 0.024         !Kchla: AOP, light attenuation due to chla
!astar490 = 0.0375     !astar490: Chla specific absorption at 490 nm
!aw490 = 0.015         !aw490: seawater absorption at 490 nm
!astarOMA = 0.01       !astarOMA: OM_A specific absorption at 490 nm
!astarOMZ = 0.01       !astarOMZ: OM_Z specific absorption at 490 nm
!astarOMR = 0.01       !astarOMR: OM_R specific absorption at 490 nm
!astarOMBC = 0.01      !astarOMBC: OM_BC specific absorption at 490 nm
!PARfac = 1.           !PARfac: Multiplies surface PAR

open(action='read',file='cgem.nml',iostat=rc,newunit=fu)
read(nml=optics,iostat=rc,unit=fu)
if (rc /= 0) write (stderr, '("Error: invalid Namelist format")')
close(fu)

write(6,*) "optics"
write(6,*) "Kw,Kcdom,Kspm,Kchla,astar490,aw490,astarOMA,astarOMZ,astarOMR,astarOMBC,PARfac"
write(6,*) Kw,Kcdom,Kspm,Kchla,astar490,aw490,astarOMA,astarOMZ,astarOMR,astarOMBC,PARfac

end program test_nml
