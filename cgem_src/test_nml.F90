program test_nml
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
implicit none

integer                          :: istat, iunit 
real, dimension(3) :: umax,CChla,alpha,beta,respg,respb,QminN,QminP,QmaxN,QmaxP
real, dimension(3) :: Kn,Kp,Ksi,KQn,KQp,nfQs,vmaxN,vmaxP,vmaxSi,aN,volcell,Qc,Athresh,sinkA,mA,A_wt
namelist /phytoplankton/ umax,CChla,alpha,beta,respg,respb,&
 QminN,QminP,QmaxN,QmaxP,Kn,Kp,Ksi,KQn,KQp,nfQs,vmaxN,vmaxP,vmaxSi,aN,volcell,Qc,Athresh,sinkA,mA,A_wt
character(len=1000) :: line

open(action='read',file='cgem.nml',iostat=istat,newunit=iunit)
read(nml=phytoplankton,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist: '//trim(line)
endif
close(iunit)

write(6,*) "phytoplankton"
write(6,*) umax,CChla,alpha,beta,respg,respb,QminN,QminP,QmaxN,QmaxP
write(6,*) Kn,Kp,Ksi,KQn,KQp,nfQs,vmaxN,vmaxP,vmaxSi,aN,volcell,Qc,Athresh,sinkA,mA,A_wt
end program test_nml
