subroutine sgrid

use cgem_vars

implicit none

!Allocates nea, km(nvrt), nospA/Z, ff
call cgem_vars_allocate
call cgem_init


!nvrt = km 
!idry_e = 0
!!idry_e(:) - wet/dry flag
!allocate(idry_e(nea))


end subroutine sgrid
