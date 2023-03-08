!-------------------------------------------------------------
  SUBROUTINE func_E( E, min_S, f_E,Which_photosynthesis,alphad,betad,nospA )   
!-------------------------------------------------------------

!  USE cgem_vars, ONLY:Which_photosynthesis,alphad,betad,nospA
    
  !--------------------------------------------------------------------------
  ! INPUT:  
  !   E = Irradiance at cell k (quanta/cm**2/sec)
  !   min_S = Minimum substrate value (of f_N, f_P, and Si)
  !
  ! OUTPUT:
  !   f_E = Dimensionless factor for light dependent phytoplankton growth rate 
  ! 
  ! REFERENCES:
  !   Which_photosynthesis  ==  1 : With photoinhibition, Platt et al. (1980)
  !                         ==  2 : Without photoinhibition
  !                         ==  3 : Nutrient dependent, Flynn (2003)
  !------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nospA,Which_photosynthesis
    REAL, INTENT(IN) :: alphad(nospA),betad(nospA)
    REAL, INTENT(IN) :: E    ! Irradiance (quanta/cm**2/sec)
    REAL, INTENT(IN),  DIMENSION(nospA) :: min_S ! Function of rate limiting nutrient
    REAL, INTENT(OUT),DIMENSION(nospA)  :: f_E   ! Growth rate factor (dimensionless) 
    REAL, parameter :: alpha = 1.93e-16


    if (Which_photosynthesis.eq.1) then         !With photoinhibition 
        f_E(1:nospA) = ( 1.0 - exp(-alphad(1:nospA) * E) ) * exp(-betad(1:nospA)*E)
    else if (Which_photosynthesis.eq.2) then    !Without photoinhibition
        f_E(1:nospA) = ( 1.0 - exp(-alphad(1:nospA) * E) )
    else if (Which_photosynthesis.eq.3) then    !Nutrient dependent
        f_E(1:nospA) = ( 1.0 - exp(-alphad(1:nospA) * E / min_S) )
    else if (Which_photosynthesis.eq.4) then !GoMDOM
        f_E(1:nospA) = tanh(alpha * E)
    else
        write(6,*) "Error in func_E"
        stop
    endif

     !write(6,*) "fE,alphad,E",f_E,alphad,E
 
    RETURN
  END SUBROUTINE func_E  
