module cgem_functions

contains

  FUNCTION Fixed_CChla(A_k,nz) RESULT(Chla_tot)

    use Model_dim
    use INPUT_VARS_CGEM, ONLY:Qc,CChla
    implicit none

    ! Input parameters
    real, intent(in) :: A_k(nospA,km)  ! A's number density, cells/m3
    integer, intent(in) :: nz   ! Number of layers
    ! Function return value
    real :: Chla_tot(km)

    ! Local variables
    integer :: k, isp

    DO k = 1, nz
       Chla_tot(k) = 0.0
       DO isp = 1, nospA
          Chla_tot(k) =  Chla_tot(k) + A_k(isp,k) * Qc(isp) * 12. *
          C(1./CChla(isp))
       ENDDO ! isp = 1, nospA
    ENDDO ! k = 1, nz
    RETURN
  END FUNCTION Fixed_CChla

end module
