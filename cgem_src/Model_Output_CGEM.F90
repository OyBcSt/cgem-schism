      Subroutine Model_Output_CGEM

      use grid_vars
      use cgem_vars
      use OUTPUT_NETCDF_CGEM

      IMPLICIT NONE

      real :: dumf(1,1,km,nf)
      integer :: k,im,jm
      im = 1
      jm = 1

          do k=1,km
            dumf(1,1,k,:) = ff(k,:)
          enddo

        CALL WRITE_DATA( 1,im, 1,jm, 1, km, istep_out, dumf)

      return

      End Subroutine
