       Subroutine Model_Finalize_CGEM

       USE OUTPUT_NETCDF_CGEM

       IMPLICIT NONE


!Just close the netCDF file
         CALL CLOSE_FILE()

       return

       End Subroutine Model_Finalize_CGEM
