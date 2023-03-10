       Subroutine Init_Output_CGEM(BASE_NETCDF_OUTPUT_FILE_NAME)

       USE OUTPUT_NETCDF_CGEM
       USE grid
       USE cgem
       USE date_time

       IMPLICIT NONE

       character(100),intent(in) :: BASE_NETCDF_OUTPUT_FILE_NAME
       character(256) :: NETCDF_OUTPUT_FILE_NAME
       real :: dumf(1,1,km,nf)
       integer :: k,im,jm
       im = 1
       jm = 1


      WRITE ( NETCDF_OUTPUT_FILE_NAME, '(A, I6.6, A)' )&
              trim(BASE_NETCDF_OUTPUT_FILE_NAME), 0, '.nc'

          CALL CREATE_FILE( trim(NETCDF_OUTPUT_FILE_NAME), &
                            im, jm, km, nstep, iYr0,       & 
                            IYRS, IMONS, IDAYS, IHRS, IMINS, ISECS, &
                            IYRE, IMONE, IDAYE, IHRE, IMINE, ISECE,dT_out )
          CALL CLOSE_FILE()


! Opens the output file for writing:
       CALL OPEN_FILE( trim(NETCDF_OUTPUT_FILE_NAME), 0 )


          do k=1,km
           dumf(1,1,k,:) = ff(k,:)
          enddo

        CALL WRITE_DATA( 1, im, 1,jm, 1, km, 0, dumf)

#ifdef DEBUG
write(6,*) "---- Init_Output_CGEM ---"
write(6,*) 
#endif

       return
       End Subroutine Init_Output_CGEM
