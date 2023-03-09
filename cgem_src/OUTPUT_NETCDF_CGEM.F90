!******************************************************************************
! PURPOSE: OUTPUT_NETCDF.F90 - Routines for outputting NCOM_GEM data to a
!          NetCDF file.
! NOTES:   Non-ADT module.
! HISTORY: 2010/04/26, Todd Plessel, plessel.todd@epa.gov, Created.
!******************************************************************************

MODULE OUTPUT_NETCDF_CGEM

  USE xnetcdf
  USE NETCDF_UTILITIES ! For CHKERR, DEFDIM, DEFVI1, CONVERT_LONGITUDES, etc.
  USE DATE_TIME ! For TOTAL_SECONDS
  USE OUTPUT
  USE cgem, only:nospA,nospZ,nf

  IMPLICIT NONE

  ! Private

  INTEGER TIME_VAR ! NetCDF ID for time array variable.
  INTEGER,SAVE :: FILE_ID ! NetCDF ID for file.
  INTEGER FILE_FIRST_TIMESTEP ! 0-based model timestep number of current file.
  INTEGER SECONDS_PER_TIMESTEP ! Output timestep size in seconds.
  INTEGER(8) SECONDS0 ! Seconds since Model_dim::iYr0.

PUBLIC CREATE_FILE, OPEN_FILE, &
       WRITE_DATA, CLOSE_FILE, FLUSH_FILE, &
       OUTPUT_NETCDF_CGEM_ALLOCATE

PRIVATE
CONTAINS

  ! Public

  ! Commands:

  ! CREATE_FILE: Create output NetCDF file with given header data.
  ! In a concurrent program, only one process should call this routine
  ! then call CLOSE_FILE then
  ! worker processes should call OPEN_FILE, WRITE_DATA, CLOSE_FILE.
  !
  SUBROUTINE CREATE_FILE( NAME, IM, JM, KM, NSTEP, IYR0, &
                          IYRS, IMONS, IDAYS, IHRS, IMINS, ISECS, &
                          IYRE, IMONE, IDAYE, IHRE, IMINE, ISECE, DT_OUT )
    USE OUTPUT 

    IMPLICIT NONE
    include 'netcdf.inc'
    CHARACTER(LEN=*),INTENT(IN):: NAME
    INTEGER,INTENT(IN):: IM, JM, KM
    INTEGER,INTENT(IN):: NSTEP
    INTEGER,INTENT(IN):: IYR0 ! Reference year (before start of model run).
    INTEGER,INTENT(IN):: IYRS, IMONS, IDAYS, IHRS, IMINS, ISECS ! Run start.
    INTEGER,INTENT(IN):: IYRE, IMONE, IDAYE, IHRE, IMINE, ISECE ! Run end.
    INTEGER,INTENT(IN):: DT_OUT ! Model timestep size in seconds.
   ! Locals:
    INTEGER IM_DIM, JM_DIM, KM_DIM, NSTEPP1_DIM
    INTEGER K, i, J, INFO
    INTEGER ERR, VARIABLE, DIM_IDS( 4 )
!    INTEGER NF_CLOBBER, NF_64BIT_OFFSET
!    EXTERNAL NF_CLOBBER, NF_64BIT_OFFSET
    REAL,DIMENSION(IM,JM):: RLON_COPY
    CHARACTER(LEN=40):: TIME_UNITS
    CHARACTER(LEN=14):: var
    INTEGER(KIND=MPI_OFFSET_KIND) :: nospAZ, nospA_m, nospZ_m
    FILE_ID = -1
    nospAZ=nospA+nospZ
    nospA_m = nospA
    nospZ_m = nospZ

    CALL MPI_INFO_CREATE( INFO, ERR )
    CALL MPI_INFO_SET( INFO, 'ind_wr_buffer_size', '16777216', ERR )
    ! Create/overwrite NetCDF output file:
    ERR = ncdf_CREATE( MPI_COMM_SELF,trim(NAME), IOR( NF_CLOBBER, NF_64BIT_OFFSET ), INFO, FILE_ID) 
    CALL CHKERR( ERR, 'create NetCDF output file ' // NAME )
    CALL MPI_INFO_FREE( INFO, ERR )

    ! Create dimensions:
    CALL DEFDIM( FILE_ID, IM_DIM, 'longitude', IM )
    CALL DEFDIM( FILE_ID, JM_DIM, 'latitude', JM )
    CALL DEFDIM( FILE_ID, KM_DIM, 'k', KM )
!!  CALL DEFDIM( FILE_ID, NSTEPP1_DIM, 'time', NSTEP )
    CALL DEFDIM( FILE_ID, NSTEPP1_DIM, 'time', 0 ) ! 0 Means UNLIMITED size.
    ! Write global scalar attributes:

    CALL DEFIAT( FILE_ID, 'nstep', NSTEP )
    CALL DEFIAT( FILE_ID, 'iYr0', IYR0 )
    CALL DEFIAT( FILE_ID, 'iYrS', IYRS )
    CALL DEFIAT( FILE_ID, 'iMonS', IMONS )
    CALL DEFIAT( FILE_ID, 'iDayS', IDAYS )
    CALL DEFIAT( FILE_ID, 'iHrS', IHRS )
    CALL DEFIAT( FILE_ID, 'iMinS', IMINS )
    CALL DEFIAT( FILE_ID, 'iSecS', ISECS )
    CALL DEFIAT( FILE_ID, 'iYrE', IYRE )
    CALL DEFIAT( FILE_ID, 'iMonE', IMONE )
    CALL DEFIAT( FILE_ID, 'iDayE', IDAYE )
    CALL DEFIAT( FILE_ID, 'iHrE', IHRE )
    CALL DEFIAT( FILE_ID, 'iMinE', IMINE )
    CALL DEFIAT( FILE_ID, 'iSecE', ISECE )
    CALL DEFIAT( FILE_ID, 'dT_out', DT_OUT )


    ! Define time array variable as each output data's seconds since IYR0:
    WRITE ( TIME_UNITS, '(A,I4.4,A)' ) &
            'Seconds since ', iYr0, '-01-01 00:00:00Z'

    ! Note: this is defined as an 8-byte real in the file
    ! because NetCDF does not support 8-byte integers.
    ! Its value is converted to/from interger during reads/writes.

    CALL DEFVD1( FILE_ID, NSTEPP1_DIM, TIME_VAR, 'time', &
                 TRIM( TIME_UNITS ), TRIM( TIME_UNITS ) )

    ! Define time-varying array variables:

    DIM_IDS( 1 ) = IM_DIM
    DIM_IDS( 2 ) = JM_DIM
    DIM_IDS( 3 ) = KM_DIM
    DIM_IDS( 4 ) = NSTEPP1_DIM

    DO VARIABLE = 1, STATE_VARIABLES
      IF ( WRITE_VARIABLE( VARIABLE ) ) THEN
        CALL DEFVR4( FILE_ID, DIM_IDS, F_VAR( VARIABLE ), &
                     TRIM( VARIABLE_NAMES( VARIABLE ) ), &
                     TRIM( VARIABLE_DESCRIPTIONS( VARIABLE ) ), &
                     TRIM( VARIABLE_UNITS( VARIABLE ) ),1 )
      END IF
    END DO

    ERR = ncdf_ENDDEF( FILE_ID )
    CALL CHKERR( ERR, 'create NetCDF output header' )
    ERR = ncdf_BEGIN_INDEP_DATA( FILE_ID )
    CALL CHKERR( ERR, 'begin independent data access mode' )


    CALL FLUSH_FILE()


    RETURN
  END SUBROUTINE CREATE_FILE



  ! OPEN_FILE: Open existing output NetCDF file for shared writing.
  ! In a concurrent program, each worker process should call
  ! OPEN_FILE, WRITE_DATA, CLOSE_FILE.
  ! FIRST_TIMESTEP is 0-based timestep number to begin writing.
  !
  SUBROUTINE OPEN_FILE( NAME, FIRST_TIMESTEP )
    USE OUTPUT 

    IMPLICIT NONE

    include 'netcdf.inc'

    CHARACTER(LEN=*),INTENT(IN):: NAME
    INTEGER,INTENT(IN):: FIRST_TIMESTEP
    ! Locals:
    INTEGER ERR, VARIABLE
    INTEGER DT_OUT
    INTEGER IYR0 ! Reference year (before start of model run).
    INTEGER IYRS, IMONS, IDAYS, IHRS, IMINS, ISECS ! Run start.
    INTEGER(8) SECONDS_FROM_YEAR0
    INTEGER BUFFER_SIZE,INFO

    CALL MPI_INFO_CREATE( INFO, ERR )
    CALL MPI_INFO_SET( INFO, 'ind_wr_buffer_size', '16777216', ERR )

    ! Open existing shared 64-bit NetCDF output file for writing:
!    ERR = ncdf_OPEN( NAME, 2565, FILE_ID )
    ERR = ncdf_OPEN( MPI_COMM_WORLD, trim(NAME), &
                      IOR( NF_WRITE, NF_64BIT_OFFSET ), INFO, FILE_ID )
    !BUFFER_SIZE = 256 * 1024
    !ERR = ncdf_OPEN( NAME, NUM, BUFFER_SIZE, FILE_ID ) ! 1+4+256+2+512 64-bit
    CALL CHKERR( ERR, 'open existing shared writable output file ' // NAME )
    CALL MPI_INFO_FREE( INFO, ERR )

    ! Get time variable id:

    ERR = ncdf_INQ_VARID( FILE_ID, 'time', TIME_VAR )
    CALL CHKERR( ERR, 'inquire NetCDF variable ID ' )

    ! Read time attributes to compute SECONDS0:
    
    DT_OUT = READIAT( FILE_ID, 'dT_out' )
    IYR0   = READIAT( FILE_ID, 'iYr0' )
    IYRS   = READIAT( FILE_ID, 'iYrS' )
    IMONS  = READIAT( FILE_ID, 'iMonS' )
    IDAYS  = READIAT( FILE_ID, 'iDayS' )
    IHRS   = READIAT( FILE_ID, 'iHrS' )
    IMINS  = READIAT( FILE_ID, 'iMinS' )
    ISECS  = READIAT( FILE_ID, 'iSecS' )
    SECONDS_FROM_YEAR0 = &
      TOTAL_SECONDS( IYR0, IYRS, IMONS, IDAYS, IHRS, IMINS, ISECS )

    FILE_FIRST_TIMESTEP = FIRST_TIMESTEP ! 0-based model timestep of this file.
    SECONDS_PER_TIMESTEP = DT_OUT
    SECONDS0 = SECONDS_PER_TIMESTEP
    SECONDS0 = SECONDS0 * FILE_FIRST_TIMESTEP
    SECONDS0 = SECONDS0 + SECONDS_FROM_YEAR0

    ! Get time-varying array variable ids:

    DO VARIABLE = 1, STATE_VARIABLES 

      IF ( WRITE_VARIABLE( VARIABLE ) ) THEN
        ERR = ncdf_INQ_VARID( FILE_ID, TRIM(VARIABLE_NAMES( VARIABLE )), &
                               F_VAR( VARIABLE ) )
        CALL CHKERR( ERR, 'inquire NetCDF variable ID ' )
      END IF
    END DO

    RETURN
  END SUBROUTINE OPEN_FILE



  ! CLOSE_FILE: Close output file.
  !
  SUBROUTINE CLOSE_FILE()
    IMPLICIT NONE
    INTEGER ERR

    ERR = ncdf_CLOSE( FILE_ID )
    CALL CHKERR( ERR, 'close NetCDF output file ' )
    FILE_ID = -1

    RETURN
  END SUBROUTINE CLOSE_FILE


  ! WRITE_DATA: Write current timestep data to output file.
  ! Called by non-concurrent programs to write all data per timestep.
  !

  SUBROUTINE WRITE_DATA( IMSTART, IM, JMSTART, JM, KMSTART, KM, TIMESTEP, F ) 
    USE OUTPUT
    IMPLICIT NONE
    INTEGER,INTENT(IN):: IMSTART, JMSTART, KMSTART
    INTEGER,INTENT(IN):: IM, JM, KM, TIMESTEP ! Model TIMESTEP is 0-based.
    REAL,DIMENSION(IM, JM, KM, STATE_VARIABLES):: F
    ! Locals:
    REAL(8):: SECONDS(1) ! TIME_VAR.
    INTEGER ERR, VARIABLE, FILE_TIMESTEP
    INTEGER(KIND=MPI_OFFSET_KIND) STARTS(4),COUNTS(4)
    INTEGER(KIND=MPI_OFFSET_KIND) STARTS1(1),COUNTS1(1)
    INTEGER REQUEST_COUNT,REQUEST
    INTEGER,DIMENSION(STATE_VARIABLES+1):: REQUESTS, STATUSES

    FILE_TIMESTEP = TIMESTEP - FILE_FIRST_TIMESTEP

    ! Write time variable as an 8-byte real since NetCDF lacks 8-byte integer:

    STARTS1( 1 ) = FILE_TIMESTEP + 1
    COUNTS1( 1 ) = 1
    SECONDS( 1 ) = SECONDS0 + FILE_TIMESTEP * SECONDS_PER_TIMESTEP
    !write(6,*) "seconds",seconds
    ERR = ncdf_PUT_VARA_DOUBLE( FILE_ID, TIME_VAR, STARTS1, COUNTS1, SECONDS, REQUESTS(1)) 
    CALL CHKERR( ERR, 'write output variable time' )
    !write(6,*) "double",FILE_ID, TIME_VAR, STARTS1, COUNTS1, SECONDS, REQUESTS(1)    

    REQUEST_COUNT = 1

    STARTS( 1 ) = IMSTART  
    STARTS( 2 ) = JMSTART 
    STARTS( 3 ) = KMSTART
    STARTS( 4 ) = FILE_TIMESTEP + 1 ! NetCDF follows FORTRAN 1-based convention
    COUNTS( 1 ) = IM
    COUNTS( 2 ) = JM
    COUNTS( 3 ) = KM
    COUNTS( 4 ) = 1


    DO VARIABLE = 1, STATE_VARIABLES 
        REQUEST_COUNT = REQUEST_COUNT + 1
        ERR = ncdf_PUT_VARA_REAL( FILE_ID, F_VAR( VARIABLE ), &
              STARTS, COUNTS, F( 1, 1, 1, VARIABLE ), REQUESTS( REQUEST_COUNT ))
        CALL CHKERR( ERR, 'write output variable ' // VARIABLE_NAMES(VARIABLE))
    END DO

    ERR = ncdf_WAIT_ALL( FILE_ID, REQUEST_COUNT, REQUESTS, STATUSES )
    CALL CHKERR( ERR, 'implement non-blocking interface' )

    DO REQUEST = 1, REQUEST_COUNT
      CALL CHKERR( STATUSES( REQUEST ), 'nonblocking call ' )
    END DO


    ! check status of each nonblocking call

    CALL FLUSH_FILE() ! Flush buffers to disk in case of crash.

    RETURN
  END SUBROUTINE WRITE_DATA

  ! FLUSH_FILE: Flush buffers to disk in case of crash.
  !
  SUBROUTINE FLUSH_FILE()
    IMPLICIT NONE
    ! External NetCDF routines:
    INTEGER ERR
    ERR = ncdf_SYNC( FILE_ID )
    CALL CHKERR( ERR, 'flush buffers to disk ' )

    RETURN
  END SUBROUTINE FLUSH_FILE


  ! Private



!******************************************************************************

Subroutine OUTPUT_NETCDF_CGEM_allocate

  USE OUTPUT 

  IMPLICIT NONE

  integer i,counter
  character(6) var
  character(19) var1
  character(100) var2

  ALLOCATE(VARIABLE_NAMES(nf))
   counter = 0
 if(nospA.le.9) then
  do i=1,nospA
    counter = counter + 1
    write(var,'(A1,i1)') 'A',i
    VARIABLE_NAMES(counter) = var
  enddo 
  do i=1,nospA
    counter = counter + 1
    write(var,'(A2,i1)') 'Qn',i
    VARIABLE_NAMES(counter) = var
  enddo
  do i=1,nospA
    counter = counter + 1
    write(var,'(A2,i1)') 'Qp',i
    VARIABLE_NAMES(counter) = var
  enddo
 elseif(nospA.le.99) then
  do i=1,nospA
    counter = counter + 1
    write(var,'(A1,i2)') 'A',i
    VARIABLE_NAMES(counter) = var
  enddo
  do i=1,nospA
    counter = counter + 1
    write(var,'(A2,i2)') 'Qn',i
    VARIABLE_NAMES(counter) = var
  enddo
  do i=1,nospA
    counter = counter + 1
    write(var,'(A2,i2)') 'Qp',i
    VARIABLE_NAMES(counter) = var
  enddo
 else
  write(6,*) "No more than 99 phytoplankton species allowed in netCDF"
  stop
 endif
 if(nospZ.le.9) then
  do i=1,nospZ
    counter = counter + 1
    write(var,'(A1,i1)') 'Z',i
    VARIABLE_NAMES(counter) = var
  enddo
 elseif(nospZ.le.99) then
  do i=1,nospZ
    counter = counter + 1
    write(var,'(A1,i2)') 'Z',i
    VARIABLE_NAMES(counter) = var
  enddo
 else
  write(6,*) "No more than 99 phytoplankton species allowed in netCDF"
  stop
 endif
    VARIABLE_NAMES(counter+1) = 'NO3   '
    VARIABLE_NAMES(counter+2) = 'NH4   '
    VARIABLE_NAMES(counter+3) = 'PO4   '
    VARIABLE_NAMES(counter+4) = 'DIC   '
    VARIABLE_NAMES(counter+5) = 'O2    '
    VARIABLE_NAMES(counter+6) = 'OM1_A '
    VARIABLE_NAMES(counter+7) = 'OM2_A '
    VARIABLE_NAMES(counter+8) = 'OM1_Z '
    VARIABLE_NAMES(counter+9) = 'OM2_Z '
    VARIABLE_NAMES(counter+10) = 'OM1_R '
    VARIABLE_NAMES(counter+11) = 'OM2_R '
    VARIABLE_NAMES(counter+12) = 'CDOM  '
    VARIABLE_NAMES(counter+13) = 'Si    '
    VARIABLE_NAMES(counter+14) = 'OM1_BC'
    VARIABLE_NAMES(counter+15) = 'OM2_BC'
    VARIABLE_NAMES(counter+16) = 'ALK   '
    VARIABLE_NAMES(counter+17) = 'Tr   '


  ALLOCATE( WRITE_VARIABLE(nf) ) 
  do i=1,nf
   WRITE_VARIABLE(i) = .TRUE.
!   WRITE_VARIABLE(i) = .FALSE.
  enddo


  ALLOCATE(VARIABLE_DESCRIPTIONS(nf))
   counter = 0
  do i=1,nospA
    counter = counter + 1
    write(var2,'(A21,i3,A15)') 'Phytoplankton group ',i," number density"
    VARIABLE_DESCRIPTIONS(counter) = var2
  enddo
  do i=1,nospA
    counter = counter + 1
    write(var2,'(A22,i3,A20)') 'Phytoplankton group ',i," nitrogen quota."
    VARIABLE_DESCRIPTIONS(counter) = var2
  enddo
  do i=1,nospA
    counter = counter + 1
    write(var2,'(A22,i3,A20)') 'Phytoplankton group ',i," phosphorus quota." 
    VARIABLE_DESCRIPTIONS(counter) = var2
  enddo
   do i=1,nospZ
    counter = counter + 1
    write(var2,'(A22,i3,A20)') 'Zooplankton group ',i," number density."
    VARIABLE_DESCRIPTIONS(counter) = var2
  enddo
    VARIABLE_DESCRIPTIONS(counter+1) = 'Nitrate.  '
    VARIABLE_DESCRIPTIONS(counter+2) = 'Ammonium.                                                           '
    VARIABLE_DESCRIPTIONS(counter+3) = 'Phosphate.                                                          '
    VARIABLE_DESCRIPTIONS(counter+4) = 'Dissolved inorganic carbon.                                         '
    VARIABLE_DESCRIPTIONS(counter+5) = 'Molecular oxygen.                                                   '
    VARIABLE_DESCRIPTIONS(counter+6) = 'Particulate organic matter derived from dead algae.                 '
    VARIABLE_DESCRIPTIONS(counter+7) = 'Dissolved organic matter derived from dead algae.                   '
    VARIABLE_DESCRIPTIONS(counter+8) = 'Particulate organic matter derived from zooplankton fecal pellets.  '
    VARIABLE_DESCRIPTIONS(counter+9) = 'Dissolved organic matter derived from zooplankton fecal pellets.    '
    VARIABLE_DESCRIPTIONS(counter+10) = 'Particulate organic matter derived from river outflow.              '
    VARIABLE_DESCRIPTIONS(counter+11) = 'Dissolved organic matter derived from river outflow.                '
    VARIABLE_DESCRIPTIONS(counter+12) = 'Colored dissolved organic matter.                                   '
    VARIABLE_DESCRIPTIONS(counter+13) = 'Silica.                                                             '
    VARIABLE_DESCRIPTIONS(counter+14) = 'Particulate organic matter in the initial and boundary conditions '
    VARIABLE_DESCRIPTIONS(counter+15) = 'Dissolved organic matter in the initial and boundary conditions '
    VARIABLE_DESCRIPTIONS(counter+16) = 'Alkalinity.                                                         '
    VARIABLE_DESCRIPTIONS(counter+17) = 'Tracer should be =1.  '

  ALLOCATE(VARIABLE_UNITS(nf))
  counter = 0
  do i=1,nospA
    counter = counter + 1
    write(var2,'(A8)') "cells/m3" 
    VARIABLE_UNITS(counter) = var2
  enddo
  do i=1,nospA
    counter = counter + 1
    write(var2,'(A11)') "mmol-N/cell"
    VARIABLE_UNITS(counter) = var2
  enddo
  do i=1,nospA
    counter = counter + 1
    write(var2,'(A11)') "mmol-P/cell"
    VARIABLE_UNITS(counter) = var2
  enddo
  do i=1,nospZ
    counter = counter + 1
    write(var2,'(A12)') "organisms/m3"
    VARIABLE_UNITS(counter) = var2
  enddo
    VARIABLE_UNITS(counter+1) = 'mmol-N/m3                       '
    VARIABLE_UNITS(counter+2) = 'mmol-N/m3                       '
    VARIABLE_UNITS(counter+3) = 'mmol-P/m3                       '
    VARIABLE_UNITS(counter+4) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+5) = 'mmol-O2/m3                      '
    VARIABLE_UNITS(counter+6) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+7) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+8) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+9) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+10) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+11) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+12) = 'ppb                             '
    VARIABLE_UNITS(counter+13) = 'mmol-Si/m3                      '
    VARIABLE_UNITS(counter+14) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+15) = 'mmol-C/m3                       '
    VARIABLE_UNITS(counter+16) = 'mmol/m3                         '
    VARIABLE_UNITS(counter+17) = 'NONE                         '

  ALLOCATE(F_VAR(nf)) ! NetCDF IDs for each variable.
  F_VAR = -9999. 

  STATE_VARIABLES = nf 

  return

END Subroutine OUTPUT_NETCDF_CGEM_allocate

END MODULE OUTPUT_NETCDF_CGEM

