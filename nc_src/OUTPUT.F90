       Module OUTPUT 

        SAVE

        CHARACTER(LEN=6),DIMENSION(:),allocatable :: VARIABLE_NAMES
        LOGICAL,DIMENSION(:), allocatable :: WRITE_VARIABLE
        CHARACTER(LEN=100),DIMENSION(:), allocatable :: VARIABLE_DESCRIPTIONS
        CHARACTER(LEN=100),DIMENSION(:), allocatable :: VARIABLE_UNITS
        INTEGER,DIMENSION(:),allocatable :: F_VAR ! NetCDF IDs for each variable.
        INTEGER,save :: STATE_VARIABLES


       END Module OUTPUT 
