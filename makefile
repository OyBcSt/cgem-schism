#
## Makefile to build the CGEM FishTank model
#
## Note:
# - These settings are for the *INTEL* compiler.
#   # - Ensure you load the appropriates modules before building (if applicable).

### =============== User Modifiable Section =============== ###

### Uncomment the next line to enable debugging
#DFLAGS = -DDEBUG
#DFLAGS = -g -warn -debug all -g -check all -ftrapuv#  -DDEBUG #-mcmodel=medium -shared-intel 
#DFLAGS = -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -D_CGEM -DRDEBUG -DDEBUG 
#DFLAGS = -g

### Build options for specific platforms. 
### LIBS has path to netCDF
#INC   = -I. -I/usr/local/apps/pnetcdf-1.9.0/intel-18.0/include -I/usr/local/apps/netcdf-4.6.3/intel-18.0/include/
#LIBS  = -L. -L/usr/local/apps/pnetcdf-1.9.0/intel-18.0/lib -lpnetcdf -L/usr/local/apps/netcdf-4.6.3/intel-18.0/lib -lnetcdff -lnetcdf

### =============== End User Modifiable Section  =============== ####

EXE = CGEM.exe

OBJ  =  \
	cgem_vars.o\
        sgrid.o\
        cgem.o\
        main.o

cgem: ${OBJ} 
	$(F90) -o $(EXE) $(FFLAGS) $(DFLAGS) $(INC) $(OBJ) $(LIBS)


#
## Pattern rules
#
$(OBJ):%.o: %.F90
	$(F90) -c $(FFLAGS) $(DFLAGS) $<

clean:
	rm -f *.o ${EXE} *.mod *genmod*

tags:
	ctags --language-force=Fortran *.F90

etags:
	etags -l fortran *.F90

