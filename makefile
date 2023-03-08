#
## Makefile to build the CGEM FishTank model
#
## Note:
# - These settings are for the *INTEL* compiler.
#   # - Ensure you load the appropriates modules before building (if applicable).

### =============== User Modifiable Section =============== ###

### Uncomment the next line to enable debugging
#DFLAGS = -DDEBUG
#DFLAGS = -g -warn -debug all -g -check all -ftrapuv -DDEBUG #-mcmodel=medium -shared-intel 
DFLAGS = -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -D_CGEM -DRDEBUG -DDEBUG 
#DFLAGS = -g

### Build options for specific platforms. 
### LIBS has path to netCDF
#INC   = -I. -I/usr/local/apps/pnetcdf-1.9.0/intel-18.0/include -I/usr/local/apps/netcdf-4.6.3/intel-18.0/include/
#LIBS  = -L. -L/usr/local/apps/pnetcdf-1.9.0/intel-18.0/lib -lpnetcdf -L/usr/local/apps/netcdf-4.6.3/intel-18.0/lib -lnetcdff -lnetcdf

### =============== End User Modifiable Section  =============== ####
include cgem_src/src_files
include moc_src/src_files
#include sdm_src/src_files
cgemdir=cgem_src
mocdir=moc_src
#sdmdir=sdm_src

F90 = gfortran

EXE = CGEM.exe

cgem: ${MOC_OBJ} ${OBJ} ${SDM_OBJ} 
	$(F90) -o $(EXE) $(FFLAGS) $(DFLAGS) $(INC) $(MOC_OBJ) $(OBJ) $(SDM_OBJ) $(LIBS)

#
## Pattern rules
#
$(OBJ):%.o: $(cgemdir)/%.F90
	$(F90) -c $(FFLAGS) $(DFLAGS) $<

$(SDM_OBJ): %.o: $(sdmdir)/%.f
	$(F90) -c $(FFLAGS_SDM) $<

$(MOC_OBJ):%.o: $(mocdir)/%.F90
	$(F90) -c $(FFLAGS) $(DFLAGS)  $<


clean:
	rm -f *.o ${EXE} *.mod *genmod*

tags:
	ctags --language-force=Fortran *.F90

etags:
	etags -l fortran *.F90

