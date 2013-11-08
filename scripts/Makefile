# Define include files

include ./configure.mpasatograds

# Define modules and associated libraries

MODS = kinds.o                                       \
       constants.o                                   \
       mpi-interface.o                               \
       namelist.o                                    \
       kdtree2.o                                     \
       sphere.o                                      \
       interpolation-interface.o                     \
       variable-interface.o                          \
       grads-interface.o                             \
       mpasatograds-interface.o

SUBS = ncep-routines.o                               \
       fishpak.o                                     \
       spherepack.o 

# Define compile time variables

NETCDFcomp = -L$(NETCDFPATH)/lib -lnetcdf -I$(NETCDFPATH)/include

# Define executables

all: clean MPASATOGRADS

MPASATOGRADS:                      Driver.f90  $(MODS) $(SUBS) 
	$(LD) $(LDFLAGS)           Driver.f90  $(MODS) $(SUBS) $(NETCDFcomp) -lm -o $@.exe

clean:
	-rm *.o *~ \#* *.mod *.bin fort.* *.exe >& /dev/null

# Define make rules

.f90.o:
	$(F90) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(OPTIMIZATION) $(DEBUG) -c $*.f90

.f.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(OPTIMIZATION) $(DEBUG) -c $*.f

.F.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(OPTIMIZATION) $(DEBUG) -c $*.F
