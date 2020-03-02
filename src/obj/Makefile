#==============================================================================#
#			   Fortran Compiler/Flags
#-O0=no optimization,-O3=full optimization
#==============================================================================#
FC=gfortran 
LIBS=

# DEBUG FLAGS
FFLAGS=-O0 -gdwarf  -Wall
LDFLAGS= 

# OPT FLAGS
FFLAGS=-O3 -march=native -mtune=native -fomit-frame-pointer -pipe
LDFLAGS=-Wl,-O2 -Wl,--as-needed

# USE PARALLELIZATION
FFLAGS+=-fopenmp
LIBS+=-lgomp

GRID_OBJS=grid_types.o grid.o grid_LJ.o space.o
SYSTEM_OBJS=system_types.o system.o system_basis.o system_basis_gaussian.o space.o 

SYSTEMS=systems.a
GRIDS=grids.a
MYLIBS=$(SYSTEMS) $(GRIDS)

SRCDIR=..
CORE=$(SRCDIR)/core
SYSTEM=$(SRCDIR)/system
GRID=$(SRCDIR)/grid

# MAIN BINARY

qrg: qrg.o
	$(FC) $(LDFLAGS) -o $(SRCDIR)/qrg qrg.o -L. $(MYLIBS) $(LIBS)

# CORE

qrg.o: $(CORE)/qrg.f90 $(MYLIBS)
	$(FC) $(FFLAGS) -c $(CORE)/qrg.f90 

systems.a: $(SYSTEM_OBJS)
	ar ru systems.a $(SYSTEM_OBJS)

grids.a: $(GRID_OBJS)
	ar ru grids.a $(GRID_OBJS)

system.o: $(CORE)/system.f90 space.o
	$(FC) $(FFLAGS) -c $(CORE)/system.f90

system_basis.o: $(CORE)/system_basis.f90 system.o
	$(FC) $(FFLAGS) -c $(CORE)/system_basis.f90

space.o: $(CORE)/space.f90
	$(FC) $(FFLAGS) -c $(CORE)/space.f90

grid.o: $(CORE)/grid.f90 system.o space.o
	$(FC) $(FFLAGS) -c $(CORE)/grid.f90

# SYSTEMS IMPLEMENTED

system_basis_gaussian.o: $(SYSTEM)/system_basis_gaussian.f90 system_basis.o
	$(FC) $(FFLAGS) -c $(SYSTEM)/system_basis_gaussian.f90

system_types.o: $(SYSTEM)/system_types.f90 system_basis_gaussian.o
	$(FC) $(FFLAGS) -c $(SYSTEM)/system_types.f90

# GRIDS IMPLEMENTED


grid_LJ.o: $(GRID)/grid_LJ.f90 grid.o
	$(FC) $(FFLAGS) -c $(GRID)/grid_LJ.f90

grid_types.o: $(GRID)/grid_types.f90 grid_LJ.o
	$(FC) $(FFLAGS) -c $(GRID)/grid_types.f90

# MISC
#
.PHONY: clean

clean:
	rm -rf *.o *.mod *.dat *.a
