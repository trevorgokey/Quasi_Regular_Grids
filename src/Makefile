#==============================================================================#
#			   Fortran Compiler/Flags
#-O0=no optimization,-O3=full optimization
#==============================================================================#
FC=gfortran 
FFLAGS=-O0 -gdwarf  -Wall#-O3 -march=native -mtune=native -fomit-frame-pointer -pipe -fopenmp
LDFLAGS= #-Wl,-O2 -Wl,--as-needed
FFLAGS=-O3 -march=native -mtune=native -fomit-frame-pointer -pipe -fopenmp
LDFLAGS=-Wl,-O2 -Wl,--as-needed

#OBJS=grid_LJ.o grid.o system.o 
#SRC=grid_LJ.f90 grid.f90 system.f90 system_3d.f90 

#OBJS=system_types.o grid_types.o space.o  qrg.o
#SRC=system_types.f90 grid_types.f90 system.f90 system_basis.f90 system_basis_gaussian.f90 space.f90 grid.f90 grid_LJ.f90 qrg.f90

GRID_OBJS=grid_types.o grid.o grid_LJ.o space.o
SYSTEM_OBJS=system_types.o system.o system_basis.o system_basis_gaussian.o space.o 

SYSTEMS=systems.a
GRIDS=grids.a
MYLIBS=$(SYSTEMS) $(GRIDS)
LIBS=-lgomp

qrg: qrg.o
	$(FC) $(LDFLAGS) -o qrg qrg.o -L. $(MYLIBS) $(LIBS)

qrg.o: qrg.f90 $(MYLIBS)
	$(FC) $(FFLAGS) -c qrg.f90

grid_LJ.o: grid_LJ.f90 grid.o
	$(FC) $(FFLAGS) -c grid_LJ.f90

grid.o: grid.f90 system.o space.o
	$(FC) $(FFLAGS) -c grid.f90


system.o: system.f90 space.o
	$(FC) $(FFLAGS) -c system.f90

system_basis_gaussian.o: system_basis_gaussian.f90 system_basis.o
	$(FC) $(FFLAGS) -c system_basis_gaussian.f90

system_basis.o: system_basis.f90 system.o
	$(FC) $(FFLAGS) -c system_basis.f90

space.o: space.f90
	$(FC) $(FFLAGS) -c space.f90

grid_types.o: grid_types.f90 grid_LJ.o
	$(FC) $(FFLAGS) -c grid_types.f90

system_types.o: system_types.f90 system_basis_gaussian.o
	$(FC) $(FFLAGS) -c system_types.f90

systems.a: $(SYSTEM_OBJS)
	ar ru systems.a $(SYSTEM_OBJS)

grids.a: $(GRID_OBJS)
	ar ru grids.a $(GRID_OBJS)

.PHONY: clean

clean:
	rm -rf *.o *.mod *.dat qrg *.a
