OBJS = Sobol S_std rej main
FC = gfortran -g -fcheck=all

all: Sobol rej main

Sobol:
	$(FC) -c sobol.f90

rej:
	$(FC) -c rejection.f90 

main:
	$(FC) -O rejection.o sobol.o 

clean:
	rm -rf *.o *.mod *.out $(OBJS)
