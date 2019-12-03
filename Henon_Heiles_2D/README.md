# Program Main:
QMC Implementation 

The user defines the number of Basis Functions to describe the monomer with (Must be less than 10), the highest excitation available to each indivdual basis function (must be less than 10), and the total excitation for the entire system. 
The program then computes the Hessian and eigenvectors for the water cluser input geometry. 
A normal mode transformation is then computed and the Potential Difference Matrix (Potential Energy Surface - Harmonic Approximation) is constructed.
This matrix is then diagonalized (using the LAPACK library) and the fundamental frequencies for the first water within the cluser are computed.

See the development directory for an implementation computing all of the monomers fundamental frequencies within the cluster. 
This implementation should be verified and uploaded as the main (4-29-18)!

## Software Overview:
A quick summary of the various pieces of code available. 

#### Makefile
A simple Makefile for compiling the program.

#### sobol_stdnormal.f90:
This fortran program takes in a vector of sobol points (0,1) and uses the Beasley-Springer-Moro algorithm to transform it to a Gaussian Distribution. 

#### run_code
This directory contains an example input file and input geometry for running a simulation.
