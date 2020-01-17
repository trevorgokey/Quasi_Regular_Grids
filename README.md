# Quasi-Regular Grids; Applications to the Vibrational Spectra Calculations
Fortran and Python implementations for generating a QRG using model Potential
Energy Functions.
See our JCP Communication [manuscript](https://doi.org/10.1063/1.5134677) for details.

Given a general distribution of interest the method generates a set of points
(grid) that optimally samples the distribution.
We demonstrate the utility of these grids by constructing a DGB
(Distributed Gaussian Basis) and computing the associated Vibrational
Eigenspectra.

## Potentials/Distributions Available
Individual scripts for studying a specific distribution of interest.

To test accuracy/utility we have various standard grid methods available for
each potential:
* Direct-Product
* Pseudo-Random with a cutoff
* Pseudo-Random with Rejection
* Quasi-Random (Sobol Sequence) with a cutoff
* Quasi-Random (Sobol Sequence) with Rejection
* Quasi-Regular

The following potentials are currently available:
* 2D Henon-Heiles Potential
* Morse Potential (sum of 1D Morse potentials)

## Analysis
Various scripts I have generated for studying/analyzing the accuracy of these
methods.

## miscellaneous
Various scripts I have used for understanding the project.

## External Codes
I have used various Fortran packages during this project.

#### Quasirandom Sequences:
The [Sobol Sequence Generator](https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html) was used for the Quasirandom Calculations (sobol.f90).
Special thanks to John Burkardt for general advice on quasi-random numbers and
their implementation.

#### Evaluation of the Potential Energy Matrix Elements:
We used Gauss-Hermite Quadrature to evaluate the Potential Energy Matrix
Elements.
Again we have used code made available by [John Burkardt](https://people.sc.fsu.edu/~jburkardt/f_src/gen_hermite_rule/gen_hermite_rule.html)
(gen_hermite_rule.f90).

Also the Monte Carlo and Quasi-Monte Carlo Wiki page ([MCQMC](http://roth.cs.kuleuven.be/wiki/Main_Page)) is a useful resource for generating and applying Monte Carlo methods.

## References/Context
Accurate RoVibrational Energy calculations have a long history in computational
chemistry/physics.
I found these references quite helpful while working on this project.

#### Distributed Gaussian Basis Sets
The project was motivated in part by the work of [Garashchuk and Light](https://aip.scitation.org/doi/abs/10.1063/1.1348022).
Special thanks to Sophya Garashcuk for discussing her previous results with us.

The work by [Poirier and Light](https://aip.scitation.org/doi/abs/10.1063/1.481787) is a good reference for using Distributed Gaussian Basis sets in the
context of RoVibrational spectroscopy.

Interested readers should also consider modern work in the field.
[Bill Poirier](https://aip.scitation.org/doi/full/10.1063/1.4769402) and
[Tucker Carrington](https://aip.scitation.org/doi/full/10.1063/1.3246593) work
are a good place to start.
Special thanks to both Bill Poirier and Tucker Carrington Jr. for their insight
and suggestions during the early phases of this project.

## Authors
Shane W. Flynn, Vladimir A. Mandelshtam. 2019. UCI
