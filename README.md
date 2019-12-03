# Quasi-Regular Grids; Applications to the Vibrational Spectra Calculations
Fortran implementation for generating a QRG using model Potential Energy
Functions.
The program generates a set of points optimally distributed for sampling a
potential of interest.
These grids are used to compute the associated vibrational eigenspectra.
For more details see our publication Shane Flynn and Vladimir Mandelshtam
(add link here).

## References
Accurate RoVibrational Energy calculations have a long history in computational
chemistry.
These references will help give context for the application we used to
demonstrate the utility of our QRGs.

#### Distributed Gaussian Basis Sets
The project was motivated in part by the work of [Garashchuk and Light](https://aip.scitation.org/doi/abs/10.1063/1.1348022).
Special thanks to Garashcuk for discussing her previous results with us and
helping compare our findings.

The work by [Poirier and Light] (https://aip.scitation.org/doi/abs/10.1063/1.481787)
is a good reference for using Distributed Gaussian Basis sets in the context of
RoVibrational spectroscopy.

Interested readers should also consider modern work in the field, the work by
[Bill Poirier] (https://aip.scitation.org/doi/full/10.1063/1.4769402) and
[Tucker Carrington] (https://aip.scitation.org/doi/full/10.1063/1.3246593) is a
good place to start.
Special thanks to both authors for their insight and suggestions during the
early phases of this project.

#### Quasirandom Sequences:
The Monte Carlo and Quasi-Monte Carlo Wiki page ([MCQMC](http://roth.cs.kuleuven.be/wiki/Main_Page)) is a useful resource for generating and applying Monte Carlo methods.
The [Sobol Sequence Generator](https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html) was used for the Quasirandom Calculations (sobol.f90).
Thanks to John Burkardt for general advice on quasi-random numbers and their
implementation.


## Authors
Shane W. Flynn, Vladimir A. Mandelshtam. 2019. UCI
