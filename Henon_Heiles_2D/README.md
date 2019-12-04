# Henon-Heiles 2D Potential:
Generate a Quasi-Regular Grid for the Henon-Heiles potential (2-Dimensions hard-coded)

## QRG_Generation
Generates the QRG according to the HH-Potential
Computes some of the lower moments for the distribution to determine global accuracy

## Spectra
Computes the eigenspectra for the H-H potential (File with gridpoints for input)

#### Distribution
Uses the underlying target distribution to determine the Gaussian Widths. 
This should be used for grids made with Pseudo-Random and Quasi-Random numbers.

#### Nearest Neighbor
Uses the distance between nearest neighbors to determine the Gaussian Widths. 
This should be used for our QRGs or grids with uniform spacing. 

## Contour
Python script for plotting gridpoints super-imposed on a contour plot of the 
potential. 
Generates a 4-panel figure for showing the scaling wrt. the number of gridpoints. 
