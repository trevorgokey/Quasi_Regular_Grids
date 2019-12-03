#!/usr/bin/env python
#=============================================================================80
#                     Gaussian Distributed Gridpoints (1D)
#=============================================================================80
#           Discussion: 
# Python 2 script for computing grid points distributed according to a 1D 
# Gaussian basis. 
#==============================================================================#
# Code requires numpy and scipy python packages. 
# Npoints should be an odd number, the gaussian is symmetric. 
# The code assumes x0 := 0, therefore Npoints / 2 are computed and both positive
# and negative values are taken for the grid. 
#==============================================================================#
#       Modified:
# 13 September 2018
#       Author:
# Shane Flynn
#==============================================================================#
from sys import argv
import numpy as np
from scipy.special import erfinv
from scipy.special import erf
#==============================================================================#
#==============================================================================#
#                               Variables
#==============================================================================#
#           Integer:
# Npoints   ==> Number of points to generate (1D)
#           Float:
# x_i       ==> lower_bound  for x dimension
# x_f       ==> upper_bound for x dimension
# area      ==> Total Area under f(x) from x_i to x_f
# bin_area  ==> Area per bin. (Area / Npoints)
# x0        ==> Starting point for the grid
# x         ==> Gaussian Grid Points (positive values)
# y         ==> Gaussian Grid Points (negative values)
# data      ==> Gaussian Grid Points (all grid points)
#==============================================================================#
#        argv = script, lower_bound, upper_bound, Number_Points
#==============================================================================#
script,x_i,x_f,Npoints = argv
x_i = float(x_i)
x_f = float(x_f)
Npoints = int(Npoints)
#==============================================================================#
#                   Area from x_i,x_f, under f(x)
#           int  dx e^{-x^2/2} = sqrt(pi/2)erf(x/sqrt(2)) + c
#==============================================================================#
area = np.sqrt(np.pi/2)*(erf(x_f/np.sqrt(2.))-erf(x_i/np.sqrt(2.)))
#print('Area under Gaussian [x_i,x_f] ==> %r')% area
#==============================================================================#
#                           Width per bin
#==============================================================================#
bin_area = area / Npoints
#print('Area/ bin ==> %r')% bin_area
#==============================================================================#
#                       starting point, x_0 := 0
#==============================================================================#
x0 = 0
x = []
x.append(x0)
#==============================================================================#
#                   Compute grid points from x0 tp x_f
#               Python does not include upper bound (add 1 to index)
#==============================================================================#
for i in range(1,(Npoints+1)/2):
    xi = np.sqrt(2)* erfinv( bin_area / (np.sqrt(np.pi/2)) + erf(x[i-1] / np.sqrt(2))  )
    x.append(xi)
#==============================================================================#
#           Gaussian is symetic, set negative values for [x_i,x0] 
#==============================================================================#
y = []
y.append(-x0)
for i in range(1,(Npoints+1)/2):
    yi = -x[i]
    y.append(yi)
#print y
#==============================================================================#
#                   Remove x0 from y to avoid redundancy
#==============================================================================#
y.remove(y[0])
#print y
#==============================================================================#
#                             Collect all points
#==============================================================================#
data = x + y
#print data
#==============================================================================#
#               Generate 2D Gaussian Distributed Grid
#           Grid is constructed as permutations of two 1D grids
#==============================================================================#
k=0
for i in range(0,(Npoints)):
    for j in range(0,(Npoints)):
        k +=1
        print('%r   %r')%(data[i],data[j])
