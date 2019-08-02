#=============================================================================80
#               Distributed Gaussian Basis (Eigen Error Analysis)
#python code for computing and sorting all 2dmorse eigenvalues in 1 script
#current code gives different results compared to fortran after 6 decimals, 
#need to figure out why there is a descrepancy, use this coe at your own risk!
#==============================================================================#
#       Discussion:
#Compute and sort eigenvalues for the 2D Morse Oscillator (analytic)
#==============================================================================#
#       Modified:
#   5 May 2019
#       Author:
#   Shane Flynn 
#==============================================================================#
import pandas as pd
import numpy as np
#==============================================================================#
#                               2D Morse Parameters
#==============================================================================#
D_morse=12
ax=0.2041241  
ay=0.18371169 
NB=23
#==============================================================================#
#           Compute 1D eigenvalues for x,y dimensions seperately
#==============================================================================#
eigx=[]
eigy=[]
eigs=[]
for i in range(0,NB):
    eigx.append(np.float(ax*np.sqrt(2.*D_morse)*(i+0.5)-((ax*np.sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))))
    eigy.append(np.float(ay*np.sqrt(2.*D_morse)*(i+0.5)-((ay*np.sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))))
print eigx
print eigy
#==============================================================================#
#                       Compute Excitation Combinatorics
#==============================================================================#
for i in range(0,NB):
    for j in range(0,NB):
        eigs.append(np.float(eigx[i]+eigy[j]))
#==============================================================================#
#                              Sort by Magnitude
#==============================================================================#
df_theory=pd.DataFrame(eigs)
df_sort=df_theory.sort_values([0])
df_sort=df_sort.reset_index(drop=True)
#==============================================================================!
#                           Write Data to CSV File
#==============================================================================!
df_sort.to_csv('sort_theory.dat',sep='\t',index=False)
