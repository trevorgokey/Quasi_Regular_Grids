#=============================================================================80
#                       nD Morse Analytic Eigenvalues
#==============================================================================#
#Compute exact eigenvalues for the nD Morse Oscillator (analytic expression)
#Need to make a general python implementation for computing the analytic
#eigenvalues for the morse oscillator, come back to this!
#==============================================================================#
#       Modified:
#   1 January 2020
#       Author:
#   Shane Flynn
#==============================================================================#
#from sys import argv
#import pandas as pd
#import numpy as np
#==============================================================================#
#                          Collect Morse Parameters
#==============================================================================#
d=int(input())
print(d)
i=0
omega=[]
while i < d:
    omega.append(float(input()))
    i+=1
print(omega)
D_morse=float(input())
print(D_morse)
#The eigenvalues will become negative due to the decreasing distance between
#energy states at higher excitations (After 24 excitations in 1D).
##NB=23
#==============================================================================#
#           Compute 1D eigenvalues for x,y dimensions seperately
#==============================================================================#
##eigx=[]
##eigy=[]
##eigs=[]
##for i in range(0,NB):
##    eigx.append(np.float(ax*np.sqrt(2.*D_morse)*(i+0.5)- \
##            ((ax*np.sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))))
##    eigy.append(np.float(ay*np.sqrt(2.*D_morse)*(i+0.5)- \
##            ((ay*np.sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))))
##print eigx
##print eigy
###==============================================================================#
###                       Compute Excitation Combinatorics
###==============================================================================#
##for i in range(0,NB):
##    for j in range(0,NB):
##        eigs.append(np.float(eigx[i]+eigy[j]))
###==============================================================================#
###                              Sort by Magnitude
###==============================================================================#
##df_theory=pd.DataFrame(eigs)
##df_sort=df_theory.sort_values([0])
##df_sort=df_sort.reset_index(drop=True)
###==============================================================================!
###                           Write Data to CSV File
###==============================================================================!
##df_sort.to_csv('sort_theory.dat',sep='\t',index=False)
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
#               Discussion:
#script             ==> Name of THIS Python script
#theory_values      ==> datafile containing 3DMorse eigenvalue combinatorics
#df_theory          ==> dataframe for theory_values
#df_sort            ==> dataframe with theory_values sorted by magnitude
#==============================================================================!
#                           Read In Theory Values
#==============================================================================!
#script,theory_values=argv
#df_theory=pd.read_csv(theory_values,delim_whitespace=True,header=None,
#dtype=np.float64)
#print('theory input')
#print(df_theory.head())
##==============================================================================!
##                           Sort Theory Eigenvalues
##==============================================================================!
#df_sort=df_theory.sort_values([0])
#df_sort=df_sort.reset_index(drop=True)
#print('Sorted Theory')
#print(df_sort.head())
##==============================================================================!
##                           Write Data to CSV File
##==============================================================================!
#df_sort.to_csv('3Dmorse_sorted.dat',sep='\t',index=False,header=None)
#print('Hello Universe')
