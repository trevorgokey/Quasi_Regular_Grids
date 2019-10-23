#=============================================================================80
#                Contour Plot for Morse Potential
#=============================================================================80
#Python Contour Plot of the 2D-Morse Potential
#super-impose grid onto the contour plot.
#Currently makes a 4-panel figure, requires 4 external data sets (x y1 y2...y_n)
#not working so far
#==============================================================================#
#       Modified:
# 10 October 2019
#       Author:
# Shane Flynn
#==============================================================================#
from sys import argv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#==============================================================================#
#       Discussion:
#2D Henon-Heiles Potential
#lambda         ==>Scaling Parameter (see literature for value)
#==============================================================================#
def morse_pot(x,y):
    omega_x=0.2041241
    omega_y=0.18371169
    D_morse=12.0
    return D_morse*((np.exp(-omega_x*x)-1)**2+(np.exp(-omega_y*y)-1)**2)
#==============================================================================#
#                             Program Main
#==============================================================================#
#                              Discussion:
#==============================================================================#
#script        ==>The name of THIS script (contour.py)
#grid_#        ==>Data file containing QR Grid to superimpose over the contour
#Nint          ==>Number of points to construct the Contour
#Ncols         ==>Number of columns in the data file (plot x,y1; x,y2;...x,y_n)
#x/y min/max   ==>Domain for drawing the contours
#x/y X/Y       ==>Grid for evaluating the contour
#Z             ==>Potential at each contour gridpoint
#==============================================================================#
#                    Read in External Grid Points to Data Frame
#==============================================================================#
#script,grid_50,grid_100,grid_200,grid_300=argv
#df_50=pd.read_csv(grid_50,delim_whitespace=True,header=None,dtype=np.float64)
#df_100=pd.read_csv(grid_100,delim_whitespace=True,header=None,dtype=np.float64)
#df_200=pd.read_csv(grid_200,delim_whitespace=True,header=None,dtype=np.float64)
#df_300=pd.read_csv(grid_300,delim_whitespace=True,header=None,dtype=np.float64)
##==============================================================================#
##                           Define Parameters
##==============================================================================#
#Ncols=len(df_50.columns)
Nint=100
xmin=-10
xmax=5
ymin=-10
ymax=5
##==============================================================================#
##                             Global Contour
##==============================================================================#
x=np.linspace(xmin,xmax,Nint)
y=np.linspace(ymin,ymax,Nint)
X,Y=np.meshgrid(x,y)
Z=morse_pot(X,Y)
##==============================================================================#
##                         Inner Contours (Hard-Coded)
##==============================================================================#
#x1=np.linspace(-7,8,Nint)
#y1=np.linspace(-7,8,Nint)
#X1,Y1=np.meshgrid(x1,y1)
#Z1=HH_pot(X1,Y1)
#level=[0.5,5,10,12.9]
##==============================================================================#
##                             Plot Grid + Contours
##==============================================================================#
#x_50=df_50[0]
#x_100=df_100[0]
#x_200=df_200[0]
#x_300=df_300[0]
##==============================================================================#
#fig=plt.figure()
##==============================================================================#
#plt.subplot(2,2,1)
#for i in range(1,Ncols):
#    plt.scatter(x_50,df_50[i],c="black",s=10)
#plt.contour(X,Y,Z,colors='royalblue')
#plt.contour(X1,Y1,Z1,colors='royalblue',levels=level)
#plt.xticks(x,[])
##==============================================================================#
#plt.subplot(2,2,2)
#for i in range(1,Ncols):
#    plt.scatter(x_100,df_100[i],c="black",s=10)
#plt.contour(X,Y,Z,colors='royalblue')
#plt.contour(X1,Y1,Z1,colors='royalblue',levels=level)
#plt.xticks(x,[])
#plt.yticks(y,[])
##==============================================================================#
#plt.subplot(2,2,3)
#for i in range(1,Ncols):
#    plt.scatter(x_200,df_200[i],c="black",s=10)
#plt.contour(X,Y,Z,colors='royalblue')
#plt.contour(X1,Y1,Z1,colors='royalblue',levels=level)
##==============================================================================#
#plt.subplot(2,2,4)
#for i in range(1,Ncols):
#    plt.scatter(x_300,df_300[i],c="black",s=10)
#plt.contour(X,Y,Z,colors='royalblue')
#plt.contour(X1,Y1,Z1,colors='royalblue',levels=level)
#plt.yticks(y,[])
##==============================================================================#
plt.contour(X,Y,Z,colors='royalblue')
plt.show()
#plt.savefig('plot.png', bbox_inches='tight')
