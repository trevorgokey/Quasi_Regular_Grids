#=============================================================================80
#                Contour Plot for Henon-Heiles Potential
#=============================================================================80
#       Discussion:
#Generate Contour Plots for the 2D-Henon-Heiles Potential
#Hard-Coded to make a 4-panel figure (requiring 4 external data sets)
#Show scaling wrt. number of grid points
#Should re-write this code and generalize the analysis...
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
#2D Henon-Heiles Potential (Hard-Coded)
#Returns a value of the potential given x,y coordinates.
#==============================================================================#
#lambda         ==>Scaling Parameter (see literature for value)
#x/y            ==>Coordinates
#==============================================================================#
def HH_pot(x,y):
    lambd=np.sqrt(0.0125)
    return 0.5*(x**2+y**2)+lambd*(x**2*y-y**3/3.)
#==============================================================================#
#==============================================================================#
#==============================================================================#
#                             Program Main
#==============================================================================#
#==============================================================================#
#==============================================================================#
#script        ==>The name of THIS script (contour.py)
#grid_#        ==>Data files containing the Grid to superimpose over the contour
#Ncont         ==>Number of points to construct the Contour
#Ncols         ==>Number of columns in the data file (plot x,y1; x,y2;...x,y_n)
#x/y min/max   ==>Domain for drawing the contours
#x/y X/Y       ==>Grid for evaluating the contour
#Z             ==>Potential at each contour point
#==============================================================================#
#                            Read in Grid Files
#==============================================================================#
script,grid_50,grid_100,grid_200,grid_300=argv
df_50=pd.read_csv(grid_50,delim_whitespace=True,header=None,dtype=np.float64)
df_100=pd.read_csv(grid_100,delim_whitespace=True,header=None,dtype=np.float64)
df_200=pd.read_csv(grid_200,delim_whitespace=True,header=None,dtype=np.float64)
df_300=pd.read_csv(grid_300,delim_whitespace=True,header=None,dtype=np.float64)
#==============================================================================#
#                           Define Parameters
#Assumes all data files have the same number of columns
#==============================================================================#
Ncols=len(df_50.columns)
Ncont=100
xmin=-8
xmax=8
ymin=-6
ymax=8.5
#==============================================================================#
#                             Outer Contours
#==============================================================================#
x=np.linspace(xmin,xmax,Ncont)
y=np.linspace(ymin,ymax,Ncont)
X,Y=np.meshgrid(x,y)
Z=HH_pot(X,Y)
#==============================================================================#
#                    Inner Contours (Hard-Coded Due to Potential)
#==============================================================================#
x1=np.linspace(-7,8,Ncont)
y1=np.linspace(-7,8,Ncont)
X1,Y1=np.meshgrid(x1,y1)
Z1=HH_pot(X1,Y1)
level=[0.5,5,10,12.9]
#==============================================================================#
#                              Ecut(Hard-Coded)
#==============================================================================#
level2=[12.9]
#==============================================================================#
#                             Plot Grid + Contours
#==============================================================================#
x_50=df_50[0]
x_100=df_100[0]
x_200=df_200[0]
x_300=df_300[0]
#==============================================================================#
#                       4 Panel Figure Using Subplots
#==============================================================================#
fig,axs=plt.subplots(2,2)
#==============================================================================#
my_labels=np.arange(-6,8,2)
my_dict={'fontsize': 30}
for i in range(1,Ncols):
    axs[0,0].scatter(x_50,df_50[i],c="black",s=14)
    axs[0,0].contour(X,Y,Z,colors='royalblue')
    axs[0,0].contour(X1,Y1,Z1,colors='royalblue',levels=level)
    axs[0,0].contour(X1,Y1,Z1,colors='red',levels=level2)
    axs[0,0].set_xticks([])
    axs[0,0].set_yticks(my_labels)
    axs[0,0].set_yticklabels(labels = my_labels, fontdict = my_dict)
#==============================================================================#
for i in range(1,Ncols):
    axs[0,1].scatter(x_100,df_100[i],c="black",s=14)
    axs[0,1].contour(X,Y,Z,colors='royalblue')
    axs[0,1].contour(X1,Y1,Z1,colors='royalblue',levels=level)
    axs[0,1].contour(X1,Y1,Z1,colors='red',levels=level2)
    axs[0,1].set_xticks([])
    axs[0,1].set_yticks([])
#==============================================================================#
for i in range(1,Ncols):
    axs[1,0].scatter(x_200,df_200[i],c="black",s=14)
    axs[1,0].contour(X,Y,Z,colors='royalblue')
    axs[1,0].contour(X1,Y1,Z1,colors='royalblue',levels=level)
    axs[1,0].contour(X1,Y1,Z1,colors='red',levels=level2)
    axs[1,0].set_xticklabels(labels = my_labels, fontdict = my_dict)
    axs[1,0].set_yticklabels(labels = my_labels, fontdict = my_dict)
    axs[1,0].set_xticks(my_labels)
    axs[1,0].set_yticks(my_labels)
#==============================================================================#
for i in range(1,Ncols):
    axs[1,1].scatter(x_300,df_300[i],c="black",s=14)
    axs[1,1].contour(X,Y,Z,colors='royalblue')
    axs[1,1].contour(X1,Y1,Z1,colors='royalblue',levels=level)
    axs[1,1].contour(X1,Y1,Z1,colors='red',levels=level2)
    axs[1,1].set_xticklabels(labels = my_labels, fontdict = my_dict)
    axs[1,1].set_xticks(my_labels)
    axs[1,1].set_yticks([])
#==============================================================================#
plt.show()
fig.savefig("plot.pdf",bbox_inches='tight')
