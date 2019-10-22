#=============================================================================80
#                Contour Plot for Henon-Heiles Potential
#=============================================================================80
#Python Contour Plot of the 2D-Henon-Heiles Potential
#super-impose grid onto the contour plot.
#Currently makes a 4-panel figure, requires 4 external data sets (x y1 y2...y_n)
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
def HH_pot(x,y):
    lambd=np.sqrt(0.0125)
    return 0.5*(x**2+y**2)+lambd*(x**2*y-y**3/3.)
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
script,grid_50,grid_100,grid_200,grid_300=argv
df_50=pd.read_csv(grid_50,delim_whitespace=True,header=None,dtype=np.float64)
df_100=pd.read_csv(grid_100,delim_whitespace=True,header=None,dtype=np.float64)
df_200=pd.read_csv(grid_200,delim_whitespace=True,header=None,dtype=np.float64)
df_300=pd.read_csv(grid_300,delim_whitespace=True,header=None,dtype=np.float64)
#==============================================================================#
#                           Define Parameters
#==============================================================================#
Ncols=len(df_50.columns)
Nint=100
xmin=-8
xmax=8
ymin=-6
ymax=8.5
#==============================================================================#
#                             Global Contour
#==============================================================================#
x=np.linspace(xmin,xmax,Nint)
y=np.linspace(ymin,ymax,Nint)
X,Y=np.meshgrid(x,y)
Z=HH_pot(X,Y)
#==============================================================================#
#                         Inner Contours (Hard-Coded)
#==============================================================================#
x1=np.linspace(-7,8,Nint)
y1=np.linspace(-7,8,Nint)
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
fig,axs=plt.subplots(2,2)
#==============================================================================#
my_labels = np.arange(-8,10,2)
my_dict = {'fontsize': 12}
for i in range(1,Ncols):
    axs[0,0].scatter(x_50,df_50[i],c="black",s=10)
axs[0,0].contour(X,Y,Z,colors='royalblue')
axs[0,0].contour(X1,Y1,Z1,colors='royalblue',levels=level)
axs[0,0].contour(X1,Y1,Z1,colors='red',levels=level2)
axs[0,0].set_xticks([])
axs[0,0].set_yticks(my_labels)
axs[0,0].set_yticklabels(labels = my_labels, fontdict = my_dict)
#==============================================================================#
for i in range(1,Ncols):
    axs[0,1].scatter(x_100,df_100[i],c="black",s=10)
axs[0,1].contour(X,Y,Z,colors='royalblue')
axs[0,1].contour(X1,Y1,Z1,colors='royalblue',levels=level)
axs[0,1].contour(X1,Y1,Z1,colors='red',levels=level2)
axs[0,1].set_xticks([])
axs[0,1].set_yticks([])
#==============================================================================#
for i in range(1,Ncols):
    axs[1,0].scatter(x_200,df_200[i],c="black",s=10)
axs[1,0].contour(X,Y,Z,colors='royalblue')
axs[1,0].contour(X1,Y1,Z1,colors='royalblue',levels=level)
axs[1,0].contour(X1,Y1,Z1,colors='red',levels=level2)
axs[1,0].set_xticklabels(labels = my_labels, fontdict = my_dict)
axs[1,0].set_yticklabels(labels = my_labels, fontdict = my_dict)
axs[1,0].set_xticks(my_labels)
axs[1,0].set_yticks(my_labels)
#==============================================================================#
for i in range(1,Ncols):
    axs[1,1].scatter(x_300,df_300[i],c="black",s=10)
axs[1,1].contour(X,Y,Z,colors='royalblue')
axs[1,1].contour(X1,Y1,Z1,colors='royalblue',levels=level)
axs[1,1].contour(X1,Y1,Z1,colors='red',levels=level2)
axs[1,1].set_xticklabels(labels = my_labels, fontdict = my_dict)
axs[1,1].set_xticks(my_labels)
axs[1,1].set_yticks([])
#==============================================================================#
plt.show()
fig.savefig("plot.pdf", bbox_inches = 'tight')
