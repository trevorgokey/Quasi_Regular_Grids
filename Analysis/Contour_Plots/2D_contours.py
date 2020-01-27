#=============================================================================80
#                     2D Contour Plot for Model Potentials
#=============================================================================80
#       Discussion:
#Given a set of grid-points and a potential, generate a contour
#Plots for the 2D-Henon-Heiles or 2D Morse Potentials
#==============================================================================#
#To do, make it run for a single grid input to make a single plot.
#test to make sure it works for both potentials
#==============================================================================#
#       Modified:
# 24 January 2020
#       Author:
# Shane Flynn
#==============================================================================#
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#==============================================================================#
def potentials(my_potential,d1,d2):
#==============================================================================#
#Model 2d potential evaluations
#Returns a value of the potential given 2d-coordinates.
#==============================================================================#
#lambda         ==>Scaling Parameter (see literature for value)
#d1/d2          ==>Cartesian Coordinates
#==============================================================================#
    if(lower(my_potential)=='hh'):
        lambd=np.sqrt(0.0125)
        return 0.5*(d1**2+d2**2)+lambd*(d1**2*d2-d2**3/3.)
    elif(lower(my_potential)=='morse'):
        omega_d1=0.2041241
        omega_d2=0.18371169
        D_morse=12.0
        return D_morse*((np.exp(-omega_d1*d1)-1)**2+(np.exp(-omega_d2*d2)-1)**2)
    else:
        sys.exit("Potential Not Recognized. See potentials function.")
#==============================================================================#
#                               Program Main                                   #
#==============================================================================#
#script        ==>The name of THIS script (contour.py)
#grid_#        ==>Data files containing the Grid to superimpose over the contour
#Ncontour      ==>Number of points to construct the Contour
#Ncols         ==>Number of columns in the data file (plot x,y1; x,y2;...x,y_n)
#x/y min/max   ==>Domain for drawing the contours
#x/y X/Y       ==>Grid for evaluating the contour
#Z             ==>Potential at each contour point
#==============================================================================#
#                             Read in Grid Files
#==============================================================================#
script,grid=argv
df_grid=pd.read_csv(grid,delim_whitespace=True,header=None,dtype=np.float64)
Ncols=len(df_grid.columns)
#==============================================================================#
#                              Define Parameters
#==============================================================================#
Ncontour=input()
xmin=input()
xmax=input()
ymin=input()
ymax=input()
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








#=============================================================================80
#                Contour Plot for Morse Potential
#=============================================================================80
#Python Contour Plot of the 2D-Morse Potential
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
script,grid_1,grid_2,grid_3,grid_4=argv
df_1=pd.read_csv(grid_1,delim_whitespace=True,header=None,dtype=np.float64)
df_2=pd.read_csv(grid_2,delim_whitespace=True,header=None,dtype=np.float64)
df_3=pd.read_csv(grid_3,delim_whitespace=True,header=None,dtype=np.float64)
df_4=pd.read_csv(grid_4,delim_whitespace=True,header=None,dtype=np.float64)
#==============================================================================#
#                           Define Parameters
#==============================================================================#
Ncols=len(df_1.columns)
Nint=100
xmin=-5
xmax=21
ymin=-5
ymax=23
#==============================================================================#
#                             Ecut Contour
#==============================================================================#
x=np.linspace(xmin,xmax,Nint)
y=np.linspace(ymin,ymax,Nint)
X,Y=np.meshgrid(x,y)
Z=morse_pot(X,Y)
my_level=[12]
#==============================================================================#
#                             Plot Grid + Contours
title_dict={'fontsize': 32}
label_dict={'fontsize': 30}
my_labels=np.arange(-5,23,10)
#==============================================================================#
fig,axs=plt.subplots(2,2, sharex=True,sharey=True)
axs[0,0].scatter(df_1[0],df_1[1],c='black',s=15)
axs[0,0].contour(X,Y,Z,colors='red',levels=my_level)
axs[0,0].set_title('Direct-product',fontdict=title_dict)
axs[0,0].text(14,20,'N=482',fontdict=label_dict)
#==============================================================================#
axs[0,1].scatter(df_2[0],df_2[1],c='black',s=15)
axs[0,1].contour(X,Y,Z,colors='red',levels=my_level)
axs[0,1].set_title('Unif. quasi-random',fontdict=title_dict)
axs[0,1].text(14,20,'N=482',fontdict=label_dict)
#==============================================================================#
axs[1,0].scatter(df_3[0],df_3[1],c='black',s=15)
axs[1,0].contour(X,Y,Z,colors='red',levels=my_level)
axs[1,0].set_title('Unif. quasi-random+rejection',fontdict=title_dict)
axs[1,0].text(14,20,'N=400',fontdict=label_dict)
#==============================================================================#
axs[1,1].scatter(df_4[0],df_4[1],c='black',s=15)
axs[1,1].contour(X,Y,Z,colors='red',levels=my_level)
axs[1,1].set_title('Quasi-regular',fontdict=title_dict)
axs[1,1].text(14,20,'N=300',fontdict=label_dict)
#==============================================================================#
for ax in axs.flat:
    ax.set(xlabel='', ylabel='')
    for tick in ax.get_xticklabels():
        tick.set_fontsize(28)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(28)
#==============================================================================#
#     Hide x labels and tick labels for top plots and y ticks for right plots
#==============================================================================#
for ax in axs.flat:
    ax.label_outer()
plt.show()
fig.savefig('plot.pdf', bbox_inches='tight')
