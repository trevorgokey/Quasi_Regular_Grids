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
