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
xmax=20
ymin=-5
ymax=22
#==============================================================================#
#                             Global Contour
#==============================================================================#
x=np.linspace(xmin,xmax,Nint)
y=np.linspace(ymin,ymax,Nint)
X,Y=np.meshgrid(x,y)
Z=morse_pot(X,Y)
#==============================================================================#
#                         Inner Contours (Hard-Coded)
#==============================================================================#
#x1=np.linspace(-7,8,Nint)
#y1=np.linspace(-7,8,Nint)
x1=np.linspace(xmin,xmax,Nint)
y1=np.linspace(ymin,ymax,Nint)
X1,Y1=np.meshgrid(x1,y1)
Z1=morse_pot(X1,Y1)
level=[12]
#==============================================================================#
#                             Plot Grid + Contours
#==============================================================================#
x_1=df_1[0]
x_2=df_2[0]
x_3=df_3[0]
x_4=df_4[0]
#==============================================================================#
fig,axs=plt.subplots(2,2)
my_labels = np.arange(-5,23,10)
my_dict = {'fontsize': 25}
#==============================================================================#
for i in range(1,Ncols):
    axs[0,0].scatter(x_1,df_1[i],c="black",s=14)
axs[0,0].contour(X1,Y1,Z1,colors='red',levels=level)
axs[0,0].set_xticks([])
axs[0,0].set_yticks(my_labels)
axs[0,0].set_yticklabels(labels=my_labels,fontdict=my_dict)
axs[0,0].text(14,20,'N=482',fontdict=my_dict)
#==============================================================================#
for i in range(1,Ncols):
    axs[0,1].scatter(x_2,df_2[i],c="black",s=14)
axs[0,1].contour(X1,Y1,Z1,colors='red',levels=level)
axs[0,1].set_xticks([])
axs[0,1].set_yticks([])
axs[0,1].text(14,20,'N=482',fontdict=my_dict)
#==============================================================================#
for i in range(1,Ncols):
    axs[1,0].scatter(x_3,df_3[i],c="black",s=14)
axs[1,0].contour(X1,Y1,Z1,colors='red',levels=level)
axs[1,0].set_xticklabels(labels=my_labels,fontdict=my_dict)
axs[1,0].set_yticklabels(labels=my_labels,fontdict=my_dict)
axs[1,0].set_xticks(my_labels)
axs[1,0].set_yticks(my_labels)
axs[1,0].text(14,20,'N=400',fontdict=my_dict)
#==============================================================================#
for i in range(1,Ncols):
    axs[1,1].scatter(x_4,df_4[i],c="black",s=14)
axs[1,1].contour(X1,Y1,Z1,colors='red',levels=level)
axs[1,1].set_xticklabels(labels=my_labels,fontdict=my_dict)
axs[1,1].set_xticks(my_labels)
axs[1,1].set_yticks([])
axs[1,1].text(14,20,'N=300',fontdict=my_dict)
#==============================================================================#
plt.show()
fig.savefig('plot.pdf', bbox_inches='tight')
