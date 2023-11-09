#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
#matplotlib.use("Agg")
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

#vmin = -2
#vmax = 30

### figure 1
nc = Dataset(sys.argv[2])
var1=nc.variables[''.join([sys.argv[1],"Atl"])]
values=var1[:]
xs = nc["lat"][:]
ys = nc["zoc"][:]
error_value = -1.e30
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
#title = sys.argv[2].split("/")[-1].replace(".nc","")
title = ''.join([sys.argv[1],"Atl"])
ylabel=var1.dimensions[0]
xlabel=var1.dimensions[1]
#clabel=var1.units

siz=16
siz_tick=14
X, Y = np.meshgrid(xs, ys)

fig, ax = plt.subplots()
plt.contour(X, Y, values,colors='k')
plt.contourf(X, Y, values,50,cmap=plt.cm.Spectral_r)
plt.gca().invert_yaxis()
plt.title(title,size=siz)
plt.tick_params(labelsize=siz_tick)

ax.set_xlabel(xlabel,size=siz)
ax.set_ylabel(ylabel,size=siz)

cbar=plt.colorbar()
#cbar.set_label(clabel,size=siz,rotation=270,labelpad=30)

file=''.join([sys.argv[1],"_Atl",".ps"])
plt.savefig(file)


### figure 2
nc = Dataset(sys.argv[2])
var1=nc.variables[''.join([sys.argv[1],"Pac"])]
values=var1[:]
xs = nc["lat"][:]
ys = nc["zoc"][:]
error_value = -1.e30
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
#title = sys.argv[2].split("/")[-1].replace(".nc","")
title = ''.join([sys.argv[1],"Pac"])
ylabel=var1.dimensions[0]
xlabel=var1.dimensions[1]

siz=16
siz_tick=14
X, Y = np.meshgrid(xs, ys)

fig, ax = plt.subplots()
plt.contour(X, Y, values,colors='k')
plt.contourf(X, Y, values,50,cmap=plt.cm.Spectral_r)
plt.gca().invert_yaxis()
plt.title(title,size=siz)
plt.tick_params(labelsize=siz_tick)

ax.set_xlabel(xlabel,size=siz)
ax.set_ylabel(ylabel,size=siz)
cbar=plt.colorbar()

file=''.join([sys.argv[1],"_Pac",".ps"])
plt.savefig(file)

plt.show()


