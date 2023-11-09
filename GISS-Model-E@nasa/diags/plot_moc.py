#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
import numpy.ma as ma
#matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

#vmin = -2
#vmax = 30

### figure 1
nc = Dataset(sys.argv[1])
xs = nc["lato2"][:]
ys = nc["zoce"][:]
values = nc["sf_Atl"][:]
error_value = -6.849315e+26
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contour(X, Y, values,colors='k')
plt.contourf(X, Y, values,50,cmap=plt.cm.Spectral_r)
plt.gca().invert_yaxis()
plt.colorbar()
title = sys.argv[1].split("/")[-1].replace(".nc","_Atl")
plt.title(title)
plt.savefig("mocsection1.ps")

### figure 2
xs = nc["lato2"][:]
ys = nc["zoce"][:]
values = nc["sf_Pac"][:]
error_value = -6.849315e+26
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contour(X, Y, values,colors='k')
plt.contourf(X, Y, values,50,cmap=plt.cm.Spectral_r)
plt.gca().invert_yaxis()
plt.colorbar()
title = sys.argv[1].split("/")[-1].replace(".nc","_Pac")
plt.title(title)
plt.savefig("mocsection2.ps")


### figure 3
xs = nc["lato2"][:]
ys = nc["zoce"][:]
values = nc["sf_Ind"][:]
error_value = -6.849315e+26
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contour(X, Y, values,colors='k')
plt.contourf(X, Y, values,50,cmap=plt.cm.Spectral_r)
plt.gca().invert_yaxis()
plt.colorbar()
title = sys.argv[1].split("/")[-1].replace(".nc","_Ind")
plt.title(title)
plt.savefig("mocsection3.ps")


### figure 4
xs = nc["lato2"][:]
ys = nc["zoce"][:]
values = nc["sf_Glo"][:]
error_value = -6.849315e+26
values = ma.masked_values(values, error_value)
ma.set_fill_value(values, 0)
plt.figure()
X, Y = np.meshgrid(xs, ys)
plt.contour(X, Y, values,colors='k')
plt.contourf(X, Y, values,50,cmap=plt.cm.Spectral_r)
plt.gca().invert_yaxis()
plt.colorbar()
title = sys.argv[1].split("/")[-1].replace(".nc","_Glo")
plt.title(title)
plt.savefig("mocsection4.ps")

plt.show()

