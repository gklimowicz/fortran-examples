#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
#matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

# Suggested values for min and max values for plotting

if sys.argv[1]=='pot_temp':
    vmin = -4
    vmax = 28
if  sys.argv[1]=='salt':
    vmin = 24 
    vmax = 40
if  sys.argv[1]=='Gas_Exchange_CO2n':
    vmin = -3
    vmax = 5
if  sys.argv[1]=='oicefr':
    vmin = 0
    vmax = 100

### figure 1
nc = Dataset(sys.argv[2])
var1=nc.variables[sys.argv[1]]
if var1.ndim==3: # 3D fields, e.g. temp and salinity 
    values=var1[:][0]
if var1.ndim==2: # 2D fields, e.g. co2 flux
    values=var1[:]
lons = nc["lon"][:]
lats = nc["lat"][:]
title=''.join([sys.argv[1],"_",nc.xlabel])
clabel=var1.units
#print(nc.xlabel)

lons, lats = np.meshgrid(lons, lats)
m = Basemap(projection='robin', lon_0=0, resolution='c')
x, y = m(lons, lats)

plt.figure()
m.drawcoastlines()
plt.contourf(x, y, values,50,cmap=plt.cm.Spectral_r)
cbar=plt.colorbar(orientation="horizontal")
plt.clim(vmin,vmax) 
siz=16
cbar.set_label(clabel,size=siz,rotation=0,labelpad=20)
plt.title(title)
file=''.join([sys.argv[1],"_",nc.xlabel,".ps"])
plt.savefig(file)

### figure 2 
nc = Dataset(sys.argv[3])
var1=nc.variables[sys.argv[1]]
if var1.ndim==3: #3D fields, e.g. temp and salinity
    values=var1[:][0]
if var1.ndim==2: #2D fields, e.g. co2 flux
    values=var1[:]
lons = nc["lon"][:]
lats = nc["lat"][:]
lons, lats = np.meshgrid(lons, lats)
title = sys.argv[3].split("/")[-1].replace(".nc","")
#clabel=var1.units

m = Basemap(projection='robin', lon_0=0, resolution='c')
x, y = m(lons, lats)
plt.figure()
m.drawcoastlines()
plt.contourf(x, y, values,50,cmap=plt.cm.Spectral_r)
cbar=plt.colorbar(orientation="horizontal")
plt.clim(vmin,vmax)
siz=16
cbar.set_label(clabel,size=siz,rotation=0,labelpad=20)
plt.title(title)
file = sys.argv[3].split("/")[-1].replace(".nc",".ps")
plt.savefig(file)

### figure 3 
nc = Dataset(sys.argv[4])
var1=nc.variables[''.join([sys.argv[1],"_diff"])]
if var1.ndim==3: #3D fields,e.g. temp and salinity
        values=var1[:][0]
if var1.ndim==2: #2D fields, e.g. co2 flux
        values=var1[:]
if np.amin(values)<-1e10: #In case diff file has a value of -1e30
    mask=np.where(values<-1e10)
    values[mask]=np.nan
lons = nc["lon"][:]
lats = nc["lat"][:]
title=''.join([sys.argv[1],"_diff"])
#clabel=var1.units
vmaxc=np.nanmax(values)
lons, lats = np.meshgrid(lons, lats)
m = Basemap(projection='robin', lon_0=0, resolution='c')
x, y = m(lons, lats)

plt.figure()
m.drawcoastlines()
plt.contourf(x, y, values,50,cmap=plt.cm.bwr)
cbar=plt.colorbar(orientation="horizontal")
plt.clim(-vmaxc,vmaxc)
siz=16
cbar.set_label(clabel,size=siz,rotation=0,labelpad=20)
plt.title(title)
file=''.join([sys.argv[1],"_diff",".ps"])
plt.savefig(file)

plt.show()
