#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
#matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

### figure 1
nc1 = Dataset(sys.argv[2])
nc2 = Dataset(sys.argv[4])
var1=nc1.variables[''.join([sys.argv[1],"Atl"])]
var2=nc2.variables[''.join([sys.argv[3],"Atl"])]
var1y=nc1.variables['lat']
var2y=nc2.variables['lat']
val1=var1[:]
val2=var2[:]
title=''.join([sys.argv[1],"Atl"])
ylabel='Lat ($^o$)'
xlabel=sys.argv[3]
siz=16

fig, ax=plt.subplots()
plt.plot(val2,var2y,linewidth=2,marker='o',label='Obs',color='black')
plt.plot(val1,var1y,linewidth=2,marker='o',label='Model',color='blue')
plt.title(title,size=siz)

ax.legend(loc='upper left')
ax.set_xlabel('',size=siz)
ax.set_ylabel(ylabel,size=siz)
ax.set_xlabel(xlabel,size=siz)
ax.grid(True)

file = sys.argv[1].split("/")[-1].replace(".nc","_Atl.ps")
plt.savefig(file)

### figure 2
nc1 = Dataset(sys.argv[2])
nc2 = Dataset(sys.argv[4])
var1=nc1.variables[''.join([sys.argv[1],"Pac"])]
var2=nc2.variables[''.join([sys.argv[3],"Pac"])]
var1y=nc1.variables['lat']
var2y=nc2.variables['lat']
val1=var1[:]
val2=var2[:]
title=''.join([sys.argv[1],"Pac"])
ylabel='Lat ($^o$)'
xlabel=sys.argv[3]
siz=16

fig, ax=plt.subplots()
plt.plot(val2,var2y,linewidth=2,marker='o',label='Obs',color='black')
plt.plot(val1,var1y,linewidth=2,marker='o',label='Model',color='blue')
plt.title(title,size=siz)
ax.legend(loc='upper left')

ax.set_xlabel('',size=siz)
ax.set_ylabel(ylabel,size=siz)
ax.set_xlabel(xlabel,size=siz)
ax.grid(True)

file = sys.argv[1].split("/")[-1].replace(".nc","_Pac.ps")
plt.savefig(file)
plt.show()
