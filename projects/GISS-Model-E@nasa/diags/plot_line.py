#!/usr/bin/python3

import numpy as np
from netCDF4 import Dataset
import matplotlib
#matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as  plt 
import sys
import os

nc = Dataset(sys.argv[2])
var1=nc.variables[sys.argv[1]]
record=np.arange(1, var1.shape[0]+1, 1)
values=var1[:]
title=''.join([sys.argv[1],"_",nc.xlabel])
ylabel=var1.units
xlabel=var1.dimensions[0]
siz=16

fig, ax=plt.subplots()
plt.plot(record,values,linewidth=2,marker='o',label='Model',color='black')
plt.title(title,size=siz)

ax.legend(loc='upper left')
ax.set_xlabel(xlabel,size=siz)
ax.set_ylabel(ylabel,size=siz)
ax.grid(True)
ax.set_xticks(record)

nc = Dataset(sys.argv[3])
var1=nc.variables[sys.argv[1]]
values=var1[:]
record=np.arange(1, var1.shape[0]+1, 1)

plt.plot(record,values,linewidth=2,marker='o',label='Obs.',color='blue')
plt.legend()
file=''.join([sys.argv[1],".ps"])
plt.savefig(file)

plt.show()
