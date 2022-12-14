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
title = sys.argv[2].split("/")[-1].replace(".nc","")
siz=16

fig, ax=plt.subplots()
plt.plot(record,values,linewidth=2,marker='o',label='model',color='black')
plt.title(title,size=siz)

ax.legend(loc='upper left')
ax.set_xlabel('Year',size=siz)
ax.grid(True)
ax.set_xticks(np.arange(1, var1.shape[0]+1, 1))

file = sys.argv[2].split("/")[-1].replace(".nc",".ps")
plt.savefig(file)
plt.show()


