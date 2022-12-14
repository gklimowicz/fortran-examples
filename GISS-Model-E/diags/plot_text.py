import matplotlib.pyplot as plt
import numpy as np
import sys

#Figure 1
siz=16
file = open(sys.argv[1], "r")
ys = []
for line in file:
    ys.append(float(line.split()[-1]))
xs=np.arange(1, len(ys)+1, 1)
fig, ax=plt.subplots()
leg1 = sys.argv[1].split("/")[-1].replace(".txt","")
plt.plot(xs,ys,linewidth=2,marker='o',label=leg1,color='black')

file = open(sys.argv[2], "r")
ys = []
for line in file:
    ys.append(float(line.split()[-1]))
xs=np.arange(1, len(ys)+1, 1)
leg1 = sys.argv[2].split("/")[-1].replace(".txt","")
plt.plot(xs,ys,linewidth=2,marker='o',label=leg1,color='blue')

file = open(sys.argv[3], "r")
ys = []
for line in file:
    ys.append(float(line.split()[-1]))
xs=np.arange(1, len(ys)+1, 1)
leg1 = sys.argv[3].split("/")[-1].replace(".txt","")
plt.plot(xs,ys,linewidth=2,marker='o',label=leg1,color='green')

ax.legend(loc='best')
ax.set_ylabel('Sv',size=siz)
ax.set_xlabel('Years',size=siz)
ax.grid(True)

#file = sys.argv[1].split("/")[-1].replace(".txt",".ps")
file='Transposrt.ps'
plt.savefig(file)
plt.title('Transport')
plt.show()

