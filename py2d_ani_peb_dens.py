#!bin/python
from readcol import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
num=1
LUNIT=1.49597871e13
rad0=readcol("rad_chart.txt")
size0=readcol("size_chart.txt");
n_rad=len(rad0)
n_size=len(size0)
dr=rad0[1]-rad0[0]
print n_rad,n_size
dens=np.ndarray(shape=(n_size,n_rad),dtype=float)
print dens
#print np.shape(dens)
#print np.shape(dens[0])
rad,size=np.meshgrid(rad0,size0)

X = 10*np.random.rand(5,3)
dens=np.random.rand(n_size,n_rad)
print dens
fig, ax = plt.subplots()
ax.imshow(X, cmap=cm.jet)#, interpolation='nearest')
print X
plt.show()
lev=np.linspace(0,0,256)
lev2=np.linspace(0,5,256)

for i in range(0,256):
	lev[i]=0.000001*np.exp(i*1.0/256*np.log(10.00/0.000001))
for i in range(0,10000,1000):
	num=i*1;
#dens[0],dens[1],dens[2],dens[3],dens[4],dens[5],dens[6],dens[7],dens[8],dens[9],dens[10]=readcol("out_sigma"+str(num)+".txt",twod=False)
	dens=readcol("out_sigma"+str(num)+".txt",twod=False)
#	dens=np.loadtxt("out_sigma"+str(num)+".txt")
	print np.shape(dens)
	mass=0.0
	for i in range(n_size):
		for j in range(n_rad):
			AREA=np.pi*((rad0[j]+dr/2)**2-(rad0[j]-dr/2.0)**2)
			if(dens[i][j]<=1e-3 or True): dens[i][j]+=0.00001
			mass+=AREA*(dens[i][j])*LUNIT*LUNIT

	print mass
#	plt.scatter(rad,size,c=dens)
#	plt.imshow(dens,cmap=plt.cm.jet,interpolation='nearest')
#	plt.yscale("log")
#	plt.xlabel("r (AU)")
#	plt.ylabel("size (cm)")
#	plt.colorbar(label="\Sigma")
#	plt.show()	
	plt.contourf(rad,size,dens,cmap=plt.cm.jet,levels=lev,norm = LogNorm())
#	plt.scatter(rad,size,c=dens,facecolor=None)
	plt.yscale("log")
        plt.xlabel("r (AU)")
        plt.ylabel("size (cm)")
        plt.colorbar(label="\Sigma",ticks=[0.00001,0.0001,0.001,0.01,0.1])
        plt.show()
