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
dens_grid=np.ndarray(shape=(n_size,n_rad),dtype=float)
dr=rad0[1]-rad0[0]
print n_rad,n_size
dens=np.ndarray(shape=(n_size,n_rad),dtype=float)
print type(dens)
#print np.shape(dens)
#print np.shape(dens[0])
rad,size=np.meshgrid(rad0,size0)

"""X = 10*np.random.rand(5,3)
dens=np.random.rand(n_size,n_rad)
print dens
fig, ax = plt.subplots()
ax.imshow(X, cmap=cm.jet)#, interpolation='nearest')
print X
plt.show()
"""
lev=np.linspace(0,0,256)
lev2=np.linspace(0,5,256)
plt.rc('text',usetex=True)
for i in range(0,256):
	lev[i]=0.000001*np.exp(i*1.0/256*np.log(1e8/0.000001))
for i in range(0,2000,100):
	num=i*1;
#dens[0],dens[1],dens[2],dens[3],dens[4],dens[5],dens[6],dens[7],dens[8],dens[9],dens[10]=readcol("out_sigma"+str(num)+".txt",twod=False)
	dens=readcol("out_sigma"+str(num)+".txt",twod=False)
#	dens=np.loadtxt("out_sigma"+str(num)+".txt")
	print type(dens)
	print np.shape(rad)
	print np.shape(size)
	mass=0.0
	for i in range(n_size):
		for j in range(n_rad):
			AREA=np.pi*((rad0[j]+dr/2)**2-(rad0[j]-dr/2.0)**2)
			dens[i][j]=dens[i][j]*20.0
			dens_grid[i][j]=dens[i][j]*20.0
			if(dens[i][j]<=1e-3 or True): 
				dens[i][j]+=0.00001
				dens_grid[i][j]+=0.00001
			mass+=AREA*(dens[i][j])*LUNIT*LUNIT

	print mass
#	plt.scatter(rad,size,c=dens)
	plt.imshow(dens,cmap=plt.cm.jet,interpolation='nearest',norm = LogNorm(1e-5,1e9),origin="lower",extent=[1,10.75,0.1,282],aspect='auto')
	plt.yscale("log")
	plt.xlabel(r"\textrm{r (AU)}")
	plt.ylabel(r"\textrm{size (cm)}")
#	plt.colorbar(label=r"$\Sigma \textrm{(g\ cm^{-2}\ dex^{-1})}$")
	plt.colorbar(label=r"\Sigma / {\rm d log} a_p ({\rm g cm^{-2} dex^{-1}})")
	plt.title("t="+str(num)+" yr")
	plt.show()	
#	norm = LogNorm(levels, ncolors=plt.cm.jet.N, clip=True)
#	plt.pcolormesh(rad,size,dens_grid,interpolation='none',cmap=plt.cm.jet,norm = LogNorm(1e-5,1e9))
	plt.scatter(rad,size,c=dens,norm = LogNorm(1e-5,1e9),edgecolor="face")
	plt.yscale("log")
        plt.xlabel(r"\textrm{r (AU)}")
	plt.xlim(min(rad0),max(rad0))
	plt.ylim(min(size0),max(size0))
        plt.ylabel(r"\textrm{size (cm)}")
        plt.colorbar(label=r"$\Sigma \textrm{(g/cm^2)}$",ticks=[0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,1e5,1e6,1e7,1e8])
	plt.title("t="+str(num)+" yr")
        plt.show()
	plt.pcolormesh(rad,size,dens_grid,cmap=plt.cm.jet,norm = LogNorm(1e-5,1e9))
	plt.yscale("log")
        plt.xlabel(r"\textrm{r (AU)}")
#        plt.xlim(min(rad0),max(rad0))
#       plt.ylim(min(size0),max(size0))
        plt.ylabel(r"\textrm{size (cm)}")
        plt.colorbar(label=r"$\Sigma \textrm{(g/cm^2)}$",ticks=[0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,1e5,1e6,1e7,1e8])
        plt.title("t="+str(num)+" yr")
        plt.show()

