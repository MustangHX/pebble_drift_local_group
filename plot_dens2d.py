#!bin/python
from readcol import *
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import sys
num=1
LUNIT=1.49597871e13

rad0=readcol("rad_chart.txt")
size0=readcol("size_chart.txt");
n_rad=len(rad0)
n_size=len(size0)
print n_rad,n_size
dens_grid=np.ndarray(shape=(n_size,n_rad),dtype=float)
dr=rad0[1]-rad0[0]
dens=np.ndarray(shape=(n_size,n_rad),dtype=float)
print dens
#print np.shape(dens)
#print np.shape(dens[0])
rad,size=np.meshgrid(rad0,size0)

lev=np.linspace(0,0,256)
lev2=np.linspace(0,10000,256)
plt.rc('text',usetex=True)
for i in range(0,256):
	lev[i]=0.00001*np.exp(i*1.0/256*np.log(0.1/0.00001))
for i in range(0000,0001,50):
	num=i*1
	num=sys.argv[1]
	dt=sys.argv[2]
#dens[0],dens[1],dens[2],dens[3],dens[4],dens[5],dens[6],dens[7],dens[8],dens[9],dens[10]=readcol("out_sigma"+str(num)+".txt",twod=False)
	dens=readcol("./out_sigma"+str(num)+".txt",twod=False)
#	dens=np.loadtxt("out_sigma"+str(num)+".txt")
	mass=0.0
	for i in range(n_size):
		for j in range(n_rad):
			AREA=np.pi*((rad0[j]+dr/2)**2-(rad0[j]-dr/2.0)**2)
			dens[i][j]=dens[i][j]*20.0
			dens_grid[i][j]=dens[i][j]*20.0
			if(dens[i][j]<=1e-3 and False): 
				dens[i][j]+=0.00000
				dens_grid[i][j]+=0.0000
			mass+=AREA*(dens[i][j])*LUNIT*LUNIT

	print mass
#	plt.figure(figsize=(12,12))
#	plt.scatter(rad,size,c=dens)
#	masked_dens = np.ma.array (dens, mask=np.ma.masked_less(dens,1e-5))
	masked_dens = np.ma.masked_less(dens,1e-2)
	cmap=plt.cm.jet
#cmap.set_bad('w',1.)
	plt.imshow(masked_dens,cmap=cmap,interpolation='nearest',origin="lower",aspect='auto',extent=[1,rad0[n_rad-1],size0[0],size0[n_size-1]],norm=LogNorm())
	plt.yscale("log")
	plt.xlabel(r"\textrm{r (AU)}")
	plt.ylabel(r"$a_p$ ({\rm cm})")
#	plt.colorbar(label=r"$\Sigma \textrm{(g\ cm^{-2}\ dex^{-1})}$")
	plt.colorbar(label=r"\Sigma / {\rm d log} a_p ({\rm g cm^{-2} dex^{-1}})")
	yr=float(num)*float(dt)
	plt.title("t="+str((yr))+" yr")
	k=num
	if k<10:para="./ani/pebble0000"+str(k)
	elif k< 100: para="./ani/pebble000"+str(k)
	elif k< 1000: para="./ani/pebble00"+str(k)
	elif k < 10000: para="./ani/pebble0"+str(k)
	elif k < 100000: para="./ani/pebble"+str(k)
	a1,a2=readcol("1mm.txt",twod=False)
	plt.plot(a1,a2)
	plt.xlim(1.0,10.5)
#	plt.savefig(para+".png")#,dpi=300)
#	plt.clf()
#	plt.close()
	plt.show()
