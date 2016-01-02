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
v_rel=np.ndarray(shape=(n_size,n_size),dtype=float)
print v_rel
#print np.shape(dens)
#print np.shape(dens[0])
rad,size=np.meshgrid(size0,size0)

lev=np.linspace(0,0,256)
lev2=np.linspace(0,10000,256)
plt.rc('text',usetex=True)
for i in range(0,256):
	lev[i]=0.00001*np.exp(i*1.0/256*np.log(0.1/0.00001))
for i in range(0000,0001,50):
	num=i*1
#	num=sys.argv[1]
#`	dt=sys.argv[2]
#dens[0],dens[1],dens[2],dens[3],dens[4],dens[5],dens[6],dens[7],dens[8],dens[9],dens[10]=readcol("out_sigma"+str(num)+".txt",twod=False)
	v_rel=readcol("./relative_velocity.txt",twod=False)
#	dens=np.loadtxt("out_sigma"+str(num)+".txt")
       	for i in range(n_size):
	       for j in range(n_size):
		       v_rel[i][j]=v_rel[i][j]*1.0
				         
	mass=0.0
	sigma_peb=np.linspace(0,0,n_rad)
	print v_rel
#	plt.figure(figsize=(12,12))
#	plt.scatter(rad,size,c=dens)
#	masked_dens = np.ma.array (dens, mask=np.ma.masked_less(dens,1e-5))
#	masked_dens = np.ma.masked_less(dens,1e-2)
	cmap=plt.cm.jet
#cmap.set_bad('w',1.)
	fig,ax1=plt.subplots()
#im=ax1.imshow(masked_dens,cmap=cmap,interpolation='nearest',origin="lower",aspect='auto',extent=[1,rad0[n_rad-1],size0[0],size0[n_size-1]],norm=LogNorm())
	im=ax1.contourf(rad,size,v_rel,cmap=plt.cm.jet)

	ax1.set_yscale("log")
	ax1.set_xscale("log")
	ax1.set_ylim(size0[0],size0[n_size-1])
	ax1.set_xlim(size0[0],size0[n_size-1])
	ax1.set_xlabel(r"$a_p$ ({\rm cm}")
	ax1.set_ylabel(r"$a_p$ ({\rm cm})")
#	plt.colorbar(label=r"$\Sigma \textrm{(g\ cm^{-2}\ dex^{-1})}$")
	fig.subplots_adjust(right=0.75)
	cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
#fig.colorbar(im,cax=cbar_ax,label=r"\Sigma / {\rm d log} a_p ({\rm g cm^{-2} dex^{-1}})")
	fig.colorbar(im,cax=cbar_ax,label=r"\delta v ({\rm cm/s})")
#	yr=float(num)*float(dt)
#	ax1.set_title("t="+str((yr))+" yr")
	k=num
	if k<10:para="./ani/pebble0000"+str(k)
	elif k< 100: para="./ani/pebble000"+str(k)
	elif k< 1000: para="./ani/pebble00"+str(k)
	elif k < 10000: para="./ani/pebble0"+str(k)
	elif k < 100000: para="./ani/pebble"+str(k)
	plt.show()
