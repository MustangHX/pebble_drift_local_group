#!bin/python
import matplotlib.pyplot as plt
from readcol import *

a1,a2=readcol("vr_estimate.txt",twod=False)
plt.plot(a1,a2,ls="--",label="estimate")

a1,a2=readcol("drift_velocity.txt",twod=False)
plt.plot(a1,a2,label="calculate")

plt.show()

