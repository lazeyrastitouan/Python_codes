from numpy import*
import h5py
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import random as rd
from pylab import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

NumT  = 5
Numd  = 1
CPATH = "/level_0/data:datatype=0"
sideL = 512
z     = 3.

f     = h5py.File("/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/CDM/zoom1_hydro/Eagle_cosmo_arepo_200_512_CDM_z3_hydro_zoom1.hdf5","r")
dset  = f[CPATH]
Temp  = dset[(NumT - 1) * sideL ** 3 : (NumT * sideL ** 3)]
Dens  = dset[(Numd - 1) * sideL ** 3 : (Numd * sideL ** 3)]
Dens  = Dens*9.06336e-31/1.6735575e-24*(1.+z)**3. #convert to atom/cm^3

#plot only a subsample for speed
lT    = len(Temp)
subi  = rd.sample(range(lT),lT/10)
subT  = [Temp[i] for i in subi] 
subD  = [Dens[i] for i in subi]

fig = figure()
plt.loglog()
plt.plot(subD,subT,'b.',markersize=4,alpha=0.5)
plt.xlabel(r'Density $[cm^{-3}]$',fontsize=25, labelpad=15)
plt.ylabel(r'Temperature [K]',fontsize=25, labelpad=15)
plt.title(r'CDM, nsfr=10',fontsize=20)
plt.show()
#savefig('/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/CDM/zoom1_hydro/T_vs_n_Halo4_z3_zoom1_CDM.pdf', bbox_inches='tight')
#close(fig)
