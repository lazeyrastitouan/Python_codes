from numpy import*
import h5py
import matplotlib.pyplot as plt
import matplotlib
import os, sys
from matplotlib.colors import*
from matplotlib.colors import LogNorm

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

##### small routine to check the relation between temperature and internal energy in gizmo #####

#load
path    = "/data/titouan/COLIBRE/L0100N1504/z3/colibre_0064/"
h5files = ['colibre_0064.000.hdf5'] #,'colibre_0064.009.hdf5','colibre_0064.019.hdf5','colibre_0064.029.hdf5','colibre_0064.039.hdf5','colibre_0064.049.hdf5','colibre_0064.059.hdf5','colibre_0064.069.hdf5','colibre_0064.079.hdf5']
h5files = sorted(h5files)
h5files = sorted(h5files,key=len)

elecAbundance_gas, u_gas, T_gas, NH = [],[],[],[]
for i in h5files:
    f     = h5py.File(path+i, "r")
    print(f)

    Xet   = array(f['PartType0']['TotalElectronNumberDensities'])
    IEt   = array(f['PartType0']['InternalEnergies'])
    NHt   = f['PartType0']['ElementMassFractions'][:,0]
    Tempt = array(f['PartType0']['Temperatures'])
    
    np = len(IEt)
    for j in range(np):
        elecAbundance_gas.append(Xet[j])
        u_gas.append(IEt[j])
        NH.append(NHt[j])
        T_gas.append(Tempt[j])
    
    del Xet, IEt, NHt, Tempt, np, f	

#convert
mproton                  = 1.6726e-24 #[g]
kB                       = 1.3806e-16 #erg/K
gamma                    = 5./3.
HFrac                    = 0.76
f                        = h5py.File(path+'colibre_0064.000.hdf5','r')
UnitMass_in_g            = f['Units'].attrs['Unit mass in cgs (U_M)'][0]
UnitLength_in_cm         = f['Units'].attrs['Unit length in cgs (U_L)'][0]
UnitTime_in_s            = f['Units'].attrs['Unit time in cgs (U_t)'][0]
UnitVelocity_in_cm_per_s = UnitLength_in_cm/UnitTime_in_s
UnitEnergy_in_cgs        = UnitMass_in_g*UnitLength_in_cm**2/UnitTime_in_s**2 
u_gas                    = array(u_gas)
elecAbundance_gas        = array(elecAbundance_gas)
T_gas                    = array(T_gas)
NH                       = array(NH)

#u_gas = u_gas*4.0/(3.*HFrac + 1.+4.*HFrac*elecAbundance_gas)*mproton/kB*(gamma - 1.)*UnitEnergy_in_cgs/UnitMass_in_g # Temperature in K
u_gas = u_gas*4.0/(3*NH + 1.+4*NH*1.1)*mproton/kB*(gamma - 1)*UnitEnergy_in_cgs/UnitMass_in_g

#do histogram
binst = linspace(min(log10(T_gas)),max(log10(T_gas)),500) 
binsu = linspace(min(log10(u_gas)),max(log10(u_gas)),500)
H, xedges, yedges = histogram2d(log10(T_gas), log10(u_gas),bins=[binst,binsu],density=True)
# some curing of H: remove nan, mask zeros
H           = array(H)
H[isnan(H)] = 0.0
Hmasked     = ma.masked_where(H==0,H) # Mask pixels with a value of 0

#plot
my_cmap = plt.cm.jet
plt.pcolormesh(binst,binsu,Hmasked.T,cmap=my_cmap,norm=LogNorm())
#plt.plot(log10(T_gas),log10(T_gas),'k--')
plt.colorbar()
plt.plot([0,8],[0,8],'k--')
plt.xlabel('Temperature [K]',fontsize=15)
plt.ylabel('Temperature from U [K]',fontsize=15)
plt.show()
