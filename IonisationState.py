from numpy import*
import sys
import os
import h5py
import argparse
import matplotlib.pyplot as plt
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# What to do
parser = argparse.ArgumentParser(description='Plot electron fraction for particles with T < Tthresh,  and neutral hydrogen columnd density for gizmo output high refinement. Args: < path to snapshot (hdf5) (str) >')
parser.add_argument('snappath',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('--Tmax',nargs=1,type=float, help='Temperature max to consider for ne')
args     = parser.parse_args()
snappath = args.snappath[0]
Tmax     = 1e9
if args.Tmax is not None:
	Tmax = args.Tmax[0]

# settings
datatype = 'HDF5' #'Gadget2'
verbose  = False
savenpy  = False
ptype    = 'PartType0'

if verbose:
    print("Datatye = "+datatype)

h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5']
h5files  = sorted(h5files)
h5files  = sorted(h5files,key=len)
if verbose:
	print(h5files)
 
numfiles   = len(h5files)
neA, temperature, mass, nhf = [],[],[],[]
all_neA, all_temperature, all_mass, all_nhf = [],[],[],[]

#load
np = 0
for i in h5files:
    f   = h5py.File(snappath+"/"+i, "r")
    print(f)

    ne = array(f[ptype]['ElectronAbundance']) # #electron/#protons averaged over cell mass 
    nH = array(f[ptype]['NeutralHydrogenAbundance']) #neutral hydrogen fraction
    M  = array(f[ptype]['Masses']) #1e10 Msun/h
    T  = array(f[ptype]['Temperature']) #Kelvin

    numparts = f['Header'].attrs['NumPart_ThisFile'][0]
    for j in range(numparts):
        all_neA.append(ne[j])
        all_temperature.append(T[j])
        all_nhf.append(nH[j])
        all_mass.append(M[j])
        if  T[j] <= Tmax:
            neA.append(ne[j])
            temperature.append(T[j])
            nhf.append(nH[j])
            mass.append(M[j])

#transforming lists to arrays shape=((numfiles,npthifile))
neA             = array(neA)
nhf             = array(nhf)
temperature     = array(temperature)
mass            = array(mass)
all_neA         = array(all_neA)
all_nhf         = array(all_nhf)
all_temperature = array(all_temperature)
all_mass        = array(all_mass)


#[T] = K, transforming density to atoms/pcm^-3
#mH      = 1.6735575e-24 #gram
#mHe     = 6.6464731e-24 #gram
#mA      = 0.76*mH+0.24*mHe
#kpctocm = 3.086e21
#Mtogram = 1.989e43
red     = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['Redshift']
#a       = 1./(1.+red)
#h       = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['HubbleParam']
#dens    = dens*a**(-3.)*h**(-2.)*6.767e-22/mA

#save
if savenpy:
	save(snappath+'/ElectronAbundance', neA)
	save(snappath+'/NeutralHydrogenabundance', nhf)

#do histogram
H     = histogram(log10(neA+1), bins=50)
all_H = histogram(log10(all_neA+1), bins=50)

#plot
plt.loglog()
plt.plot(10**all_H[1][:-1]-1,all_H[0],'r--',label=r"All gas")
plt.plot(10**H[1][:-1]-1,H[0],'b-',label=r"$T<$"+str(Tmax)+" gas")
plt.xlabel(r'Electron Abundance [??]',fontsize=25, labelpad=15)
plt.ylabel(r'Count',fontsize=25, labelpad=15)
plt.title(r'Halo4, z = '+str(round(red,2)),fontsize=20)
plt.legend(loc=0,prop={'size':14})
plt.tight_layout()
plt.show()
