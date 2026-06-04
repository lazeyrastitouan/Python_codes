from numpy import*
import sys
import os
import h5py
import argparse
import glob
import matplotlib.pyplot as plt
'''
h5files = glob.glob('/data/titouan/COLIBRE/L0050N0752/z3/colibre_0064/colibre_0064.*.hdf5')

U_M       = float(h5py.File(h5files[0],'r')['Units'].attrs['Unit mass in cgs (U_M)'][0])
U_L       = float(h5py.File(h5files[0],'r')['Units'].attrs['Unit length in cgs (U_L)'][0])
U_t       = float(h5py.File(h5files[0],'r')['Units'].attrs['Unit time in cgs (U_t)'][0])
z         = h5py.File(h5files[0],'r')['Header'].attrs['Redshift'][0]
a         = 1./(1.+z)
a_exp_D   = float(h5py.File(h5files[0],'r')['PartType0/Densities'].attrs['a-scale exponent'][0])
U_L_exp_D = float(h5py.File(h5files[0],'r')['PartType0/Densities'].attrs['U_L exponent'][0])
U_M_exp_D = float(h5py.File(h5files[0],'r')['PartType0/Densities'].attrs['U_M exponent'][0])
D_factor  = a**a_exp_D * U_L**U_L_exp_D * U_M**U_M_exp_D
mA        = 1.6735575e-24 * 0.76 + 6.6464731e-24 * 0.24

dust, dens = empty(0),empty(0)
for f in h5files:
    sim   = h5py.File(f,'r')    
    dustt = array(sim['PartType0/TotalDustMassFractions'])
    dust  = concatenate((dust,dustt))
    denst = array(sim['PartType0/Densities']) * D_factor / mA 
    
    dens  = concatenate((dens,denst))

inds = dust.argsort()

dust_sorted = dust[inds]
dens_sorted = dens[inds]

filter = dust_sorted > 0
dust_toplot = dust_sorted[filter][0:10000000]
dens_toplot = dens_sorted[filter][0:10000000]

plt.loglog()
plt.plot(dens_toplot,dust_toplot,'k.')
plt.show()
'''

#save('/data/titouan/COLIBRE/L0050N0752/z3/colibre_0064/dust.npy',dust)
dust = load('/data/titouan/COLIBRE/L0050N0752/z3/colibre_0064/dust.npy')
H, E = histogram(log10(dust[dust>0]),bins=400,density=True)

plt.stairs(H,E,fill=False)
plt.plot([-7,-7],[0,1],'k--')
plt.yscale('log')
plt.savefig('/data/titouan/COLIBRE/L0050N0752/z3/colibre_0064/dust.pdf',bbox_inches='tight')
