# ==========================================================================
# Theoretical size function of cosmic voids
# ==========================================================================
# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function
 
import numpy as np
import os
import matplotlib.pyplot as plt
import CosmoBolognaLib as cbl
import matplotlib
import voids_functions as vf

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

#Ricciardelli parameters
alpha = [0.0,0.5] #DM, gal -0.19
beta  = [1.86,4.16] #DM, gal

# define a cosmological model, using Planck 2018 parameters #
cosm = cbl.Cosmology(cbl.CosmologicalModel__Planck18_)
cosm.set_parameter(cbl.CosmologicalParameter__sigma8_,0.79)
cosm.set_parameter(cbl.CosmologicalParameter__Omega_matter_,0.29)
cosm.set_parameter(cbl.CosmologicalParameter__Omega_DE_,0.71)
cosm.set_parameter(cbl.CosmologicalParameter__Omega_neutrinos_,0.0)
s8 = cosm.value(cbl.CosmologicalParameter__sigma8_)

# Minimum and maximum of effective void radii
R_0   = 1e-3
R_min = 60.
R_max = 100.

# number of radii at which the size function is computed
n_val = 1000.

# list of effective void radii with logarithmic binning
RR = np.logspace(np.log10(R_0), np.log10(R_max), n_val , endpoint=True)

#position inside void as a fraction of Rv for CDF
vfrac = 1.

# redshift of the samples: 0.3 for LOWZ, 0.55 for CMASS
zz = 0.3
h  = cosm.value(cbl.CosmologicalParameter__hh_)

#Nonlinear void density threshold and effective bias, matching Nadathur BOSSDR11 
vdn_params = {'b_eff':2./s8, 'slope':2.5, 'offset':0., 'deltav_NL':-0.795, 'delta_c':1.686}

#BOSS DR12 volume Veff = 5.1 Gpc3 for CMASS and 2.3 Gpc3 for LOWZ
volume   = 2.3e9*h**3 #(Mpc/h)^3

fv     = []
VSF = np.empty(n_val)
for i in range(n_val):
    # Jennings +13 VSF
    VSF[i] = cosm.size_function(RR[i], zz, "Vdn", vdn_params['b_eff'],vdn_params['slope'],vdn_params['offset'],vdn_params['deltav_NL'],vdn_params['delta_c']) 
VSF=VSF*volume/RR #/RR to go from d/dlnR to d/dR
#normalization
ncumul = np.trapz(VSF*vol(RR),RR,dx=0.01)

#CDF at vfrac*Rv[i]
filter = RR >= R_min
if Vfrac == 1. :
    CDF = 1
else:
    CDF_dm = dprofile_Hamaus(    

#now let's integrate from minimum R
fv.append(vol(RR[filter])*VSF[filter]/RR[filter]/ncumul)
#print(str(ds[j])+' volume filling factor = '+str(fv/ncumul)+' for R = '+str(RR[filter])+' cMpc/h')

    #plt.plot(RR,VSF/RR,'--',label=r''+ds[j])
    #plt.plot(NV[0],NV[1],'b.')

np.savez(fout,R=RR[filter],fv_DR12=fv[0],fv_DR11=fv[1],fv_CDM=fv[2])

#plt.legend()
#plt.show()
