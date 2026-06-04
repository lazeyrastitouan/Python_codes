from numpy import*
import sys
import os
import h5py
import argparse
#from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# What to do
parser = argparse.ArgumentParser(description='Compute and plot the stellar mass to halo mass relation. < path to groups (hdf5, str) >')
parser.add_argument('group',nargs=1,type=str, help='Path to group files (hdf5)')
args      = parser.parse_args()
grouppath = args.group[0]

verbose = False
hubble  = 0.6777

h5groups = [f for f in os.listdir(grouppath+"/") if os.path.isfile(os.path.join(grouppath+"/", f)) and f[-4:] == 'hdf5']
h5groups = sorted(h5groups)
h5groups = sorted(h5groups,key=len)
if verbose:
	print("Groups "+str(h5groups))

fg0    = h5py.File(grouppath+'/'+h5groups[0],'r')
attrsg = fg0['Header'].attrs
red    = round(attrsg['Redshift'],1)
out    = 'logMh_logMstar_z'+str(red)+'.dat'

Mh     = []
Mstar  = []
nsub   = []
Mssub  = []
Msh    = []
Mc     = []
Mcen   = []

#loop over group files to get DM and stellar masses of groups
print("Loading halos...")
for i in h5groups:
	print("File "+str(h5groups.index(i)+1)+"/"+str(len(h5groups)))
	
	fg    = h5py.File(grouppath+'/'+i,'r')
	Mht   = list(fg['Group/Group_vsDelta_Mass'])
	Nhtf  = len(Mht)
	Mhst  = list(fg['Group/Group_vsDelta_StellarMass'])

	for h in range(Nhtf):
		Mh.append(Mht[h]*1e10/hubble)
		nsub.append(fg['Group/GroupNsubs'][h])
		Msh.append(Mhst[h]*1e10/hubble)
	
	Mssubtemp = fg['Subhalo/SMST'][:,4]*1e10/hubble
	Mcentemp  = fg['Subhalo/SMST'][:,1]*1e10/hubble
	for j in range(len(Mssubtemp)):
		Mssub.append(Mssubtemp[j])
		Mcen.append(Mcentemp[j])

	del Nhtf, Mht, Mssubtemp


print("Done!")
print("Keeping only central galaxies")

Nh = len(Mh)

sub = 0
for h in range(Nh):
	if nsub[h] != 0:
		Mstar.append(Mssub[sub])
		Mc.append(Mcen[sub])
		sub = sub + nsub[h]
	else:
		Mstar.append(Msh[h])
		Mc.append(Mh[h])

print("Done!")

Mh = array(Mh)
Mstar = array(Mstar)

#Engine
engine = loadtxt('/data/titouan/Eagle_cosmo_engine/woAGN/haloMass_stellarStuff_NoAGNL0050N1504_snap12.txt',delimiter=';',usecols=(1,3)).T
#Eagle
eagle  = loadtxt('/data/titouan/Eagle_cosmo_eagle/haloMass_stellarStuff_RefL0100N1504_snap12.txt',delimiter=';',usecols=(1,3)).T

#write
#fout = open(out,'w')
#for i in range(Nh):
#	fout.write(str(Mh[i])+"\t"+str(Mstar[i])+"\n") 
#fout.close()

#model Moster et al. 2009
M1_0  = 10.**11.88
mMo_0 = 0.0282
b_0   = 1.06
b_1   = 0.17
g_0   = 0.556
g_1   = -0.26
nu    = -0.72
mu    = 0.019

z    = 3.
M1z  = 10**(log10(M1_0)*(1.+z)**mu)
mMoz = mMo_0*(z+1.)**nu
gz   = g_0*(z+1.)**g_1
bz   = b_1*z+b_0

Mhm = 10.**linspace(log10(1e7),log10(1e14),100)
Ms  = 2.*Mhm*mMoz/((Mhm/M1z)**(-bz)+(Mhm/M1z)**gz)

#plot Mstar vs Mh
fb = 0.157
ax1 = plt.subplot(111)
#ax1.plot(log10(Mh),log10(Mstar),'r.')
ax1.plot(log10(Mc),log10(Mstar),'r.')
ax1.plot(log10(10.0**array([7.,14.])),log10(fb*10**array([7.,14.])),'k--',label=r'Cosmic baryon fraction')
ax1.plot(log10(Mhm),log10(Ms),'b-',label=r'Moster+ 2009')
ax1.set_xlabel(r'$\log(M) \; [M_\odot]$',fontsize=25,labelpad=15)
ax1.set_ylabel(r'$\log(M_\star) \; [M_\odot]$',fontsize=25,labelpad=15)
ax1.set_title(r'Gizmo, $z = $'+str(red),fontsize=20)
ax1.legend(loc=0,frameon=False,prop={'size':14})
ax1.set_xlim(7.3,13.4)
ax1.set_ylim(5.4,12.2)
#ax2 = plt.subplot(132) 
#ax2.plot(log10(eagle[0]),log10(eagle[1]),'m.')
#ax2.plot(log10(10.0**array([7.,14.])),log10(fb*10**array([7.,14.])),'k--',label=r'Cosmic baryon fraction')
#ax2.plot(log10(Mhm),log10(Ms),'b-',label=r'Moster+ 2009')
#ax2.set_xlabel(r'$\log(M) \; [M_\odot]$',fontsize=25,labelpad=15)
#ax2.set_ylabel(r'$\log(M_\star) \; [M_\odot]$',fontsize=25,labelpad=15)
#ax2.set_title(r'Eagle, $z = 3.017$',fontsize=20)
#ax2.legend(loc=0,frameon=False,prop={'size':14})
#ax2.set_xlim(7.3,13.4)
#ax2.set_ylim(5.4,12.2)
#ax3 = plt.subplot(133) 
#ax3.plot(log10(engine[0]),log10(engine[1]),'y.')
#ax3.plot(log10(10.0**array([7.,14.])),log10(fb*10**array([7.,14.])),'k--',label=r'Cosmic baryon fraction')
#ax3.plot(log10(Mhm),log10(Ms),'b-',label=r'Moster+ 2009')
#ax3.set_xlabel(r'$\log(M) \; [M_\odot]$',fontsize=25,labelpad=15)
#ax3.set_ylabel(r'$\log(M_\star) \; [M_\odot]$',fontsize=25,labelpad=15)
#ax3.set_title(r'Engine, $z = 3.017$',fontsize=20)
#ax3.legend(loc=0,frameon=False,prop={'size':14})
#ax3.set_xlim(7.3,13.4)
#ax3.set_ylim(5.4,12.2)
#plt.savefig('Mstar_Mh_z3_zoom1p5_WDM3p5_FOF.pdf',format='pdf',bbox_inches='tight')
plt.show()
