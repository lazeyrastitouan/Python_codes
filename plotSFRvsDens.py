from numpy import*
import h5py
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')

snap  = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1_hydro_fire2/snapdir_004/snap_004.'
z     = 3.
a     = 1./(1.+z)
mh    = 1.6735e-24
h     = 0.6777
tocgs = 2.9349269588313416e+43/(4.5531621661502137e+21**3.)*h**2.

sfr  = []
dens = []
mgas = []
for i in range(16):
	f=h5py.File(snap+str(i)+".hdf5",'r')
	print(i)
	denst = list(f['PartType0/Density'])
	sfrt  = list(f['PartType0/StarFormationRate'])
	mgast = list(f['PartType0/Masses'])
	n = len(sfrt)
	for j in range(n):
		sfr.append(sfrt[j])
		dens.append(denst[j])
		mgas.append(mgast[j])
	del n, f, sfrt, denst

mgas = array(mgas)
sfr  = array(sfr)
dens = array(dens)*tocgs/a**3./mh #assumes all atoms are hydrogen

plt.loglog()
plt.plot(dens,sfr,'b.')
plt.xlim(1e0,1e2)
plt.ylabel('SFR',fontsize=20)
plt.xlabel('Density [atom/cm^3]',fontsize=20)
plt.title('z=3, Fire-2',fontsize=25)
plt.savefig('/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1_hydro_fire2/snapdir_004/sfr_vs_density_z3_zoom1_hydro_fire2.pdf',bbox_inches='tight')
#plt.show()
