from numpy import*
import h5py
import matplotlib.pyplot as plt

f = h5py.File("200_512/WDM_3p5/zoom1p5_hydro/groups_004/fof_subhalo_tab_004.0.hdf5","r")

#radii
R   = f['Group/Group_R_Mean200'][0]
RT  = f['Group/Group_R_TopHat200'][0]
Rc  = f['Group/Group_R_Crit200'][0] 
Rc1 = f['Group/Group_R_Crit500'][0]

#centre
centre = f['Group/GroupPos'][0]

#Subhalos
subpos = f['Subhalo/SubhaloPos'][:,0:813] 

#plot
circle = plt.Circle((centre[0], centre[1]), R, color='k', fill=False)
cT     = plt.Circle((centre[0], centre[1]), RT,color='g', fill=False)
cc     = plt.Circle((centre[0], centre[1]), Rc,color='r', fill=False)
cc1    = plt.Circle((centre[0], centre[1]), Rc1,color='m', fill=False)

fig    = plt.figure()
ax     = fig.add_subplot(111)
ax.plot(subpos[:,0],subpos[:,1],'rx')  
ax.plot(centre[0], centre[1],'ko')
ax.add_patch(circle) 
ax.add_patch(cT) 
ax.add_patch(cc) 
ax.add_patch(cc1)
ax.set_xlim(66000,72500)
ax.set_ylim(66000,72500) 
ax.set_xlabel("X [kpc/h]")
ax.set_ylabel("Y [kpc/h]")
ax.set_aspect('equal', adjustable='box')
	#####plt.savefig('Halo4_and_subhalos.pdf', bbox_inches='tight')
plt.show()  
