from numpy import*
import matplotlib.pyplot as plt

root='/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo11/sftresh_cool_low_fb_boost/AHF_000_allparts_tmp'
h = 0.6777
M=loadtxt(root+"/200_snap_z4.000_halos_sorted_main_hires",usecols=(3,))/h
M=M/h
logm=log10(M)
xyz=loadtxt(root+"/200_snap_z4.000_halos_sorted_main_hires",usecols=(5,6,7))/1000.
logmhi=logm[where(logm>10.0)]
xyzhi=xyz[where(logm>10.0)]

c2=[103.713,92.897,96.55]

plt.figure(1)
hxy=plt.scatter(xyzhi.T[0], xyzhi.T[1], c=logmhi, s=10, cmap='brg')
#circlexy1 = plt.Circle((c1[0], c1[1]), 0.25, color='k',fill=False)
#circlexy2 = plt.Circle((c2[0], c2[1]), 0.25, color='k',fill=False)
#circlexy3 = plt.Circle((c3[0], c3[1]), 0.25, color='k',fill=False)
#plt.gca().add_patch(circlexy1)
#plt.gca().add_patch(circlexy2)
#plt.gca().add_patch(circlexy3)
cbar=plt.colorbar(hxy)
cbar.set_label('logM')
plt.xlabel('X [cMpc]')
plt.ylabel('Y [cMpc]')

plt.figure(2)
hxz=plt.scatter(xyzhi.T[0], xyzhi.T[2], c=logmhi, s=10, cmap='brg')
#circlexz1 = plt.Circle((c1[0], c1[2]), 0.25, color='k',fill=False)
#circlexz2 = plt.Circle((c2[0], c2[2]), 0.25, color='k',fill=False)
#circlexz3 = plt.Circle((c3[0], c3[2]), 0.25, color='k',fill=False)
#plt.gca().add_patch(circlexz1)
#plt.gca().add_patch(circlexz2)
#plt.gca().add_patch(circlexz3)
cbar=plt.colorbar(hxz)
cbar.set_label('logM')
plt.xlabel('X [cMpc]')
plt.ylabel('Z [cMpc]')

plt.figure(3)
hyz=plt.scatter(xyzhi.T[1], xyzhi.T[2], c=logmhi, s=10, cmap='brg')
#circleyz1 = plt.Circle((c1[1], c1[2]), 0.25, color='k',fill=False)
#circleyz2 = plt.Circle((c2[1], c2[2]), 0.25, color='k',fill=False)
#circleyz3 = plt.Circle((c3[1], c3[2]), 0.25, color='k',fill=False)
#plt.gca().add_patch(circleyz1)
#plt.gca().add_patch(circleyz2)
#plt.gca().add_patch(circleyz3)
cbar=plt.colorbar(hyz)
cbar.set_label('logM')
plt.xlabel('Y [cMpc]')
plt.ylabel('Z [cMpc]')

plt.show()
