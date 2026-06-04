from numpy import*
import matplotlib.pyplot as plt
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

#PS
PLCDM  = loadtxt("/data/titouan/Eagle_cosmo_gizmo/Eagle_cosmo_matterpower_z0.dat").T
PLWDM3 = loadtxt("/data/titouan/Eagle_cosmo_gizmo/Eagle_cosmo_WDM_3p5kev_matterpower_z0.dat").T
PLWDM1 = loadtxt("/data/titouan/Eagle_cosmo_gizmo/Eagle_cosmo_WDM_1p5kev_matterpower_z0.dat").T
#PCDM  = loadtxt("200_512/zoom1/snapdir_004/Pmm_snap_004",usecols=(0,1)).T
#PWDM  = loadtxt("200_512/WDM_3p5/zoom1/snapdir_004/Pmm_snap_004",usecols=(0,1)).T 
#PCDM4 = loadtxt("200_512/zoom1/snapdir_000/Pmm_snap_000",usecols=(0,1)).T
#PWDM4 = loadtxt("200_512/WDM_3p5/zoom1/snapdir_000/Pmm_snap_000",usecols=(0,1)).T

plt.figure(1)
plt.loglog()
plt.plot(PLWDM3[0],PLWDM3[1],color='orange',ls='-',label=r'WDM 3.5 keV')
plt.plot(PLWDM1[0],PLWDM1[1],'r-',label=r'WDM 1.5 keV')
plt.plot(PLCDM[0],PLCDM[1],'k-',label=r'CDM')
plt.legend(loc=0,prop={'size':14})
plt.xlim([1e-2,8e1])
plt.ylim([4e-6,4e4])
plt.xlabel(r'$k \; [h/\textrm{Mpc}]$',fontsize=25,labelpad=15)
plt.ylabel(r'$P_L(k) \; [(h/\textrm{Mpc})^3]$',fontsize=25,labelpad=15)
plt.tight_layout()
plt.savefig('/data/titouan/Eagle_cosmo_gizmo/PL_CDM_WDM_z0.pdf')

ratioL = PLWDM3[1]/PLCDM[1]-1.
#ratio  = PWDM[1]/PCDM[1]-1.
#ratio4 = PWDM4[1]/PCDM4[1]-1.

#plt.figure(2)
#plt.xscale('log')
#plt.plot([1e-5,1e2],[0,0],'k-')
#plt.plot(PLCDM[0],ratioL,'b-',label=r"Linear")
#plt.plot(PCDM4[0],ratio4,'m:',label=r'$z=4$')
#plt.plot(PCDM[0],ratio,'g--',label=r'$z=3$')
#plt.legend(loc=0,prop={'size':14})
#plt.xlabel(r'$k \; [h/\textrm{Mpc}]$',fontsize=25,labelpad=15)
#plt.ylabel(r'$P^{wdm}(k)/P^{cdm} - 1$',fontsize=25,labelpad=15)

plt.show()
