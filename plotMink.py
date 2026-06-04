from numpy import*
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import*

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

cdm = loadtxt("/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/CDM/zoom1p5_hydro/minkowski_z3_Halo4_XY_w512_s2p_NVox10.dat").T
wdm = loadtxt("/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/WDM_3p5/zoom1p5_hydro/minkowski_z3_Halo4_XY_w512_s2p_NVox10.dat").T

plt.figure(1)
plt.xscale('log')
plt.plot(cdm[0],cdm[1],'b-',marker='o',label=r'CDM')
plt.plot(wdm[0],wdm[1],'b--',marker='o',label=r'WDM')
plt.ylabel(r'F (Area)',fontsize=25, labelpad=15)
plt.xlabel(r"SB threshold [erg/(s cm$^2$ arcsec$^2$)]",fontsize=21, labelpad=15)
plt.legend(loc=0,prop={'size': 14})
savefig('/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/MF1_z3_Halo4_XY_w512_s2p_NVox10.pdf', bbox_inches='tight')
plt.figure(2)
plt.xscale('log')
plt.plot(cdm[0],cdm[2],'r-',marker='o',label=r'CDM')
plt.plot(wdm[0],wdm[2],'r--',marker='o',label=r'WDM')
plt.ylabel(r'U (Boundary)',fontsize=25, labelpad=15)
plt.legend(loc=0,prop={'size': 14})
plt.xlabel(r"SB threshold [erg/(s cm$^2$ arcsec$^2$)]",fontsize=21, labelpad=15)
savefig('/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/MF2_z3_Halo4_XY_w512_s2p_NVox10.pdf', bbox_inches='tight')
plt.figure(3)
plt.xscale('log')
plt.plot(cdm[0],cdm[3],'g-',marker='o',label=r'CDM')
plt.plot(wdm[0],wdm[3],'g--',marker='o',label=r'WDM')
plt.ylabel(r'$\chi$ (Euler characteristic)',fontsize=25, labelpad=15)
plt.xlabel(r"SB threshold [erg/(s cm$^2$ arcsec$^2$)]",fontsize=21, labelpad=15)
plt.legend(loc=0,prop={'size': 14})
savefig('/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/MF3_z3_Halo4_XY_w512_s2p_NVox10.pdf', bbox_inches='tight')
plt.show()
close(figure(1))
close(figure(2))
close(figure(3))
