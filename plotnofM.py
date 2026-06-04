from numpy import*
import matplotlib.pyplot as plt
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

#n200z4	 = loadtxt("200_512/zoom1/groups_000/logM_NofM_fof_subhalo_tab_000.dat").T
n200z3p7 = loadtxt("200_512/CDM/base/groups_001/logM_NofM_fof_subhalo_tab_001.dat").T
n200z3p5 = loadtxt("200_512/CDM/base/groups_002/logM_NofM_fof_subhalo_tab_002.dat").T
n200z3p2 = loadtxt("200_512/CDM/base/groups_003/logM_NofM_fof_subhalo_tab_003.dat").T
n200z3   = loadtxt("200_512/CDM/zoom1p5_hydro/groups_004/logM_NofM_fof_subhalo_tab.dat").T
nWDMz3   = loadtxt("200_512/WDM_3p5/zoom1/groups_004/logM_NofM_fof_subhalo_tab_004.dat").T
#nWDMz4   = loadtxt("200_512/WDM_3p5/zoom1/groups_000/logM_NofM_fof_subhalo_tab_000.dat").T
nAHF     = loadtxt("200_512/CDM/zoom1p5_hydro/AHF/logM_NofM_200_snap_z3.000_halos_sorted.dat").T

plt.yscale("log")
plt.plot(n200z3[0],n200z3[1],'r-',label=r"FOF, $z=3$")
#plt.plot(n200z3p2[0],n200z3p2[1],'y-',label=r"$200 \,\textrm{Mpc}, \, z=3.2$")
#plt.plot(n200z3p5[0],n200z3p5[1],'b-',label=r"$200 \,\textrm{Mpc}, \, z=3.5$")
#plt.plot(n200z3p7[0],n200z3p7[1],'g-',label=r"$200 \,\textrm{Mpc}, \, z=3.7$")
plt.plot(nAHF[0],nAHF[1],'r--',label=r"AHF, $z=3$")
#plt.plot(nWDMz3[0],nWDMz3[1],'m-',label=r"WDM, $z=3$")
#plt.plot(nWDMz4[0],nWDMz4[1],'m--',label=r"WDM, $z=4$")
plt.legend(loc=0,prop={'size':14})
plt.xlabel(r"$\log M \; [M_\odot]$",fontsize=24, labelpad=14)
plt.ylabel(r"$N(\log M)$",fontsize=26,labelpad=15)
plt.show()
