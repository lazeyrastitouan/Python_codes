from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

CDMl    = loadtxt("/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo11/sftresh_cool_low_fb_boost/AHF_000_allparts/logM_MsqorhomdndM_sMsqorhomdndM_200_snap_z4.000_halos_sorted_main_hires.dat").T
CDMh    = loadtxt("/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo11/AHF_001_allparts/logM_MsqorhomdndM_sMsqorhomdndM_200_snap_z4.000_halos_sorted_main_hires.dat").T
WDM1p5l = loadtxt("/data/titouan/Eagle_cosmo_gizmo/200_512/WDM_1p5keV/zoom1p5_hydro_fire2_lvl7/Halo11/sftresh_cool_low_fb_boost/AHF_000_allparts/logM_MsqorhomdndM_sMsqorhomdndM_200_snap_z4.000_halos_sorted_main_hires.dat").T
WDM1p5h = loadtxt("/data/titouan/Eagle_cosmo_gizmo/200_512/WDM_1p5keV/zoom2_hydro_fire2_lvl7/Halo11/AHF_001_allparts/logM_MsqorhomdndM_sMsqorhomdndM_200_snap_z4.000_halos_sorted_main_hires.dat").T
h         = 0.6777
Om        = 0.307
rho_crit  = 27.75e10  #h^2 Msun/Mpc^3
CDMl[1]    = CDMl[1]/10.**CDMl[0]*rho_crit*Om*h**3.
CDMh[1]    = CDMh[1]/10.**CDMh[0]*rho_crit*Om*h**3.
WDM1p5l[1] = WDM1p5l[1]/10.**WDM1p5l[0]*rho_crit*Om*h**3.
WDM1p5h[1] = WDM1p5h[1]/10.**WDM1p5h[0]*rho_crit*Om*h**3.
CDMl[0]    = 10.**CDMl[0]/h
CDMh[0]    = 10.**CDMh[0]/h
WDM1p5l[0] = 10.**WDM1p5l[0]/h
WDM1p5h[0] = 10.**WDM1p5h[0]/h

fout="/data/titouan/Eagle_cosmo_gizmo/200_512/logM_MsqorhomdndM_CDM_WDM1p5_hydro_fire2_lvl7_Halo11_z4.pdf"
plt.plot(log10(CDMl[0]), log10(CDMl[1]), 'r-',label=r'CDM $m_b=6.7 \cdot 10^5 \, M_\odot$')
plt.plot(log10(CDMh[0]), log10(CDMh[1]), 'r--',label=r'CDM $m_b=10^5 \, M_\odot$')
plt.plot(log10(WDM1p5l[0]), log10(WDM1p5l[1]), 'g-',label=r'WDM 1.5keV $m_b=6.7 \cdot 10^5 \, M_\odot$')
plt.plot(log10(WDM1p5h[0]), log10(WDM1p5h[1]), 'g--',label=r'WDM 1.5keV $m_b=10^5 \, M_\odot$')
#plt.errorbar(WDM[0], log10(WDM[1]), xerr=0.0, yerr=log10((WDM[1]+WDM[2])/WDM[1]),fmt='b.')
plt.xlabel(r'$\log M \, [M_\odot]$',fontsize=20)
plt.ylabel(r'$\log (dn/d\log M) \, [\textrm{Mpc}^{-3}$]', fontsize=20)
plt.legend(loc=0,prop={'size':14})
plt.savefig(fout, bbox_inches='tight')
plt.show()
