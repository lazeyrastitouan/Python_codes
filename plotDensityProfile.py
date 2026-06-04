from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import argparse

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

# ---
#LOAD
H12 = load("/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/2cMpc/P2C_grid_1_256_DensityProfile_n2_Tmin_100000.0_Tmax_1000000000.0.npz") 
H13 = load("/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/2cMpc/P2C_grid_2_256_DensityProfile_n2_Tmin_100000.0_Tmax_1000000000.0.npz")

z = 3.0

# ---
#BETA MODEL
beta = 3.136
rc   = 26.366  #pkpc
ne0  = 5.45 #cm^{-3}
r    = logspace(log10(10.0),log10(50.0),400)
n_e  = ne0 * (1. + (r/rc)**2.)**(-3.*beta/2.)

# ---
#PLOT #1.24 factor to correct for helium fraction when getting ne
fout = "/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/2cMpc/DensityProfile_H1312_zoom1p5_hydro_fire2_lvl7_CDM_z3.pdf"
plt.figure(1)
plt.loglog()
plt.title(r'$m_b \, = \, 7.18 \cdot 10^5 \; M_\odot$',fontsize=14)
plt.vlines([15,23],0,10,colors='grey',ls='--')
plt.plot(r,n_e,'m-',label=r'$\beta \approx 3$ model') #$n_e=5.45\cdot(1.+\left(\frac{r}{26.366}\right)^2)^{-\frac{3\cdot3.136}{2}}$')
#plt.plot(H12['binmean']/(1.+z), H12['dsqrt']*1.24, 'b-',label=r'$M=1.11\cdot10^{13}\,M_\odot$')
#plt.fill_between(H12['binmean']/(1.+z), H12['p16sqrt']*1.24, H12['p84sqrt']*1.24, color='blue',alpha=0.3)
#plt.plot(H13['binmean']/(1.+z), H13['dsqrt']*1.24, 'g-',label=r'$M=6.96\cdot10^{12}\,M_\odot$')
#plt.fill_between(H13['binmean']/(1.+z), H13['p16sqrt']*1.24, H13['p84sqrt']*1.24, color='green',alpha=0.3)
plt.plot(H12['binmean']/(1.+z), H12['density']*1.24, 'b-',label=r'$M=1.11\cdot10^{13}\,M_\odot$')
plt.fill_between(H12['binmean']/(1.+z), H12['p16']*1.24, H12['p84']*1.24, color='blue',alpha=0.3)
plt.plot(H13['binmean']/(1.+z), H13['density']*1.24, 'g-',label=r'$M=6.96\cdot10^{12}\,M_\odot$')
plt.fill_between(H13['binmean']/(1.+z), H13['p16']*1.24, H13['p84']*1.24, color='green',alpha=0.3)
plt.text(40,0.1,r'$z=3.0$',fontsize=12)
#for legend
#plt.plot([0.001,0.001],[1e3,1e3],'k-',label=r'$\sqrt{<n^2>}$')
#plt.plot([0.001,0.001],[1e3,1e3],'k--',label=r'$<n>$')
plt.xlim([H12['binmean'][0],50])
plt.ylim([0.002,3.0])
plt.xlabel(r'R [pkpc]',fontsize=16)
plt.ylabel(r'$n_e$ [cm$^{-3}$]', fontsize=16)
plt.legend(frameon=False,loc=1,prop={'size':12})
plt.savefig(fout, bbox_inches='tight')

plt.show()
