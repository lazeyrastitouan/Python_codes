from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import argparse
import matplotlib.ticker as mticker

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

zfac      = (1.+3.5)**4.

# ---
#LOAD
root = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM"

#Dalya
#ns = 10, Tmin = 0, zoom1p5, fbboost
CDMloz3p5 = load(root+"/gizmo/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")

#COLIBRE
#ns = 10, Tmin = 0, Thermal
colloz3p5 = load(root+"/COLIBRE/L0050N0752/z3p5/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")

#HELLO
Hlo = load(root+"/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")

#Andrea MAGG high-z
SBMHeII = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Andrea_MAGG_data/SB_intr_MAGGwHeII_z3.76_r_sbAndrea_med2575_r_sbstephanie_med2575.npz")

#PLOT
fout = root+"/SBprofile_dalya_zoom1p5_HELLO_Colibre_m6_CDM_obs_singleobj_NB_CompMAGGHiz.pdf"
fig, axs = plt.subplots(1,1)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_ylim([4e-17,5e-14])
axs.set_xlim([40,400])
axs.set_xticks([40,50,70,100,200,300,400],labels=[r'$40$',r'$50$',r'70',r'$100$',r'200',r'300',r'$400$'],minor=True)
axs.xaxis.set_minor_formatter(mticker.ScalarFormatter())
axs.xaxis.set_major_formatter(mticker.ScalarFormatter())
axs.tick_params(bottom=True)
axs.set_xlabel(r'R [ckpc]',fontsize=18)
axs.set_ylabel(r'SB$_{\rm intr}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
axs.text(200,2e-14,r'$n_\star = 10$ cm$^{-3}$',fontsize=18)
#Observations
axs.plot(SBMHeII['rA'],SBMHeII['SBA'],'k-')#,label=r'Observations')
axs.fill_between(SBMHeII['rA'],SBMHeII['p25A'],SBMHeII['p75A'],color='black',alpha=0.2)
#Dalya
axs.plot(CDMloz3p5['binmean'],CDMloz3p5['sbprofile']*zfac, 'b-')
axs.errorbar(CDMloz3p5['binmean'],CDMloz3p5['sbprofile']*zfac,xerr=0.0,yerr=(CDMloz3p5['sbprofile']-0.8*CDMloz3p5['sbprofile'])*zfac,uplims=True,fmt='b')
#Colibre
axs.plot(colloz3p5['binmean'],colloz3p5['sbprofile']*zfac, 'r-')
axs.errorbar(colloz3p5['binmean'],colloz3p5['sbprofile']*zfac,xerr=0.0,yerr=(colloz3p5['sbprofile']-0.8*colloz3p5['sbprofile'])*zfac,uplims=True,fmt='r')
#HELLO
axs.plot(Hlo['binmean'],Hlo['sbprofileX']*zfac, 'g-',alpha=0.5)
axs.errorbar(Hlo['binmean'],Hlo['sbprofileX']*zfac,xerr=0.0,yerr=(Hlo['sbprofileX']-0.8*Hlo['sbprofileX'])*zfac,fmt='g',uplims=True,alpha=0.5)
axs.plot(Hlo['binmean'],Hlo['sbprofileY']*zfac, 'g-',alpha=0.5)
axs.plot(Hlo['binmean'],Hlo['sbprofileZ']*zfac, 'g-',alpha=0.5)
	
print("Saved figure to "+fout)
plt.savefig(fout, bbox_inches='tight')
	
plt.show()
