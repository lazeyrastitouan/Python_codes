from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import argparse
import matplotlib.ticker as mticker

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)

checkAGN  = False
checkns   = True
checkzm   = False
checkT    = False
checkmet  = True
zfac      = (1.+3.5)**4.

# ---
#LOAD
root = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM"

#Dalya
#ns = 10, Tmin = 0, zoom1p5, fbboost
CDMloz3   = load(root+"/gizmo/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
CDMhiz3   = load(root+"/gizmo/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.5_Mhigh_13.1.npz")
CDMloz3p5 = load(root+"/gizmo/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
CDMhiz3p5 = load(root+"/gizmo/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
if checkns:
	#ns = 3
	CDMns3loz3    = load(root+"/gizmo/nstar3/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	CDMns3hiz3    = load(root+"/gizmo/nstar3/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.5_Mhigh_13.1.npz")
	CDMns3loz3p5  = load(root+"/gizmo/nstar3/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	CDMns3hiz3p5  = load(root+"/gizmo/nstar3/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	#ns = 1
	CDMns1loz3    = load(root+"/gizmo/nstar1/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	CDMns1hiz3    = load(root+"/gizmo/nstar1/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.5_Mhigh_13.1.npz")
	CDMns1loz3p5  = load(root+"/gizmo/nstar1/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	CDMns1hiz3p5  = load(root+"/gizmo/nstar1/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	#ns = 0.1
	CDMns01loz3   = load(root+"/gizmo/nstar01/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
	CDMns01hiz3   = load(root+"/gizmo/nstar01/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-20_Mmin_12.5_Mhigh_13.1.npz")
	CDMns01loz3p5 = load(root+"/gizmo/nstar01/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
	CDMns01hiz3p5 = load(root+"/gizmo/nstar01/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.5_Mhigh_13.0.npz")
if checkT:
	CDMTlo = load(root+"/gizmo/Tmin1e3/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	CDMThi = load(root+"/gizmo/Tmin1e3/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
if checkAGN:
	CDMfb  = load(root+"/gizmo/sftresh_cool_low/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_3.0_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
if checkzm:
	CDMzm2 = load(root+"/gizmo/SBprofile_zoom2_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_3.0_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
	if checkDM:
		WDM1p5zm2  = load(root+"/gizmo/SBprofile_zoom2_hydro_fire2_lvl7_WDM_1p5keV_observed_cubex_neb_z_3.0_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
	if checkns:
		CDMzm2ns01 = load(root+"/gizmo/nstar01/SBprofile_zoom2_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_3.0_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
if checkmet:
	CDMmetloz3p5 = load(root+'/gizmo/metal005/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz')


#COLIBRE
#ns = 10, Tmin = 0, Thermal
colloz3   = load(root+"/COLIBRE/L0050N0752/z3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.2_Mhigh_12.5.npz")
colhiz3   = load(root+"/COLIBRE/L0050N0752/z3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
colloz3p5 = load(root+"/COLIBRE/L0050N0752/z3p5/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
colhiz3p5 = load(root+"/COLIBRE/L0050N0752/z3p5/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
if checkns:
	#ns = 3
	colns3loz3    = load(root+"/COLIBRE/L0050N0752/z3/nstar3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.2_Mhigh_12.5.npz")
	colns3hiz3    = load(root+"/COLIBRE/L0050N0752/z3/nstar3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	colns3loz3p5  = load(root+"/COLIBRE/L0050N0752/z3p5/nstar3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	colns3hiz3p5  = load(root+"/COLIBRE/L0050N0752/z3p5/nstar3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	#ns = 1
	colns1loz3    = load(root+"/COLIBRE/L0050N0752/z3/nstar1/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.2_Mhigh_12.5.npz")
	colns1hiz3    = load(root+"/COLIBRE/L0050N0752/z3/nstar1/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	colns1loz3p5  = load(root+"/COLIBRE/L0050N0752/z3p5/nstar1/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	colns1hiz3p5  = load(root+"/COLIBRE/L0050N0752/z3p5/nstar1/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	#ns = 0.1
	colns01loz3   = load(root+"/COLIBRE/L0050N0752/z3/nstar01/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-20_Mmin_12.2_Mhigh_12.5.npz")
	colns01hiz3   = load(root+"/COLIBRE/L0050N0752/z3/nstar01/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-20_Mmin_12.5_Mhigh_13.0.npz")
	colns01loz3p5 = load(root+"/COLIBRE/L0050N0752/z3p5/nstar01/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
	colns01hiz3p5 = load(root+"/COLIBRE/L0050N0752/z3p5/nstar01/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.5_Mhigh_13.0.npz")
	#dust-7
	coldm7loz3    = load(root+"/COLIBRE/L0050N0752/z3/dust-7/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-20_Mmin_12.2_Mhigh_12.5.npz") 
	coldm7hiz3    = load(root+"/COLIBRE/L0050N0752/z3/dust-7/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.0_s1_n5e-20_Mmin_12.5_Mhigh_13.0.npz") 
	coldm5loz3p5  = load(root+"/COLIBRE/L0050N0752/z3p5/dust-5/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
	coldm5hiz3p5  = load(root+"/COLIBRE/L0050N0752/z3p5/dust-5/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.5_Mhigh_13.0.npz")
if checkT:
	colTlo        = load(root+"/COLIBRE/L0050N0752/z3p5/Tmin1e3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	colThi        = load(root+"/COLIBRE/L0050N0752/z3p5/Tmin1e3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
if checkmet:
	colmetloz3p5  = load(root+'/COLIBRE/L0050N0752/z3p5/metal005/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.5_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz')

#HELLO
Hlo = load(root+"/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
Hhi = load(root+"/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
if checkns:
	#ns = 3
	Hns3lo  = load(root+"/HELLO/CDM/2cMpc/nstar3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	Hns3hi  = load(root+"/HELLO/CDM/2cMpc/nstar3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	#ns = 1
	Hns1lo  = load(root+"/HELLO/CDM/2cMpc/nstar1/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	Hns1hi  = load(root+"/HELLO/CDM/2cMpc/nstar1/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
	#ns = 0.1
	Hns01lo = load(root+"/HELLO/CDM/2cMpc/nstar01/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz")
	Hns01hi = load(root+"/HELLO/CDM/2cMpc/nstar01/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-20_Mmin_12.5_Mhigh_13.0.npz")
if checkT:
	HTlo   = load(root+"/HELLO/CDM/2cMpc/Tmin1e3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.0_Mhigh_12.5.npz")
	HThi   = load(root+"/HELLO/CDM/2cMpc/Tmin1e3/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-19_Mmin_12.5_Mhigh_13.0.npz")
if checkmet:
	Hmetlo = load(root+'/HELLO/CDM/2cMpc/metal005/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n5e-20_Mmin_12.0_Mhigh_12.5.npz')

#Borisova/MAGG/MUSEUM
SBMz4   = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_4.1.npz")
SBMz3p5 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_3.5.npz")
SBb     = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Borisova_SBprofile_observed_fromMariono_z_3.0.npz")
SBMz3   = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MUSEUM_FAB_SBprofile_observed_z_3.2.npz")
SBMHeII = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Andrea_MAGG_data/SB_intr_MAGGwHeII_z3.76_r_sbAndrea_med2575_r_sbstephanie_med2575.npz")

# ---
#PLOT - Effect of nstar for skipISM
if checkns and not checkT and not checkmet:
	fout = root+"/SBprofile_dalya_zoom1p5_HELLO_Colibre_m6_CDM_obs_singleobj_NB_Impact_nstar.pdf"
	fig, axs = plt.subplots(2,3,sharey=True,sharex=True,figsize=(14,7))
	plt.subplots_adjust(wspace=0.07,hspace=0.07)
	for ax in fig.get_axes():
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylim([5e-20,2e-16])
	for ax in fig.get_axes()[0:5]:
		ax.set_xlim([40,400])
		ax.set_xticks([40,50,70,100,200,300,400],labels=[r'$40$',r'$50$',r'70',r'$100$',r'200',r'300',r'$400$'],minor=True)
		ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
		ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
		ax.tick_params(bottom=True)
	for ax in fig.get_axes()[0:3]:
		ax.text(43,6e-20,r'$z=3.5$',fontsize=12)
		#MAGG
		if ax == axs[0,1]:
			ax.plot(SBMz3p5['r'],SBMz3p5['sbprofile'],'k:',label=r'MAGG')
		else:
			ax.plot(SBMz3p5['r'],SBMz3p5['sbprofile'],'k:')
		ax.fill_between(SBMz3p5['r'],SBMz3p5['p25'],SBMz3p5['p75'],color='black',alpha=0.2)
	for ax in fig.get_axes()[3:5]:
		ax.text(43,6e-20,r'$z=3$',fontsize=12)
		#MUSEUM & Borisova
		if ax == axs[1,1]:
			ax.plot(SBMz3['r'],SBMz3['sbprofile'],'k:',label=r'MUSEUM FAB')
			ax.plot(SBb['r'],SBb['sbprofile'],'k--',label=r'Borisova+ 2016')
		else:
			ax.plot(SBMz3['r'],SBMz3['sbprofile'],'k:')
			ax.plot(SBb['r'],SBb['sbprofile'],'k--')
		ax.fill_between(SBMz3['r'],SBMz3['p25'],SBMz3['p75'],color='black',alpha=0.2)
		ax.fill_between(SBb['r'],SBb['p25'],SBb['p75'],color='black',alpha=0.2)
	fig.supylabel(r'SB$_{\rm obs}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=16)
	axs[1,2].set_xlabel(r'R [ckpc]',fontsize=14)
	axs[1,0].set_xlabel(r'R [ckpc]',fontsize=14)
	axs[1,1].set_xlabel(r'R [ckpc]',fontsize=14)
	axs[1,2].set_xlabel(r'$n_\star$ [cm$^{-3}$]',fontsize=14)
	axs[0,0].set_title(r'DaLya',fontsize=14)
	axs[0,1].set_title(r'COLIBRE',fontsize=14)
	axs[0,2].set_title(r'HELLO',fontsize=14)

	#axs[0,:] z = 3.5
	#axs[0,0] Dalya
	axs[0,0].plot(CDMloz3p5['binmean'],CDMloz3p5['sbprofile'], 'b-')
	axs[0,0].plot(CDMns3loz3p5['binmean'],CDMns3loz3p5['sbprofile'], 'b--')
	axs[0,0].plot(CDMns1loz3p5['binmean'],CDMns1loz3p5['sbprofile'], 'b-.')
	axs[0,0].plot(CDMns01loz3p5['binmean'],CDMns01loz3p5['sbprofile'], 'b:')
	#axs[0,1] COLIBRE
	axs[0,1].plot(colloz3p5['binmean'],colloz3p5['sbprofile'], 'r-')
	axs[0,1].plot(colns3loz3p5['binmean'],colns3loz3p5['sbprofile'], 'r--')
	axs[0,1].plot(colns1loz3p5['binmean'],colns1loz3p5['sbprofile'], 'r-.')
	axs[0,1].plot(colns01loz3p5['binmean'],colns01loz3p5['sbprofile'], 'r:')
	axs[0,1].plot(coldm7loz3p5['binmean'],coldm7loz3p5['sbprofile'],color='orange',ls=':')
	#axs[0,2] HELLO
	axs[0,2].plot(Hlo['binmean'],Hlo['sbprofileX'], 'g-',alpha=0.5)
	axs[0,2].plot(Hlo['binmean'],Hlo['sbprofileY'], 'g-',alpha=0.5)
	axs[0,2].plot(Hlo['binmean'],Hlo['sbprofileZ'], 'g-',alpha=0.5)
	axs[0,2].plot(Hns3lo['binmean'],Hns3lo['sbprofileX'], 'g--',alpha=0.5)
	axs[0,2].plot(Hns3lo['binmean'],Hns3lo['sbprofileY'], 'g--',alpha=0.5)
	axs[0,2].plot(Hns3lo['binmean'],Hns3lo['sbprofileZ'], 'g--',alpha=0.5)
	axs[0,2].plot(Hns1lo['binmean'],Hns1lo['sbprofileX'], 'g-.',alpha=0.5)
	axs[0,2].plot(Hns1lo['binmean'],Hns1lo['sbprofileY'], 'g-.',alpha=0.5)
	axs[0,2].plot(Hns1lo['binmean'],Hns1lo['sbprofileZ'], 'g-.',alpha=0.5)
	axs[0,2].plot(Hns01lo['binmean'],Hns01lo['sbprofileX'], 'g:',alpha=0.5)	
	axs[0,2].plot(Hns01lo['binmean'],Hns01lo['sbprofileY'], 'g:',alpha=0.5)
	axs[0,2].plot(Hns01lo['binmean'],Hns01lo['sbprofileZ'], 'g:',alpha=0.5)
	#axs[1,0:2] z=3
	#axs[1,0] Dalya
	axs[1,0].plot(CDMloz3['binmean'],CDMloz3['sbprofile'], 'b-')
	axs[1,0].plot(CDMns3loz3['binmean'],CDMns3loz3['sbprofile'], 'b--')
	axs[1,0].plot(CDMns1loz3['binmean'],CDMns1loz3['sbprofile'], 'b-.')
	axs[1,0].plot(CDMns01loz3['binmean'],CDMns01loz3['sbprofile'], 'b:')
	#ax[1,1] COLIBRE
	axs[1,1].plot(colloz3['binmean'],colloz3['sbprofile'], 'r-')
	axs[1,1].plot(colns3loz3['binmean'],colns3loz3['sbprofile'], 'r--')
	axs[1,1].plot(colns1loz3['binmean'],colns1loz3['sbprofile'], 'r-.')
	axs[1,1].plot(colns01loz3['binmean'],colns01loz3['sbprofile'], 'r:')
	axs[1,1].plot(coldm7loz3['binmean'],coldm7loz3['sbprofile'],color='orange',ls=':')
	#axs[1,2] SB(100) as a function of ns
	
	#for legends
	axs[0,0].plot([50,100],[1e-20,1e-20],'k-',label=r'$n_\star = 10$ cm$^{-3}$')
	axs[0,0].plot([50,100],[1e-20,1e-20],'k--',label=r'$n_\star = 3$ cm$^{-3}$')
	axs[0,0].plot([50,100],[1e-20,1e-20],'k-.',label=r'$n_\star = 1$ cm$^{-3}$')
	axs[0,0].plot([50,100],[1e-20,1e-20],'k:',label=r'$n_\star = 0.1$ cm$^{-3}$')
	axs[0,1].plot([50,100],[1e-20,1e-20],ls=':',color='orange',label=r'$M_{frac}$ dust = $10^{-7}$')
	axs[0,0].legend(loc="upper right",prop={'size':12},frameon=False)
	axs[0,1].legend(loc="upper right",prop={'size':12},frameon=False)
	axs[1,1].legend(loc="upper right",prop={'size':12},frameon=False)

	axs[1,2].remove()
	print("Saved figure to "+fout) 
	plt.savefig(fout, bbox_inches='tight')
	
	plt.show()

#PLOT - Effect of nstar+metallicity for skipISM
if checkns and checkmet:
	fout = root+"/SBprofile_dalya_zoom1p5_HELLO_Colibre_m6_CDM_obs_singleobj_NB_Impact_nstar_metal.pdf"
	fig, axs = plt.subplots(1,3,sharey=True,sharex=True,figsize=(10,3))
	plt.subplots_adjust(hspace=0.07)
	for ax in fig.get_axes():
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylim([5e-20*zfac,5e-16*zfac])
		ax.set_xlim([40,400])
		ax.set_xticks([40,50,70,100,200,300,400],labels=[r'$40$',r'$50$',r'70',r'$100$',r'200',r'300',r'$400$'],minor=True)
		ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
		ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
		ax.tick_params(bottom=True)
		ax.set_xlabel(r'R [ckpc]',fontsize=14)
	for ax in fig.get_axes()[0:3]:
		ax.text(43,6e-20*zfac,r'$z=3.5$',fontsize=12)
		#MAGG
		if ax == axs[1]:
			ax.plot(SBMHeII['rA'],SBMHeII['SBA'],'k:',label=r'MAGG w/ HeII')
		else:
			ax.plot(SBMHeII['rA'],SBMHeII['SBA'],'k:')
		ax.fill_between(SBMHeII['rA'],SBMHeII['p25A'],SBMHeII['p75A'],color='black',alpha=0.2)
	fig.supylabel(r'SB$_{\rm intr}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=13)
	axs[0].set_title(r'DaLya',fontsize=14)
	axs[1].set_title(r'COLIBRE',fontsize=14)
	axs[2].set_title(r'HELLO',fontsize=14)
	
	#z = 3.5
	#axs[0] Dalya
	axs[0].plot(CDMloz3p5['binmean'],CDMloz3p5['sbprofile']*zfac, 'b-')
	axs[0].errorbar(CDMloz3p5['binmean'],CDMloz3p5['sbprofile']*zfac,xerr=0.0,yerr=(CDMloz3p5['sbprofile']-0.5*CDMloz3p5['sbprofile'])*zfac,uplims=True,fmt='b')
	axs[0].plot(CDMns3loz3p5['binmean'],CDMns3loz3p5['sbprofile']*zfac, 'b--')
	axs[0].plot(CDMns1loz3p5['binmean'],CDMns1loz3p5['sbprofile']*zfac, 'b-.')
	axs[0].plot(CDMns01loz3p5['binmean'],CDMns01loz3p5['sbprofile']*zfac, 'b:')
	axs[0].plot(CDMmetloz3p5['binmean'],CDMmetloz3p5['sbprofile']*zfac, 'm--')
	#axs[1] COLIBRE
	axs[1].plot(colloz3p5['binmean'],colloz3p5['sbprofile']*zfac, 'r-')
	axs[1].errorbar(colloz3p5['binmean'],colloz3p5['sbprofile']*zfac,xerr=0.0,yerr=(colloz3p5['sbprofile']-0.5*colloz3p5['sbprofile'])*zfac,uplims=True,fmt='r')
	axs[1].plot(colns3loz3p5['binmean'],colns3loz3p5['sbprofile']*zfac, 'r--')
	axs[1].plot(colns1loz3p5['binmean'],colns1loz3p5['sbprofile']*zfac, 'r-.')
	axs[1].plot(colns01loz3p5['binmean'],colns01loz3p5['sbprofile']*zfac, 'r:')
	axs[1].plot(coldm5loz3p5['binmean'],coldm5loz3p5['sbprofile']*zfac,color='orange',ls=':')
	axs[1].plot(colmetloz3p5['binmean'],colmetloz3p5['sbprofile']*zfac,color='magenta',ls='--')
	#axs[2] HELLO
	axs[2].plot(Hlo['binmean'],Hlo['sbprofileX']*zfac, 'g-',alpha=0.5)
	axs[2].errorbar(Hlo['binmean'],Hlo['sbprofileX']*zfac,xerr=0.0,yerr=(Hlo['sbprofileX']-0.5*Hlo['sbprofileX'])*zfac,fmt='g',uplims=True,alpha=0.5)
	axs[2].plot(Hlo['binmean'],Hlo['sbprofileY']*zfac, 'g-',alpha=0.5)
	axs[2].plot(Hlo['binmean'],Hlo['sbprofileZ']*zfac, 'g-',alpha=0.5)
	axs[2].plot(Hns3lo['binmean'],Hns3lo['sbprofileX']*zfac, 'g--',alpha=0.5)
	axs[2].plot(Hns3lo['binmean'],Hns3lo['sbprofileY']*zfac, 'g--',alpha=0.5)
	axs[2].plot(Hns3lo['binmean'],Hns3lo['sbprofileZ']*zfac, 'g--',alpha=0.5)
	axs[2].plot(Hns1lo['binmean'],Hns1lo['sbprofileX']*zfac, 'g-.',alpha=0.5)
	axs[2].plot(Hns1lo['binmean'],Hns1lo['sbprofileY']*zfac, 'g-.',alpha=0.5)
	axs[2].plot(Hns1lo['binmean'],Hns1lo['sbprofileZ']*zfac, 'g-.',alpha=0.5)
	axs[2].plot(Hns01lo['binmean'],Hns01lo['sbprofileX']*zfac, 'g:',alpha=0.5)
	axs[2].plot(Hns01lo['binmean'],Hns01lo['sbprofileY']*zfac, 'g:',alpha=0.5)
	axs[2].plot(Hns01lo['binmean'],Hns01lo['sbprofileZ']*zfac, 'g:',alpha=0.5)
	axs[2].plot(Hmetlo['binmean'],Hmetlo['sbprofileX']*zfac, 'm--',alpha=0.5)
	axs[2].plot(Hmetlo['binmean'],Hmetlo['sbprofileY']*zfac, 'm--',alpha=0.5)
	axs[2].plot(Hmetlo['binmean'],Hmetlo['sbprofileZ']*zfac, 'm--',alpha=0.5)
	
	#for legends
	axs[0].plot([50,100],[1e-20,1e-20],'b-',label=r'$n_\star = 10$ cm$^{-3}$')
	axs[0].plot([50,100],[1e-20,1e-20],'b--',label=r'$n_\star = 3$ cm$^{-3}$')
	axs[0].plot([50,100],[1e-20,1e-20],'b-.',label=r'$n_\star = 1$ cm$^{-3}$')
	axs[0].plot([50,100],[1e-20,1e-20],'b:',label=r'$n_\star = 0.1$ cm$^{-3}$')
	axs[1].plot([50,100],[1e-20,1e-20],ls=':',color='orange',label=r'$M_{frac}$ dust = $10^{-5}$')
	axs[2].plot([50,100],[1e-20,1e-20],ls='--',color='magenta',label=r'$Z \leq 0.05 Z_\odot$')
	axs[0].legend(loc="upper right",prop={'size':10},frameon=False)
	axs[1].legend(loc="upper right",prop={'size':10},frameon=False)
	axs[2].legend(loc="upper right",prop={'size':10},frameon=False)
	
	print("Saved figure to "+fout)
	plt.savefig(fout, bbox_inches='tight')
	
	plt.show()

#PLOT - Effect of nstar for skipISM + Tmin
if checkns and checkT:
	fout2 = root+"/SBprofile_dalya_zoom1p5_HELLO_Colibre_m6_CDM_obs_singleobj_NB_z_"+str(red)+"_n"+noise+"_Impact_nstar_Tmin.pdf"
	fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharey=True,figsize=(17,5))
	plt.subplots_adjust(wspace=0.07)
	for ax in fig.get_axes():
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlabel(r'R [ckpc]',fontsize=14)
		ax.set_xlim([40,400])
		ax.set_xticks([40,50,70,100,200,300,400],labels=[r'$40$',r'$50$',r'70',r'$100$',r'200',r'300',r'$400$'],minor=True)
		ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
		ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
		ax.tick_params(bottom=True)
		ax.set_ylim([1e-19,2e-16])
        	#MAGG
		if ax == ax3:
			ax.plot(SBMz3p5['r'],SBMz3p5['sbprofile'],'k:',label=r'MAGG, z=3.5')
		else:
			ax.plot(SBMz3p5['r'],SBMz3p5['sbprofile'],'k:')
		ax.fill_between(SBMz3p5['r'],SBMz3p5['p25'],SBMz3p5['p75'],color='black',alpha=0.2)
	ax1.set_ylabel(r'SB$_{\rm obs}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=14)
	#ax1 DaLya
	if 'sbprofile' in CDMloz3p5.files:
		ax1.plot(CDMloz3p5['binmean'], CDMloz3p5['sbprofile'], 'b-')#,label=r'$n_\star = 10 \; cm^{-3}$')
		ax1.fill_between(CDMloz3p5['binmean'], CDMloz3p5['p16'], CDMloz3p5['p84'], color='blue',alpha=0.2)
	else:
		ax1.plot(CDMloz3p5['binmean'], CDMloz3p5['sbprofileX'], 'b-',alpha=0.5)#,label=r'$n_\star = 10 \; cm^{-3}$')
		ax1.plot(CDMloz3p5['binmean'], CDMloz3p5['sbprofileY'], 'b-',alpha=0.5)
		ax1.plot(CDMloz3p5['binmean'], CDMloz3p5['sbprofileZ'], 'b-',alpha=0.5)
	if 'sbprofile' in CDMns1loz3p5.files:
		ax1.plot(CDMns1loz3p5['binmean'], CDMns1loz3p5['sbprofile'], 'b--')#,label=r'$n_\star = 1 \; cm^{-3}$')
		ax1.fill_between(CDMns1loz3p5['binmean'], CDMns1loz3p5['p16'], CDMns1loz3p5['p84'], color='blue',alpha=0.2)
	else:
		ax1.plot(CDMns1loz3p5['binmean'], CDMns1loz3p5['sbprofileX'], 'b--',alpha=0.5)#,label=r'$n_\star = 1 \; cm^{-3}$')
		ax1.plot(CDMns1loz3p5['binmean'], CDMns1loz3p5['sbprofileY'], 'b--',alpha=0.5)
		ax1.plot(CDMns1loz3p5['binmean'], CDMns1loz3p5['sbprofileZ'], 'b--',alpha=0.5)
	if 'sbprofile' in CDMTloz3p5.files:
		ax1.plot(CDMTloz3p5['binmean'], CDMTloz3p5['sbprofile'], 'b-.')
	else:
		ax1.plot(CDMTloz3p5['binmean'], CDMTloz3p5['sbprofileX'], 'b-.',alpha=0.5)
		ax1.plot(CDMTloz3p5['binmean'], CDMTloz3p5['sbprofileY'], 'b-.',alpha=0.5)
		ax1.plot(CDMTloz3p5['binmean'], CDMTloz3p5['sbprofileZ'], 'b-.',alpha=0.5)
	#ax2 COLIBRE
	if 'sbprofile' in colloz3p5.files:
		ax2.plot(colloz3p5['binmean'], colloz3p5['sbprofile'], 'r-')#,label=r'$n_\star = 10 \; cm^{-3}$')
		ax2.fill_between(colloz3p5['binmean'], colloz3p5['p16'], colloz3p5['p84'], color='red',alpha=0.2)
	else:
		ax2.plot(colloz3p5['binmean'], colloz3p5['sbprofileX'], 'r-',alpha=0.5)#,label=r'$n_\star = 10 \; cm^{-3}$')
		ax2.plot(colloz3p5['binmean'], colloz3p5['sbprofileY'], 'r-',alpha=0.5)
		ax2.plot(colloz3p5['binmean'], colloz3p5['sbprofileZ'], 'r-',alpha=0.5)
	if 'sbprofile' in colns1loz3p5.files:
		ax2.plot(colns1loz3p5['binmean'], colns1loz3p5['sbprofile'], 'r--')#,label=r'$n_\star = 1 \; cm^{-3}$')
		ax2.fill_between(colns1loz3p5['binmean'], colns1loz3p5['p16'], colns1loz3p5['p84'], color='red',alpha=0.2)
	else:
		ax2.plot(colns1loz3p5['binmean'], colns1loz3p5['sbprofileX'], 'r--',alpha=0.5)#,label=r'$n_\star = 1 \; cm^{-3}$')
		ax2.plot(colns1loz3p5['binmean'], colns1loz3p5['sbprofileY'], 'r--',alpha=0.5)
		ax2.plot(colns1loz3p5['binmean'], colns1loz3p5['sbprofileZ'], 'r--',alpha=0.5)
	if 'sbprofile' in colTloz3p5.files:
		ax2.plot(colTloz3p5['binmean'], colTloz3p5['sbprofile'], 'r-.')
		ax2.fill_between(colTloz3p5['binmean'], colTloz3p5['p16'], colTlo['p84'], color='red',alpha=0.2)
	else:
		ax2.plot(colTloz3p5['binmean'], colTloz3p5['sbprofileX'], 'r-.',alpha=0.5)
		ax2.plot(colTloz3p5['binmean'], colTloz3p5['sbprofileY'], 'r-.',alpha=0.5)
		ax2.plot(colTloz3p5['binmean'], colTloz3p5['sbprofileZ'], 'r-.',alpha=0.5)
	#ax3 HELLO
	if 'sbprofile' in Hlo.files:
		ax3.plot(Hlo['binmean'], Hlo['sbprofile'], 'g--')#,label=r'$n_\star = 10 \; cm^{-3}$')
		ax3.fill_between(Hlo['binmean'], Hlo['p16'], Hlo['p84'], color='green',alpha=0.2)
	else:
		ax3.plot(Hlo['binmean'], Hlo['sbprofileX'], 'g-',alpha=0.5)#,label=r'$n_\star = 10 \; cm^{-3}$')
		ax3.plot(Hlo['binmean'], Hlo['sbprofileY'], 'g-',alpha=0.5)
		ax3.plot(Hlo['binmean'], Hlo['sbprofileZ'], 'g-',alpha=0.5)
	if 'sbprofile' in Hns1lo.files:
		ax3.plot(Hns1lo['binmean'], Hns1lo['sbprofile'], 'g--')#,label=r'$n_\star = 1 \; cm^{-3}$')
		ax3.fill_between(Hns1lo['binmean'], Hns1lo['p16'], Hns1lo['p84'], color='green',alpha=0.2)
	else:
		ax3.plot(Hns1lo['binmean'], Hns1lo['sbprofileX'], 'g--',alpha=0.5)#,label=r'$n_\star = 1 \; cm^{-3}$')
		ax3.plot(Hns1lo['binmean'], Hns1lo['sbprofileY'], 'g--',alpha=0.5)
		ax3.plot(Hns1lo['binmean'], Hns1lo['sbprofileZ'], 'g--',alpha=0.5)
	if 'sbprofile' in HTlo:
		ax3.plot(HTlo['binmean'], HTlo['sbprofile'], 'g-.')
		ax3.fill_between(HTlo['binmean'], HTlo['p16'], HTlo['p84'], color='green',alpha=0.2)
	else:
		ax3.plot(HTlo['binmean'], HTlo['sbprofileX'], 'g-.',alpha=0.5)
		ax3.plot(HTlo['binmean'], HTlo['sbprofileY'], 'g-.',alpha=0.5)
		ax3.plot(HTlo['binmean'], HTlo['sbprofileZ'], 'g-.',alpha=0.5)
	ax2.set_title(r'$z=3.5$, $\log M/M_\odot$ = 12.0 - 12.5',fontsize=14)
	ax1.text(50,2e-19,r'DaLya',fontsize=14)
	ax2.text(50,2e-19,r'COLIBRE',fontsize=14)
	ax3.text(50,2e-19,r'HELLO',fontsize=14)
	#for legends
	ax1.plot([50,100],[1e-20,1e-20],'k-',label=r'$n_\star = 10 \; cm^{-3}$')
	ax1.plot([50,100],[1e-20,1e-20],'k--',label=r'$n_\star = 1 \; cm^{-3}$')
	ax1.plot([50,100],[1e-20,1e-20],'k-.',label=r'$T_{min} = 10^3 K$')
	ax1.legend(loc="upper right",prop={'size':12},frameon=False)
	#ax2.legend(loc="upper right",prop={'size':12},frameon=False)
	ax3.legend(loc="upper right",prop={'size':12},frameon=False)
	print("Saving figure in "+fout2)
	plt.savefig(fout2, bbox_inches='tight')
	plt.show()
'''
#PLOT - Effect of resolution 
if checkzm:
	#Effect of resolution
	fout4 = root+"/SBprofile_cubex_neb_hydro_fire2_lvl7_CDM_z"+str(red)+"_n"+noise+"_Mlow_"+str(mlow)+"_Mhigh_"+str(mhigh)+"_Impact_resolution.pdf"
	plt.figure(4)
	plt.loglog()

	if 'sbprofile' in CDM.files:
		plt.plot(CDM['binmean'], CDM['sbprofile'], 'b-',label=r'$m_b = 7.18 \cdot 10^5 M_\odot$')
		plt.fill_between(CDM['binmean'], CDM['p16'], CDM['p84'], color='blue',alpha=0.2)
	else:
		plt.plot(CDM['binmean'], CDM['sbprofileX'], 'b-',alpha=0.5,label=r'$m_b = 7.18 \cdot 10^5 M_\odot$')
		plt.plot(CDM['binmean'], CDM['sbprofileY'], 'b-',alpha=0.5)
		plt.plot(CDM['binmean'], CDM['sbprofileZ'], 'b-',alpha=0.5)	
	if 'sbprofile' in CDMzm2.files:	
		plt.plot(CDMzm2['binmean'], CDMzm2['sbprofile'], color='teal',ls='-',label=r'$m_b = 8.97 \cdot 10^4 M_\odot$')
		plt.fill_between(CDMzm2['binmean'], CDMzm2['p16'], CDMzm2['p84'], color='teal',alpha=0.2)
	else:
		plt.plot(CDMzm2['binmean'], CDMzm2['sbprofileX'], color='teal',ls='-',alpha=0.5,label=r'$m_b = 8.97 \cdot 10^4 M_\odot$')
		plt.plot(CDMzm2['binmean'], CDMzm2['sbprofileY'], color='teal',ls='-',alpha=0.5)
		plt.plot(CDMzm2['binmean'], CDMzm2['sbprofileZ'], color='teal',ls='-',alpha=0.5)

	#Eventually add ns01
	if checkns:
		if 'sbprofile' in CDMns01.files:
			plt.plot(CDMns01['binmean'],CDMns01['sbprofile'],'b:') 
			plt.fill_between(CDMns01['binmean'],CDMns01['p16'],CDMns01['p84'],color='blue',alpha=0.2)
		else:
			plt.plot(CDMns01['binmean'],CDMns01['sbprofileX'],'b:',alpha=0.5)
			plt.plot(CDMns01['binmean'],CDMns01['sbprofileY'],'b:',alpha=0.5)
			plt.plot(CDMns01['binmean'],CDMns01['sbprofileZ'],'b:',alpha=0.5)

		if 'sbprofile' in CDMzm2ns01.files:
			plt.plot(CDMzm2ns01['binmean'],CDMzm2ns01['sbprofile'],color='teal',ls=':')
			plt.fill_between(CDMzm2ns01['binmean'],CDMzm2ns01['p16'],CDMzm2ns01['p84'],color='teal',alpha=0.2)		
		else:
			plt.plot(CDMzm2ns01['binmean'],CDMzm2ns01['sbprofileX'],color='teal',ls=':',alpha=0.5)
			plt.plot(CDMzm2ns01['binmean'],CDMzm2ns01['sbprofileY'],color='teal',ls=':',alpha=0.5)
			plt.plot(CDMzm2ns01['binmean'],CDMzm2ns01['sbprofileZ'],color='teal',ls=':',alpha=0.5)

	#for labels
	plt.plot([100,110],[1e-22,1e-22],'k-',label=r'$n_\star = 10$ cm$^{-3}$')
	if checkns:
		plt.plot([100,110],[1e-22,1e-22],'k:',label=r'$n_\star = 0.1$ cm$^{-3}$')	

	plt.text(15,5e-19,r'$\log M/M_\odot$ = '+str(milieu),fontsize=12)
	plt.text(15,1e-18,r'$z=$ '+str(red),fontsize=12)
	plt.xlim([CDM['binmean'][0],400])
	plt.ylim([1e-20,2e-16])
	plt.xlabel(r'R [ckpc]',fontsize=14)
	plt.ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=14)
	plt.legend(bbox_to_anchor=(1,1),prop={'size':11},ncols=1,frameon=False)
	print("Saving figure in "+fout4)
	plt.savefig(fout4, bbox_inches='tight')

#PLOT - Effect of AGN (Engine)
if checkAGN:
	fout2 = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_z"+str(red)+"_n"+noise+"_Mlow_"+str(mlow)+"_Mhigh_"+str(mhigh)+"_Impact_AGN.pdf"
	plt.figure(3)
	plt.loglog()

	if mlow == 12.0 and red != 4.0:
		#Engine
		plt.plot(SBRECAL['r'],SBRECAL['sbprofile'],color='orange',ls=':',label=r'Engine')
		plt.fill_between(SBRECAL['r'],SBRECAL['p25'],SBRECAL['p75'],color='orange',alpha=0.2)
		plt.plot(SBNoAGN['r'],SBNoAGN['sbprofile'],'m:',label=r'Engine, no AGN')
		plt.fill_between(SBNoAGN['r'],SBNoAGN['p25'],SBNoAGN['p75'],color='magenta',alpha=0.2)

	#Dalya
	if 'sbprofile' in CDM.files:
		plt.plot(CDM['binmean'], CDM['sbprofile'], 'b-',label=r'DaLya, Fire-2 SN FB x2')
		plt.fill_between(CDM['binmean'], CDM['p16'], CDM['p84'], color='blue',alpha=0.2)
	else:
		plt.plot(CDM['binmean'], CDM['sbprofileX'], 'b-',alpha=0.5,label=r'DaLya, Fire-2 SN FB x2')
		plt.plot(CDM['binmean'], CDM['sbprofileY'], 'b-',alpha=0.5)
		plt.plot(CDM['binmean'], CDM['sbprofileZ'], 'b-',alpha=0.5)

	if 'sbprofile' in CDMfb.files:
		plt.plot(CDMfb['binmean'], CDMfb['sbprofile'], 'g:',label=r'DaLya, Fire-2 SN FB')
		plt.fill_between(CDMfb['binmean'], CDMfb['p16'], CDMfb['p84'], color='green',alpha=0.2)
	else:
		plt.plot(CDMfb['binmean'], CDMfb['sbprofileX'], 'g:',alpha=0.5,label=r'DaLya, Fire-2 SN FB')
		plt.plot(CDMfb['binmean'], CDMfb['sbprofileY'], 'g:',alpha=0.5)
		plt.plot(CDMfb['binmean'], CDMfb['sbprofileZ'], 'g:',alpha=0.5)

	plt.title(r'$\log M/M_\odot$ = '+str(milieu),fontsize=18)
	plt.text(20,1e-19,r'$z=$ '+str(red),fontsize=14)
	plt.xlim([CDM['binmean'][0],200])
	plt.ylim([1e-19,5e-16])
	plt.xlabel(r'R [ckpc]',fontsize=16)
	plt.ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=16)
	plt.legend(loc=0,prop={'size':10})
	plt.savefig(fout2, bbox_inches='tight')
'''
