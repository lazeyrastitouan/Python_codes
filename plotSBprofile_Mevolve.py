from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import argparse
import matplotlib.ticker as mticker

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

noise = '5e-20'

root  = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo"
rootC = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/COLIBRE/L0050N0752"
z     = [4.0,3.5,3.0]

#prepare figure
fout = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/SBprofile_dalya_zoom1p5_HELLO_Colibre_m6_CDM_singleobj_NB_n"+noise+"_Mevolve.pdf"
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
	ax.set_ylim([1e-19,9e-16])
ax1.set_ylabel(r'SB$_{\rm obs}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=14)

for red in z:
	i  = z.index(red)	
	ax = fig.get_axes()[i]
	ax.text(50,4e-19,r'$z=$ '+str(red),fontsize=14)

	#Dalya
	fCDM1 = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.0_Mhigh_12.5.npz"
	CDM1  = load(fCDM1)

	if red == 3.0:
		#DaLya
		fCDM2 = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.5_Mhigh_13.1.npz"
		CDM2  = load(fCDM2)
		if 'sbprofile' in CDM1.files:
			ax.plot(CDM1['binmean'], CDM1['sbprofile'], color='cyan', ls='-',marker='',label=r'DaLya, $\log M = 12.3$')
			ax.fill_between(CDM1['binmean'], CDM1['p16'], CDM1['p84'], color='cyan',alpha=0.2)
		else:
			ax.plot(CDM1['binmean'], CDM1['sbprofileX'], color='cyan', ls='-',marker='',alpha=0.5,label=r'DaLya, $\log M = 12.3$')
			ax.plot(CDM1['binmean'], CDM1['sbprofileY'], color='cyan', ls='-',marker='',alpha=0.5)
			ax.plot(CDM1['binmean'], CDM1['sbprofileZ'], color='cyan', ls='-',marker='',alpha=0.5)
		if 'sbprofile' in CDM2.files:
			ax.plot(CDM2['binmean'], CDM2['sbprofile'], 'b-',label=r'DaLya, $\log M = 12.7$')
			ax.fill_between(CDM2['binmean'], CDM2['p16'], CDM2['p84'], color='blue',alpha=0.2)
		else:
			ax.plot(CDM2['binmean'], CDM2['sbprofileX'], 'b-',alpha=0.5,label=r'DaLya, $\log M = 12.7$')
			ax.plot(CDM2['binmean'], CDM2['sbprofileY'], 'b-',alpha=0.5)
			ax.plot(CDM2['binmean'], CDM2['sbprofileZ'], 'b-',alpha=0.5)

		#Colibre
		fC1 = rootC+"/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.2_Mhigh_12.5.npz"
		C1  = load(fC1)
		fC2 = rootC+"/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz"
		C2  = load(fC2)

		if 'sbprofile' in C1.files:
			ax.plot(C1['binmean'], C1['sbprofile'], color='orange', ls='-',marker='',label=r'Colibre, $\log M = 12.3$')
			ax.fill_between(C1['binmean'], C1['p16'], C1['p84'], color='orange',alpha=0.2)
		else:
			ax.plot(C1['binmean'], C1['sbprofileX'], color='orange', ls='-',marker='',alpha=0.5,label=r'Colibre, $\log M = 12.3$')
			ax.plot(C1['binmean'], C1['sbprofileY'], color='orange', ls='-',marker='',alpha=0.5)
			ax.plot(C1['binmean'], C1['sbprofileZ'], color='orange', ls='-',marker='',alpha=0.5)
		if 'sbprofile' in C2.files:
			ax.plot(C2['binmean'], C2['sbprofile'], color='red', ls='-',marker='',label=r'Colibre, $\log M = 12.6$')
			ax.fill_between(C2['binmean'], C2['p16'], C2['p84'], color='red',alpha=0.2)
		else:
			ax.plot(C1['binmean'], C2['sbprofileX'], color='red', ls='-',marker='',alpha=0.5,label=r'Colibre, $\log M = 12.6$')			
			ax.plot(C2['binmean'], C2['sbprofileY'], color='red', ls='-',marker='',alpha=0.5)
			ax.plot(C2['binmean'], C2['sbprofileZ'], color='red', ls='-',marker='',alpha=0.5)
		
		#Borisova and MUSEUM
		#SBb = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Borisova_SBprofile_observed_fromMariono_z_3.0.npz") 
		#SBM = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MUSEUM_FAB_SBprofile_observed_z_3.2.npz")
		#ax.plot(SBM['r'],SBM['sbprofile'],color='grey',ls=':',marker='',label=r'MUSEUM FAB, z = 3.2')
		#ax.fill_between(SBM['r'],SBM['p25'],SBM['p75'],color='grey',alpha=0.2)  
		#ax.plot(SBb['r'],SBb['sbprofile'],'k--',label=r'Borisova+ 2016, z = 3.0') 	
		#ax.fill_between(SBb['r'],SBb['p25'],SBb['p75'],color='black',alpha=0.2)

	elif red ==3.5:
		#dalya 
		fCDM2 = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz"
		CDM2  = load(fCDM2)
		if 'sbprofile' in CDM1.files:
			ax.plot(CDM1['binmean'], CDM1['sbprofile'], color='cyan', ls='-',marker='',label=r'DaLya, $\log M = 12.3$')
			ax.fill_between(CDM1['binmean'], CDM1['p16'], CDM1['p84'], color='cyan',alpha=0.2)
		else:
			ax.plot(CDM1['binmean'], CDM1['sbprofileX'], color='cyan', ls='-',marker='',alpha=0.5,label=r'DaLya, $\log M = 12.3$')
			ax.plot(CDM1['binmean'], CDM1['sbprofileY'], color='cyan', ls='-',marker='',alpha=0.5)
			ax.plot(CDM1['binmean'], CDM1['sbprofileZ'], color='cyan', ls='-',marker='',alpha=0.5)
		if 'sbprofile' in CDM2.files:
			ax.plot(CDM2['binmean'], CDM2['sbprofile'], 'b-',label=r'DaLya, $\log M = 12.6$')
			ax.fill_between(CDM2['binmean'], CDM2['p16'], CDM2['p84'], color='blue',alpha=0.2)
		else:   
			ax.plot(CDM2['binmean'], CDM2['sbprofileX'], 'b-',alpha=0.5,label=r'DaLya, $\log M = 12.6$')
			ax.plot(CDM2['binmean'], CDM2['sbprofileY'], 'b-',alpha=0.5)
			ax.plot(CDM2['binmean'], CDM2['sbprofileZ'], 'b-',alpha=0.5)
	
		#Colibre
		fC1 = rootC+"/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.0_Mhigh_12.5.npz"
		C1  = load(fC1)
		fC2 = rootC+"/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz"	
		C2  = load(fC2)

		if 'sbprofile' in C1.files:
			ax.plot(C1['binmean'], C1['sbprofile'], color='orange', ls='-',marker='',label=r'Colibre, $\log M = 12.1$')
			ax.fill_between(C1['binmean'], C1['p16'], C1['p84'], color='orange',alpha=0.2)
		else:
			ax.plot(C1['binmean'], C1['sbprofileX'], color='orange', ls='-',marker='',alpha=0.5,label=r'Colibre, $\log M = 12.1$')
			ax.plot(C1['binmean'], C1['sbprofileY'], color='orange', ls='-',marker='',alpha=0.5)
			ax.plot(C1['binmean'], C1['sbprofileZ'], color='orange', ls='-',marker='',alpha=0.5)
		if 'sbprofile' in C2.files:
			ax.plot(C2['binmean'], C2['sbprofile'], color='red', ls='-',marker='',label=r'Colibre, $\log M = 12.5$')
			ax.fill_between(C2['binmean'], C2['p16'], C2['p84'], color='red',alpha=0.2)
		else:
			ax.plot(C2['binmean'], C2['sbprofileX'], color='red', ls='-',marker='',alpha=0.5,label=r'Colibre, $\log M = 12.5$')
			ax.plot(C2['binmean'], C2['sbprofileY'], color='red', ls='-',marker='',alpha=0.5)
			ax.plot(C2['binmean'], C2['sbprofileZ'], color='red', ls='-',marker='',alpha=0.5)

		# HELLO
		HLO1 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n"+noise+"_Mmin_12.0_Mhigh_12.5.npz")
		HLO2 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_3.596_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz")	
		if 'sbprofile' in HLO1.files:		
			ax.plot(HLO1['binmean'],HLO1['sbprofile'],color='lime',ls='-',marker='',label=r'HELLO, $\log M = 12.0$')
			ax.fill_between(HLO1['binmean'],HLO1['p16'],HLO1['p84'],color='lime',alpha=0.2)
		else:
			ax.plot(HLO1['binmean'],HLO1['sbprofileX'],color='lime',ls='-',marker='',alpha=0.5,label=r'HELLO, $\log M = 12.0$')
			ax.plot(HLO1['binmean'],HLO1['sbprofileY'],color='lime',ls='-',marker='',alpha=0.5)
			ax.plot(HLO1['binmean'],HLO1['sbprofileZ'],color='lime',ls='-',marker='',alpha=0.5)
		if 'sbprofile' in HLO2.files:
			ax.plot(HLO2['binmean'],HLO2['sbprofile'],'g-',label=r'HELLO, $\log M = 12.9$')
			ax.fill_between(HLO1['binmean'],HLO2['p16'],HLO2['p84'],color='green',alpha=0.2)
		else:
			ax.plot(HLO2['binmean'],HLO2['sbprofileX'],'g-',alpha=0.5,label=r'HELLO, $\log M = 12.9$')
			ax.plot(HLO2['binmean'],HLO2['sbprofileY'],'g-',alpha=0.5)
			ax.plot(HLO2['binmean'],HLO2['sbprofileZ'],'g-',alpha=0.5)

		#MAGG
		#SBM = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_3.5.npz")
		#ax.plot(SBM['r'],SBM['sbprofile'],'k:',label=r'MAGG, z = 3.5')
		#ax.fill_between(SBM['r'],SBM['p25'],SBM['p75'],color='black',alpha=0.2)
		#ax1.plot(SBb['r'],SBb['sbpezzulli'],'k--',label=r'Borisova+ 2016, z=3.4')

	elif red == 4.0:
		#dalya
		fCDM2 = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz"
		CDM2  = load(fCDM2)
		if 'sbprofile' in CDM1.files:
			ax.plot(CDM1['binmean'], CDM1['sbprofile'], color='cyan', ls='-',marker='',label=r'DaLya, $\log M = 12.2$')
			ax.fill_between(CDM1['binmean'], CDM1['p16'], CDM1['p84'], color='cyan',alpha=0.2)
		else:
			ax.plot(CDM1['binmean'], CDM1['sbprofileX'], color='cyan', ls='-',marker='',alpha=0.5,label=r'DaLya, $\log M = 12.2$')
			ax.plot(CDM1['binmean'], CDM1['sbprofileY'], color='cyan', ls='-',marker='',alpha=0.5)
			ax.plot(CDM1['binmean'], CDM1['sbprofileZ'], color='cyan', ls='-',marker='',alpha=0.5)
		if 'sbprofile' in CDM2.files:
			ax.plot(CDM2['binmean'], CDM2['sbprofile'], 'b-',label=r'DaLya, $\log M = 12.6$')
			ax.fill_between(CDM['binmean'], CDM2['p16'], CDM['p84'], color='blue',alpha=0.2)
		else:
			ax.plot(CDM2['binmean'], CDM2['sbprofileX'], 'b-',alpha=0.5,label=r'DaLya, $\log M = 12.6$')
			ax.plot(CDM2['binmean'], CDM2['sbprofileY'], 'b-',alpha=0.5)
			ax.plot(CDM2['binmean'], CDM2['sbprofileZ'], 'b-',alpha=0.5)
		
		#Colibre
		fC1 = rootC+"/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_12.0_Mhigh_12.5.npz"
		C1  = load(fC1)
		
		if 'sbprofile' in C1.files:
			ax.plot(C1['binmean'], C1['sbprofile'], color='orange', ls='-',marker='',label=r'Colibre, $\log M = 12.1$')
			ax.fill_between(C1['binmean'], C1['p16'], C1['p84'], color='orange',alpha=0.2)
		else:
			ax.plot(C1['binmean'], C1['sbprofileX'], color='orange', ls='-',marker='',alpha=0.5,label=r'Colibre, $\log M = 12.1$')
			ax.plot(C1['binmean'], C1['sbprofileY'], color='orange', ls='-',marker='',alpha=0.5)
			ax.plot(C1['binmean'], C1['sbprofileZ'], color='orange', ls='-',marker='',alpha=0.5)

		# HELLO
		HLO1 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_4.022_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz")
		if 'sbprofile' in HLO1.files:
			ax.plot(HLO1['binmean'],HLO1['sbprofile'],'g-',label=r'HELLO, $\log M = 12.8$')
			ax.fill_between(HLO1['binmean'],HLO1['p16'],HLO1['p84'],color='green',alpha=0.2)
		else:
			ax.plot(HLO1['binmean'],HLO1['sbprofileX'],'g-',alpha=0.5,label=r'HELLO, $\log M = 12.8$')
			ax.plot(HLO1['binmean'],HLO1['sbprofileY'],'g-',alpha=0.5)
			ax.plot(HLO1['binmean'],HLO1['sbprofileZ'],'g-',alpha=0.5)		

		#MAGG
		#SBM = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_4.1.npz")
		#ax.plot(SBM['r'],SBM['sbprofile'],'k:',label=r'MAGG, z = 4.1')
		#ax.fill_between(SBM['r'],SBM['p25'],SBM['p75'],color='black',alpha=0.2)


	ax.legend(loc="upper right",prop={'size':12},frameon=False)

print("Saving figure to "+fout)
#plt.savefig(fout, bbox_inches='tight')
plt.show()
