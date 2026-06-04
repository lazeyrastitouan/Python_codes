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

red   = 4.0 #HELLO has only all sims at z=4 
noise = '5e-19'
root  = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo"

#prepare figure
fout="/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/SBprofile_dalya_zoom1p5_HELLO_CDM_WDM_z"+str(red)+"_n"+noise+".pdf"
fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,sharey=True,figsize=(5,10))
plt.subplots_adjust(hspace=0.05)
for ax in fig.get_axes():
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xticks([40,50,70,100,200,300,400],labels=[r'$40$',r'$50$',r'70',r'$100$',r'200',r'300',r'$400$'],minor=True)
	ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
	ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
	ax.set_xlim([40,400])
	ax.set_ylim([1e-19,6e-16])
	ax.set_ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=14)
ax2.set_xlabel(r'R [ckpc]',fontsize=14)
ax1.text(50,3e-16,r'$z=$ '+str(red),fontsize=12)
ax1.text(50,4e-19,r'$\log M = 12.2$ ',fontsize=12)
ax2.text(50,4e-19,r'$\log M = 12.6$ (DaLya)',fontsize=12)
ax2.text(50,2e-19,r'$\log M = 12.8$ (HELLO)',fontsize=12)
#for legend
ax2.plot([10,20],[1e-20,1e-20],'k-',label=r'DaLya')
ax2.plot([10,20],[1e-20,1e-20],'k--',label=r'HELLO')
ax1.plot([10,20],[1e-20,1e-20],'b-',label=r'CDM')
ax1.plot([10,20],[1e-20,1e-20],'r-',label=r'WDM 1.5keV')
ax1.plot([10,20],[1e-20,1e-20],'g-',label=r'WDM 3.5keV')

#HELLO (only 1 halo)
HLO       = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_singleobj_NB_z_4.022_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz")
HLOWDM1p5 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/WDM_1p5/2cMpc/SBprofile_WDM_1.5_observed_cubex_neb_singleobj_NB_z_4.022_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz")
HLOWDM3p5 = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/WDM_3p5/2cMpc/SBprofile_WDM_3.5_observed_cubex_neb_singleobj_NB_z_4.022_s1_n"+noise+"_Mmin_12.5_Mhigh_13.0.npz")

ax2.plot(HLO['binmean'],HLO['sbprofileX'],'b--',alpha=0.5)
ax2.plot(HLO['binmean'],HLO['sbprofileY'],'b--',alpha=0.5)
ax2.plot(HLO['binmean'],HLO['sbprofileZ'],'b--',alpha=0.5)
ax2.plot(HLOWDM1p5['binmean'],HLOWDM1p5['sbprofileX'],'r--',alpha=0.5)
ax2.plot(HLOWDM1p5['binmean'],HLOWDM1p5['sbprofileY'],'r--',alpha=0.5)
ax2.plot(HLOWDM1p5['binmean'],HLOWDM1p5['sbprofileZ'],'r--',alpha=0.5)
ax2.plot(HLOWDM3p5['binmean'],HLOWDM3p5['sbprofileX'],'g--',alpha=0.5)
ax2.plot(HLOWDM3p5['binmean'],HLOWDM3p5['sbprofileY'],'g--',alpha=0.5)
ax2.plot(HLOWDM3p5['binmean'],HLOWDM3p5['sbprofileZ'],'g--',alpha=0.5)

#DaLya
m = [12.0,12.5]
for i in range(len(m)):
	ax=fig.get_axes()[i]

	CDM    = load(root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(m[i])+"_Mhigh_"+str(m[i]+0.5)+".npz")
	WDM1p5 = load(root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_1p5keV_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(m[i])+"_Mhigh_"+str(m[i]+0.5)+".npz")
	WDM3p5 = load(root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_3p5keV_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(m[i])+"_Mhigh_"+str(m[i]+0.5)+".npz")
	
	if 'sbprofile' in CDM.files:
		ax.plot(CDM['binmean'], CDM['sbprofile'], 'b-')
		ax.fill_between(CDM['binmean'], CDM['p16'], CDM['p84'], color='blue',alpha=0.2)
		ax.plot(WDM3p5['binmean'], WDM3p5['sbprofile'], 'g-')
		ax.fill_between(WDM3p5['binmean'],WDM3p5['p16'],WDM3p5['p84'],color='green',alpha=0.2)
		ax.plot(WDM1p5['binmean'], WDM1p5['sbprofile'], 'r-',marker='')
		ax.fill_between(WDM1p5['binmean'],WDM1p5['p16'],WDM1p5['p84'],color='red',alpha=0.2)
	else:
        	ax.plot(CDM['binmean'], CDM['sbprofileX'], 'b-',alpha=0.5)
        	ax.plot(CDM['binmean'], CDM['sbprofileY'], 'b-',alpha=0.5)
        	ax.plot(CDM['binmean'], CDM['sbprofileZ'], 'b-',alpha=0.5)
        	ax.plot(WDM3p5['binmean'], WDM3p5['sbprofileX'], 'g-',alpha=0.5)
        	ax.plot(WDM3p5['binmean'], WDM3p5['sbprofileY'], 'g-',alpha=0.5)
        	ax.plot(WDM3p5['binmean'], WDM3p5['sbprofileZ'], 'g-',alpha=0.5)
        	ax.plot(WDM1p5['binmean'], WDM1p5['sbprofileX'], 'r-',alpha=0.5)
        	ax.plot(WDM1p5['binmean'], WDM1p5['sbprofileY'], 'r-',alpha=0.5)
        	ax.plot(WDM1p5['binmean'], WDM1p5['sbprofileZ'], 'r-',alpha=0.5)

	ax.legend(loc=0,frameon=False,prop={'size':12})
print("Saving figure to "+fout)
plt.savefig(fout, bbox_inches='tight')
plt.show()
