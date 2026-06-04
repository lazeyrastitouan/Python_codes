from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import argparse

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# input parameters
parser  = argparse.ArgumentParser(description='Plot SB profiles for 3 DM for DaLya and HELLO. Args: < DaLya Reshift (float) >, < Mlow (float) ')
parser.add_argument('red',nargs=1,type=float,help='Redshift (float)')
parser.add_argument('mlow',nargs=1,type=float,help='Log 10 of lower mass bin (float)')
parser.add_argument('noise',nargs=1,type=str,help='MockOb noise (0 or 5e-20, str)')
args    = parser.parse_args()
red     = args.red[0]
mlow    = args.mlow[0]
noise   = args.noise[0]
mhigh   = mlow+0.5
milieu  = (mlow+mhigh)/2.
if red ==3.5:
	redH = 3.596
	snap = '128'
elif red == 4.0:
	redH = 4.022
	snap = '112'
else:
	exit("No valid HELLO redshift!")

checkns  = False

# ---
#LOAD
#DaLya
root    = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo"
fCDM    = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"
#fWDM1p5 = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_1p5keV_observed_z_"+str(red)+"_s1_n0_t0_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"
#fWDM3p5 = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_WDM_3p5keV_observed_z_"+str(red)+"_s1_n0_t0_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"

CDM    = load(fCDM)
#WDM1p5 = load(fWDM1p5)
#WDM3p5 = load(fWDM3p5)

#Hello
rootH    = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO"
fCDMH    = rootH+"/CDM/2cMpc/SBprofile_CDM_observed_cubex_neb_z_"+str(redH)+"_s1_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"
#fWDM3p5H = rootH+"/WDM_3p5/2cMpc/SBprofile_observed_s1_n0_t0_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+"_00"+snap+"_z"+str(redH)+".npz" 
#fFDMH    = rootH+"/FDM/2cMpc/SBprofile_observed_s1_n0_t0_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+"_00"+snap+"_z"+str(redH)+".npz"

CDMH    = load(fCDMH)
#WDM3p5H = load(fWDM3p5H)
#FDMH    = load(fFDMH)

#if red == 4.0:
#	fWDM1p5H = rootH+"/WDM_1p5/2cMpc/SBprofile_observed_s1_n0_t0_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+"_00"+snap+"_z"+str(redH)+".npz"
#	WDM1p5H  = load(fWDM1p5H)

if checkns:
	fCDMns01  = root+"/nstar01/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"
	CDMns01   = load(fCDMns01)
	fCDMHns01 = rootH+"/CDM/2cMpc/nstar01/SBprofile_CDM_observed_cubex_neb_z_"+str(redH)+"_s1_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"
	CDMHns01  = load(fCDMHns01)

#Borisova/MAGG
if red == 3.5:
        SBb = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Borisova_SBprofile_exp_power_observed_z_3.36.npz")
        SBM = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_3.5.npz")
elif red == 4.0:
	SBM = load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/MAGG_SBprofile_observed_z_4.1.npz")

# ---
#PLOT
fout = rootH+"/SBprofile_HELLO_zoom1p5_hydro_fire2_lvl7_CDM_z"+str(red)+"_Mlow_"+str(mlow)+"_Mhigh_"+str(mhigh)+".pdf"
print("Saving figure to "+fout)
plt.figure(1)
plt.loglog()
#plt.title(r'$\log M = \, $'+str(mlow),fontsize=18)

if red == 4.0:
	plt.plot(SBM['r'],SBM['sbprofile'],'k:',label=r'MAGG, z=4.1')
	plt.fill_between(SBM['r'],SBM['p25'],SBM['p75'],color='black',alpha=0.2)
elif red == 3.5:
	plt.plot(SBM['r'],SBM['sbprofile'],'k:',label=r'MAGG, z=3.5')
	plt.fill_between(SBM['r'],SBM['p25'],SBM['p75'],color='black',alpha=0.2)
	plt.plot(SBb['r'],SBb['sbpezzulli'],'k--',label=r'Borisova+ 2016, z=3.2') #should be 3.4

#DaLya
#Treat particularly mass bins with only 1 halo
if 'sbprofile' not in CDM.files:
	plt.plot(CDM['binmean'], CDM['sbprofileX'], 'b-',alpha=0.5,label=r'DaLya')
	plt.plot(CDM['binmean'], CDM['sbprofileY'], 'b-',alpha=0.5)
	plt.plot(CDM['binmean'], CDM['sbprofileZ'], 'b-',alpha=0.5)
#	plt.plot(WDM3p5['binmean'], WDM3p5['sbprofileX'], 'g-',alpha=0.5,label=r'WDM 3.5keV')
#	plt.plot(WDM3p5['binmean'], WDM3p5['sbprofileY'], 'g-',alpha=0.5)
#	plt.plot(WDM3p5['binmean'], WDM3p5['sbprofileZ'], 'g-',alpha=0.5)
#	plt.plot(WDM1p5['binmean'], WDM1p5['sbprofileX'], 'r-',alpha=0.5,label=r'WDM 1.5keV')
#	plt.plot(WDM1p5['binmean'], WDM1p5['sbprofileY'], 'r-',alpha=0.5)
#	plt.plot(WDM1p5['binmean'], WDM1p5['sbprofileZ'], 'r-',alpha=0.5)
else:
	plt.plot(CDM['binmean'], CDM['sbprofile'], 'b-',label=r'DaLya')
	plt.fill_between(CDM['binmean'], CDM['p16'], CDM['p84'], color='blue',alpha=0.2)
#	plt.plot(WDM3p5['binmean'], WDM3p5['sbprofile'], 'g-',label=r'WDM 3.5keV')
#	plt.fill_between(WDM3p5['binmean'],WDM3p5['p16'],WDM3p5['p84'],color='green',alpha=0.2)
#	plt.plot(WDM1p5['binmean'], WDM1p5['sbprofile'], 'r-',label=r'WDM 1.5keV')
#	plt.fill_between(WDM1p5['binmean'],WDM1p5['p16'],WDM1p5['p84'],color='red',alpha=0.2)

#HELLO
if 'sbprofile' in CDMH.files:
	plt.plot(CDMH['binmean'], CDMH['sbprofile'], 'b--', label=r'HELLO')
	plt.fill_between(CDMH['binmean'], CDMH['p16'], CDMH['p84'], color='blue',alpha=0.2)
else:
	plt.plot(CDMH['binmean'], CDMH['sbprofileX'], 'b--',alpha=0.5, label=r'HELLO')
	plt.plot(CDMH['binmean'], CDMH['sbprofileY'], 'b--',alpha=0.5)
	plt.plot(CDMH['binmean'], CDMH['sbprofileZ'], 'b--',alpha=0.5)
#if 'sbprofile' in WDM3p5H.files:
#        plt.plot(WDM3p5H['binmean'], WDM3p5H['sbprofile'], 'g--')
#        plt.fill_between(WDM3p5H['binmean'], WDM3p5H['p16'], WDM3p5H['p84'], color='green',alpha=0.2)
#else:
#        plt.plot(WDM3p5H['binmean'], WDM3p5H['sbprofileX'], 'g--',alpha=0.5)
#        plt.plot(WDM3p5H['binmean'], WDM3p5H['sbprofileY'], 'g--',alpha=0.5)
#        plt.plot(WDM3p5H['binmean'], WDM3p5H['sbprofileZ'], 'g--',alpha=0.5)
#if 'sbprofile' in FDMH.files:
#        plt.plot(FDMH['binmean'], FDMH['sbprofile'], 'm--',label=r'FDM $2\cdot10^{-22}$ eV/$c^2$')
#        plt.fill_between(FDMH['binmean'], FDMH['p16'], FDMH['p84'], color='magenta',alpha=0.2)
#else:
#        plt.plot(FDMH['binmean'], FDMH['sbprofileX'], 'm--',alpha=0.5,label=r'FDM $2\cdot10^{-22}$ eV/$c^2$')
#        plt.plot(FDMH['binmean'], FDMH['sbprofileY'], 'm--',alpha=0.5)
#        plt.plot(FDMH['binmean'], FDMH['sbprofileZ'], 'm--',alpha=0.5)
#if red == 4.0:
#	if 'sbprofile' in WDM1p5H.files:
#        	plt.plot(WDM1p5H['binmean'], WDM1p5H['sbprofile'], 'r--')
#        	plt.fill_between(WDM1p5H['binmean'], WDM1p5H['p16'], WDM1p5H['p84'], color='red',alpha=0.2)
#	else:
#        	plt.plot(WDM1p5H['binmean'], WDM1p5H['sbprofileX'], 'r--',alpha=0.5)
#        	plt.plot(WDM1p5H['binmean'], WDM1p5H['sbprofileY'], 'r--',alpha=0.5)
#        	plt.plot(WDM1p5H['binmean'], WDM1p5H['sbprofileZ'], 'r--',alpha=0.5)

#for legend:
#plt.plot([0,1],[1e-23,1e-23],'k-',label=r'DaLya')# $m_b \, = \, 7.18 \cdot 10^5 \; M_\odot$')
#plt.plot([0,1],[1e-23,1e-23],'k--',label=r'HELLO')
plt.text(30,1e-18,r'$z=$ '+str(red),fontsize=14)
if mlow == 12.0:
	plt.text(30,6e-19,r'$\log M = \, $'+str(12.3),fontsize=14)
elif mlow == 12.5:
	plt.text(30,6e-19,r'$\log M = \, $'+str(12.7),fontsize=14)
plt.xlim([20,400])
plt.ylim([1e-19,1e-15])
plt.xlabel(r'R [ckpc]',fontsize=14)
plt.ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=14)
plt.legend(frameon=False,loc=0,prop={'size':12})
plt.savefig(fout, bbox_inches='tight')
#plt.close()

if checkns:
	plt.figure(2)
	fout2 = rootH+"/SBprofile_HELLO_zoom1p5_hydro_fire2_lvl7_CDM_z"+str(red)+"_Mlow_"+str(mlow)+"_Mhigh_"+str(mhigh)+"_Impact_nstar.pdf"
	print("Saving figure to "+fout2)
	plt.figure(1)
	plt.loglog()
	plt.title(r'$\log M = \, $'+str(mlow),fontsize=18)

	if mlow == 12.0:
		plt.plot(SBM['r'],SBM['sbprofile'],'k:',label=r'MAGG, z=3.5')
		plt.fill_between(SBM['r'],SBM['p25'],SBM['p75'],color='black',alpha=0.2)
		plt.plot(SBb['r'],SBb['sbpezzulli'],'k--',label=r'Borisova+ 2016, z=3.4')

	#DaLya
	#Treat particularly mass bins with only 1 halo
	if 'sbprofile' not in CDM.files:
        	plt.plot(CDM['binmean'], CDM['sbprofileX'], 'b-',alpha=0.5,label=r'$n_\star = 10$ cm$^{-3}$')
        	plt.plot(CDM['binmean'], CDM['sbprofileY'], 'b-',alpha=0.5)
        	plt.plot(CDM['binmean'], CDM['sbprofileZ'], 'b-',alpha=0.5)
	else:
        	plt.plot(CDM['binmean'], CDM['sbprofile'], 'b-',label=r'$n_\star = 10$ cm$^{-3}$')
        	plt.fill_between(CDM['binmean'], CDM['p16'], CDM['p84'], color='blue',alpha=0.2)
	if 'sbprofile' not in CDMns01.files:
		plt.plot(CDMns01['binmean'], CDMns01['sbprofileX'], 'b:',alpha=0.5,label=r'$n_\star = 0.1$ cm$^{-3}$')
		plt.plot(CDMns01['binmean'], CDMns01['sbprofileY'], 'b:',alpha=0.5)
		plt.plot(CDMns01['binmean'], CDMns01['sbprofileZ'], 'b:',alpha=0.5)
	else:
		plt.plot(CDMns01['binmean'], CDMns01['sbprofile'], 'b:',label=r'$n_\star = 0.1$ cm$^{-3}$')
		plt.fill_between(CDMns01['binmean'], CDMns01['p16'], CDMns01['p84'], color='blue',alpha=0.2)

	#HELLO
	if 'sbprofile' in CDMH.files:
		plt.plot(CDMH['binmean'], CDMH['sbprofile'], color='orange', ls='-')
		plt.fill_between(CDMH['binmean'], CDMH['p16'], CDMH['p84'], color='blue',alpha=0.2)
	else:
		plt.plot(CDMH['binmean'], CDMH['sbprofileX'], color='orange',ls='-',alpha=0.5)
		plt.plot(CDMH['binmean'], CDMH['sbprofileY'], color='orange',ls='-',alpha=0.5)
		plt.plot(CDMH['binmean'], CDMH['sbprofileZ'], color='orange',ls='-',alpha=0.5)
	if 'sbprofile' in CDMHns01.files:
		plt.plot(CDMHns01['binmean'], CDMHns01['sbprofile'], color='orange', ls=':')
		plt.fill_between(CDMHns01['binmean'], CDMHns01['p16'], CDMHns01['p84'], color='orange',alpha=0.2)
	else:
		plt.plot(CDMHns01['binmean'], CDMHns01['sbprofileX'], color='orange',ls=':',alpha=0.5)
		plt.plot(CDMHns01['binmean'], CDMHns01['sbprofileY'], color='orange',ls=':',alpha=0.5)
		plt.plot(CDMHns01['binmean'], CDMHns01['sbprofileZ'], color='orange',ls=':',alpha=0.5)

	#for legend:
	plt.plot([0,1],[1e-23,1e-23],'b-',label=r'DaLya')# $m_b \, = \, 7.18 \cdot 10^5 \; M_\odot$')
	plt.plot([0,1],[1e-23,1e-23],color='orange',ls='-',label=r'HELLO')
	plt.text(20,1e-18,r'$z=$ '+str(red),fontsize=12)
	plt.xlim([CDM['binmean'][0],400])
	plt.ylim([1e-20,5e-15])
	plt.xlabel(r'R [ckpc]',fontsize=16)
	plt.ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=16)
	plt.legend(frameon=False,loc=0,prop={'size':12})
	plt.savefig(fout2, bbox_inches='tight')

plt.show()
