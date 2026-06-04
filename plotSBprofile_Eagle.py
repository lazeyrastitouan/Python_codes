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
parser  = argparse.ArgumentParser(description='Plot SB profiles for 3 DM for DaLya and Eagle. Args: < DaLya Reshift (float) >, < Mlow (float) ')
parser.add_argument('red',nargs=1,type=float,help='Redshift (float)')
parser.add_argument('mlow',nargs=1,type=float,help='Log 10 of lower mass bin (float)')
parser.add_argument('noise',nargs=1,type=str,help='MockOb noise (0 or 5e-20, str)')
args    = parser.parse_args()
red     = args.red[0]
mlow    = args.mlow[0]
noise   = args.noise[0]
mhigh   = mlow+0.6
milieu  = (mlow+mhigh)/2.
if red ==3.0:
	redE = 3.017
else:
	exit("No valid Eagle redshift!")

# ---
#LOAD
#DaLya
root     = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo"
fCDM     = root+"/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"
CDM      = load(fCDM)
fCDMns01 = root+"/nstar01/SBprofile_zoom1p5_hydro_fire2_lvl7_CDM_observed_cubex_neb_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhigh)+".npz"
CDMns01  = load(fCDMns01)

#Eagle
rootE = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Eagle/RefL0100N1504/z3/"
fCDME = rootE+"/SBprofile_Eagle_observed_cubex_neb_z_"+str(redE)+"_s1_n"+noise+"_Mmin_12.75_Mhigh_13.25.npz"
CDME  = load(fCDME)

# ---
#PLOT
fout = rootE+"/SBprofile_Eagle_zoom1p5_hydro_fire2_lvl7_CDM_z"+str(red)+"_M12p7.pdf"
print("Saving figure to "+fout)
plt.figure(1)
plt.loglog()
plt.title(r'$\log M = \, $'+str(mlow),fontsize=18)
#DaLya
#Treat particularly mass bins with only 1 halo
if 'sbprofile' not in CDM.files:
	plt.plot(CDM['binmean'], CDM['sbprofileX'], 'b--',alpha=0.5,label=r'DaLya, $n_\star=10$ cm$^{-3}$')
	plt.plot(CDM['binmean'], CDM['sbprofileY'], 'b--',alpha=0.5)
	plt.plot(CDM['binmean'], CDM['sbprofileZ'], 'b--',alpha=0.5)
else:
	plt.plot(CDM['binmean'], CDM['sbprofile'], 'b--',label=r'DaLya, $n_\star=10$ cm$^{-3}$')
	plt.fill_between(CDM['binmean'], CDM['p16'], CDM['p84'], color='blue',alpha=0.2)
if 'sbprofile' not in CDM.files:
        plt.plot(CDMns01['binmean'], CDMns01['sbprofileX'], 'b-',alpha=0.5,label=r'DaLya, $n_\star=0.1$ cm$^{-3}$')
        plt.plot(CDMns01['binmean'], CDMns01['sbprofileY'], 'b-',alpha=0.5)
        plt.plot(CDMns01['binmean'], CDMns01['sbprofileZ'], 'b-',alpha=0.5)
else:
        plt.plot(CDMns01['binmean'], CDMns01['sbprofile'], 'b-',label=r'DaLya, $n_\star=0.1$ cm$^{-3}$')
        plt.fill_between(CDMns01['binmean'], CDMns01['p16'], CDMns01['p84'], color='blue',alpha=0.2)
#Eagle
if 'sbprofile' in CDME.files:
	plt.plot(CDME['binmean'], CDME['sbprofile'], 'r-',label=r'Eagle, $n_\star=0.1$ cm$^{-3}$')
	plt.fill_between(CDME['binmean'], CDME['p16'], CDME['p84'], color='blue',alpha=0.2)
else:
	plt.plot(CDME['binmean'], CDME['sbprofileX'], 'r-',alpha=0.5, label=r'Eagle, $n_\star=0.1$ cm$^{-3}$')
	plt.plot(CDME['binmean'], CDME['sbprofileY'], 'r-',alpha=0.5)
	plt.plot(CDME['binmean'], CDME['sbprofileZ'], 'r-',alpha=0.5)

#for legend:
plt.text(20,1e-18,r'$z=$ '+str(red),fontsize=12)
plt.xlim([CDM['binmean'][0],400])
plt.ylim([1e-20,1e-15])
plt.xlabel(r'R [ckpc]',fontsize=16)
plt.ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=16)
plt.legend(frameon=False,loc=0,prop={'size':12})
plt.savefig(fout, bbox_inches='tight')
plt.show()
