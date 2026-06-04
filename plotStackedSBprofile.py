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
parser  = argparse.ArgumentParser(description='Stack 3 LOS SB profiles and plot + ratios. Args: < SB threshold (float) >')
parser.add_argument('threshold',nargs=1,type=str, help='SB threshold of profiles to stack')
args   = parser.parse_args()
thresh = args.threshold[0]

root   = '/data/titouan/Eagle_cosmo_gizmo/200_512/'
folder = '/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/4cMpc/'
rootlo = root+'/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low/snapdir_002/4cMpc/'

# ---

#LOAD
#GIZMO 
CDM1X         = load(root+'/CDM/'+folder+'mockCubeX_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
CDM1Y         = load(root+'/CDM/'+folder+'mockCubeY_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
CDM1Z         = load(root+'/CDM/'+folder+'mockCubeZ_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM1p51X      = load(root+'/WDM_1p5keV/'+folder+'mockCubeX_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM1p51Y      = load(root+'/WDM_1p5keV/'+folder+'mockCubeY_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM1p51Z      = load(root+'/WDM_1p5keV/'+folder+'mockCubeZ_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM3p51X      = load(root+'/WDM_3p5keV/'+folder+'mockCubeX_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM3p51Y      = load(root+'/WDM_3p5keV/'+folder+'mockCubeY_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM3p51Z      = load(root+'/WDM_3p5keV/'+folder+'mockCubeZ_NoNoise_gn1_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
CDM2X         = load(root+'/CDM/'+folder+'mockCubeX_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
CDM2Y         = load(root+'/CDM/'+folder+'mockCubeY_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
CDM2Z         = load(root+'/CDM/'+folder+'mockCubeZ_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM1p52X      = load(root+'/WDM_1p5keV/'+folder+'mockCubeX_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM1p52Y      = load(root+'/WDM_1p5keV/'+folder+'mockCubeY_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM1p52Z      = load(root+'/WDM_1p5keV/'+folder+'mockCubeZ_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM3p52X      = load(root+'/WDM_3p5keV/'+folder+'mockCubeX_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM3p52Y      = load(root+'/WDM_3p5keV/'+folder+'mockCubeY_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
WDM3p52Z      = load(root+'/WDM_3p5keV/'+folder+'mockCubeZ_NoNoise_gn2_s1_n0_SB_observed_SBprofile_t'+thresh+'.npz')
#CDMfbloX     = load(rootlo+'mockCubeX_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#CDMfbloY     = load(rootlo+'mockCubeY_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#CDMfbloZ     = load(rootlo+'mockCubeZ_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#CDMhiresX    = load(root+'/CDM/zoom2/mockCubeX_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#CDMhiresY    = load(root+'/CDM/zoom2/mockCubeY_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#CDMhiresZ    = load(root+'/CDM/zoom2/mockCubeZ_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#WDM1p5hiresX = load(root+'/WDM_1p5/zoom2/mockCubeX_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#WDM1p5hiresY = load(root+'/WDM_1p5/zoom2/mockCubeY_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')
#WDM1p5hiresZ = load(root+'/WDM_1p5/zoom2/mockCubeZ_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'.npz')

#HELLO
#rootH    = "/data/titouan/HELLO/"
#CDMHX    = load(rootH+'CDM/g9.20e12/4cMpc/mockCubeX_NoNoise_intrinsicSB_z4.022_Ng512_SBprofile.npz')
#CDMHY    = load(rootH+'CDM/g9.20e12/4cMpc/mockCubeY_NoNoise_intrinsicSB_z4.022_Ng512_SBprofile.npz')
#CDMHZ    = load(rootH+'CDM/g9.20e12/4cMpc/mockCubeY_NoNoise_intrinsicSB_z4.022_Ng512_SBprofile.npz')
#WDM3p5HX = load(rootH+'WDM/g9.20e12/hydro/4cMpc/mockCubeX_NoNoise_intrinsicSB_z4.022_Ng512_SBprofile.npz')
#WDM3p5HY = load(rootH+'WDM/g9.20e12/hydro/4cMpc/mockCubeY_NoNoise_intrinsicSB_z4.022_Ng512_SBprofile.npz')
#WDM3p5HZ = load(rootH+'WDM/g9.20e12/hydro/4cMpc/mockCubeZ_NoNoise_intrinsicSB_z4.022_Ng512_SBprofile.npz')

# ---

#COMPUTE MEAN AND ERRORBAR FOR MEDIAN AND MEAN 
#GIZMO
R       = CDM1X['binmean'] #same for all
#-
CDM1     = nanmedian([CDM1X['sbprofile'],CDM1Y['sbprofile'],CDM1Z['sbprofile']],axis=0)
CDM184   = nanmedian([CDM1X['p84'],CDM1Y['p84'],CDM1Z['p84']],axis=0)
CDM116   = nanmedian([CDM1X['p16'],CDM1Y['p16'],CDM1Z['p16']],axis=0)
#CDM1mean = mean([CDM1X['sbmean'],CDM1Y['sbmean'],CDM1Z['sbmean']],axis=0)
#CDM1err  = sqrt(CDM1X['sberr']**2.+CDM1Y['sberr']**2.+CDM1Z['sberr']**2.)/3.
#-
CDM2     = nanmedian([CDM2X['sbprofile'],CDM2Y['sbprofile'],CDM2Z['sbprofile']],axis=0)
CDM284   = nanmedian([CDM2X['p84'],CDM2Y['p84'],CDM2Z['p84']],axis=0)
CDM216   = nanmedian([CDM2X['p16'],CDM2Y['p16'],CDM2Z['p16']],axis=0)
#CDM2mean = mean([CDM2X['sbmean'],CDM2Y['sbmean'],CDM2Z['sbmean']],axis=0)
#CDM2err  = sqrt(CDM2X['sberr']**2.+CDM2Y['sberr']**2.+CDM2Z['sberr']**2.)/3.
#-
#CDMfblo     = nanmedian([CDMfbloX['sbprofile'],CDMfbloY['sbprofile'],CDMfbloZ['sbprofile']],axis=0)
#CDMfblo84   = nanmedian([CDMfbloX['p84'],CDMfbloY['p84'],CDMfbloZ['p84']],axis=0)
#CDMfblo16   = nanmedian([CDMfbloX['p16'],CDMfbloY['p16'],CDMfbloZ['p16']],axis=0)
#CDMfblomean = mean([CDMfbloX['sbmean'],CDMfbloY['sbmean'],CDMfbloZ['sbmean']],axis=0)
#CDMfbloerr  = sqrt(CDMfbloX['sberr']**2.+CDMfbloY['sberr']**2.+CDMfbloZ['sberr']**2.)/3.
#-
WDM3p51     = nanmedian([WDM3p51X['sbprofile'],WDM3p51Y['sbprofile'],WDM3p51Z['sbprofile']],axis=0)
WDM3p5184   = nanmedian([WDM3p51X['p84'],WDM3p51Y['p84'],WDM3p51Z['p84']],axis=0)
WDM3p5116   = nanmedian([WDM3p51X['p16'],WDM3p51Y['p16'],WDM3p51Z['p16']],axis=0)
#WDM3p51mean = mean([WDM3p51X['sbmean'],WDM3p51Y['sbmean'],WDM3p51Z['sbmean']],axis=0)
#WDM3p51err  = sqrt(WDM3p51X['sberr']**2.+WDM3p51Y['sberr']**2.+WDM3p51Z['sberr']**2.)/3.
#-
WDM3p52     = nanmedian([WDM3p52X['sbprofile'],WDM3p52Y['sbprofile'],WDM3p52Z['sbprofile']],axis=0)
WDM3p5284   = nanmedian([WDM3p52X['p84'],WDM3p52Y['p84'],WDM3p52Z['p84']],axis=0)
WDM3p5216   = nanmedian([WDM3p52X['p16'],WDM3p52Y['p16'],WDM3p52Z['p16']],axis=0)
#WDM3p52mean = mean([WDM3p52X['sbmean'],WDM3p52Y['sbmean'],WDM3p52Z['sbmean']],axis=0)
#WDM3p52err  = sqrt(WDM3p52X['sberr']**2.+WDM3p52Y['sberr']**2.+WDM3p52Z['sberr']**2.)/3.
#-
WDM1p51     = nanmedian([WDM1p51X['sbprofile'],WDM1p51Y['sbprofile'],WDM1p51Z['sbprofile']],axis=0)
WDM1p5184   = nanmedian([WDM1p51X['p84'],WDM1p51Y['p84'],WDM1p51Z['p84']],axis=0)
WDM1p5116   = nanmedian([WDM1p51X['p16'],WDM1p51Y['p16'],WDM1p51Z['p16']],axis=0)
#WDM1p51mean = mean([WDM1p51X['sbmean'],WDM1p51Y['sbmean'],WDM1p51Z['sbmean']],axis=0)
#WDM1p51err  = sqrt(WDM1p51X['sberr']**2.+WDM1p51Y['sberr']**2.+WDM1p51Z['sberr']**2.)/3.
#-
WDM1p52     = nanmedian([WDM1p52X['sbprofile'],WDM1p52Y['sbprofile'],WDM1p52Z['sbprofile']],axis=0)
WDM1p5284   = nanmedian([WDM1p52X['p84'],WDM1p52Y['p84'],WDM1p52Z['p84']],axis=0)
WDM1p5216   = nanmedian([WDM1p52X['p16'],WDM1p52Y['p16'],WDM1p52Z['p16']],axis=0)
#WDM1p52mean = mean([WDM1p52X['sbmean'],WDM1p52Y['sbmean'],WDM1p52Z['sbmean']],axis=0)
#WDM1p52err  = sqrt(WDM1p52X['sberr']**2.+WDM1p52Y['sberr']**2.+WDM1p52Z['sberr']**2.)/3.
#-
#CDMhires     = nanmedian([CDMhiresX['sbprofile'],CDMhiresY['sbprofile'],CDMhiresZ['sbprofile']],axis=0)
#CDMhires84   = nanmedian([CDMhiresX['p84'],CDMhiresY['p84'],CDMhiresZ['p84']],axis=0)
#CDMhires16   = nanmedian([CDMhiresX['p16'],CDMhiresY['p16'],CDMhiresZ['p16']],axis=0)
#CDMhiresmean = mean([CDMhiresX['sbmean'],CDMhiresY['sbmean'],CDMhiresZ['sbmean']],axis=0)
#CDMhireserr  = sqrt(CDMhiresX['sberr']**2.+CDMhiresY['sberr']**2.+CDMhiresZ['sberr']**2.)/3.
#-
#WDM1p5hires     = nanmedian([WDM1p5hiresX['sbprofile'],WDM1p5hiresY['sbprofile'],WDM1p5hiresZ['sbprofile']],axis=0)
#WDM1p5hires16   = nanmedian([WDM1p5hiresX['p84'],WDM1p5hiresY['p84'],WDM1p5hiresZ['p84']],axis=0)
#WDM1p5hires84   = nanmedian([WDM1p5hiresX['p16'],WDM1p5hiresY['p16'],WDM1p5hiresZ['p16']],axis=0)
#WDM1p5hiresmean = mean([WDM1p5hiresX['sbmean'],WDM1p5hiresY['sbmean'],WDM1p5hiresZ['sbmean']],axis=0)
#WDM1p5hireserr  = sqrt(WDM1p5hiresX['sberr']**2.+WDM1p5hiresY['sberr']**2.+WDM1p5hiresZ['sberr']**2.)/3.

#HELLO
#CDMH    = mean([CDMHX['sbprofile'],CDMHY['sbprofile'],CDMHZ['sbprofile']],axis=0)
#WDM3p5H = mean([WDM3p5HX['sbprofile'],WDM3p5HY['sbprofile'],WDM3p5HZ['sbprofile']],axis=0)

# ---

#PLOT
#zoom1p5, fb_boost, CDM vs WDM gn1 (Halo12)
fout1=root+"/SBprofile_Halo1312_gn1_zoom1p5_hydro_fire2_lvl7_CDM_WDM_z3_t"+str(thresh)+".pdf"
fig1, (ax1, ax2) = plt.subplots(2,1,sharex=True,figsize=(6,7),row_heights=[3,1])#height_ratios=[3, 1])
for ax in fig1.get_axes():
	ax.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title(r'$m_b \, = \, 7.18 \cdot 10^5 \; M_\odot$',fontsize=18)
ax1.plot(R, CDM1, 'b-',label=r'CDM')
ax1.fill_between(R, CDM116, CDM184, color='blue',alpha=0.2)
ax1.text(20,1e-20,r'$z=3$',fontsize=16)
#plt.plot(R, CDMmean, 'b--', label=r'CDM, mean')
#plt.fill_between(R, CDMmean-CDMerr, CDMmean+CDMerr,color='blue',alpha=0.3)
ax1.plot(R, WDM3p51, 'g-',label=r'WDM 3.5keV')
ax1.fill_between(R,WDM3p5116,WDM3p5184,color='green',alpha=0.2)
#plt.plot(R,WDM3p5mean,'g--',label=r'WDM 3.5keV, mean')
#plt.fill_between(R,WDM3p5mean-WDM3p5err,WDM3p5mean+WDM3p5err,color='green',alpha=0.3)
ax1.plot(R, WDM1p51, 'r-',label=r'WDM 1.5keV')
ax1.fill_between(R,WDM1p5116,WDM1p5184,color='red',alpha=0.2)
#plt.plot(R,WDM1p5mean,'r--',label=r'WDM 1.5keV, mean')
#plt.fill_between(R,WDM1p5mean-WDM1p5err,WDM1p5mean+WDM1p5err,color='red',alpha=0.3)
#plt.plot(R, CDMH, 'b-',label=r'CDM, HELLO')
#plt.errorbar(R, CDM[1], xerr=0.0, yerr=log10((CDM[1]+CDM[2])/CDM[1]),fmt='r.')
#plt.plot(R, WDM3p5H, 'b:',label=r'WDM 3.5keV, HELLO')
ax2.plot([min(R),max(R)],[1,1],'k-',linewidth=0.2)
ax2.plot(R,WDM3p51/CDM1,'g-')
ax2.plot(R,WDM1p51/CDM1,'r-')
for ax in fig1.get_axes():
	ax.set_xlim([R[0],1000])
ax1.set_ylim([1e-22,3e-16])
ax2.set_ylim([0.7,1.5])
ax2.set_xlabel(r'R [ckpc]',fontsize=18)
ax1.set_ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
ax2.set_ylabel(r'SB/SB$_\textrm{CDM}$', fontsize=18)
ax1.legend(loc=0,prop={'size':14})
plt.savefig(fout1, bbox_inches='tight')

#zoom1p5, fb_boost, CDM vs WDM gn1 (Halo12)
fout2=root+"/SBprofile_Halo1312_gn2_zoom1p5_hydro_fire2_lvl7_CDM_WDM_z3_t"+str(thresh)+".pdf"
fig2, (ax3, ax4) = plt.subplots(2,1,sharex=True,figsize=(6,7),height_ratios=[3, 1])
for ax in fig1.get_axes():
    ax.set_xscale('log')
ax3.set_yscale('log')
ax3.set_title(r'$m_b \, = \, 7.18 \cdot 10^5 \; M_\odot$',fontsize=18)
ax3.plot(R, CDM2, 'b-',label=r'CDM')
ax3.fill_between(R, CDM216, CDM284, color='blue',alpha=0.2)
ax3.text(20,1e-20,r'$z=3$',fontsize=16)
#plt.plot(R, CDMmean, 'b--', label=r'CDM, mean')
#plt.fill_between(R, CDMmean-CDMerr, CDMmean+CDMerr,color='blue',alpha=0.3)
ax3.plot(R, WDM3p52, 'g-',label=r'WDM 3.5keV')
ax3.fill_between(R,WDM3p5216,WDM3p5284,color='green',alpha=0.2)
#plt.plot(R,WDM3p5mean,'g--',label=r'WDM 3.5keV, mean')
#plt.fill_between(R,WDM3p5mean-WDM3p5err,WDM3p5mean+WDM3p5err,color='green',alpha=0.3)
ax3.plot(R, WDM1p52, 'r-',label=r'WDM 1.5keV')
ax3.fill_between(R,WDM1p5216,WDM1p5284,color='red',alpha=0.2)
#plt.plot(R,WDM1p5mean,'r--',label=r'WDM 1.5keV, mean')
#plt.fill_between(R,WDM1p5mean-WDM1p5err,WDM1p5mean+WDM1p5err,color='red',alpha=0.3)
#plt.plot(R, CDMH, 'b-',label=r'CDM, HELLO')
#plt.errorbar(R, CDM[1], xerr=0.0, yerr=log10((CDM[1]+CDM[2])/CDM[1]),fmt='r.')
#plt.plot(R, WDM3p5H, 'b:',label=r'WDM 3.5keV, HELLO')
ax4.plot([min(R),max(R)],[1,1],'k-',linewidth=0.2)
ax4.plot(R,WDM3p52/CDM2,'g-')
ax4.plot(R,WDM1p52/CDM2,'r-')
for ax in fig1.get_axes():
    ax.set_xlim([R[0],1000])
ax3.set_ylim([1e-22,3e-16])
ax4.set_ylim([0.7,1.5])
ax4.set_xlabel(r'R [ckpc]',fontsize=18)
ax3.set_ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
ax4.set_ylabel(r'SB/SB$_\textrm{CDM}$', fontsize=18)
ax3.legend(loc=0,prop={'size':14})
plt.savefig(fout1, bbox_inches='tight')

#zoom2, fb_boost, CDM vs WDM
#fout2=root+"/SBprofile_Halo10_zoom2_hydro_fire2_lvl7_CDM_WDM_z4_t"+str(thresh)+".pdf"
#fig2, (ax3, ax4) = plt.subplots(2,1,sharex=True,figsize=(6,7),height_ratios=[3, 1])
#for ax in fig2.get_axes():
#        ax.set_xscale('log')
#ax3.set_yscale('log')
#ax3.set_title(r'$m_b \, = \, 8.97 \cdot 10^4 \; M_\odot$',fontsize=18)
#ax3.plot(R, CDMhires, 'b-',label=r'CDM')
#ax3.fill_between(R, CDMhires16, CDMhires84, color='blue', alpha=0.2)
#plt.plot(R, CDMhiresmean, 'b--', label=r'CDM, mean')
#plt.fill_between(R, CDMhiresmean-CDMhireserr, CDMhiresmean+CDMhireserr,color='blue',alpha=0.3)
#ax3.plot(R, WDM1p5hires, 'r-',label=r'WDM 1.5keV')
#ax3.fill_between(R,WDM1p5hires16,WDM1p5hires84,color='red',alpha=0.2)
#plt.plot(R,WDM1p5hiresmean,'r--',label=r'WDM 1.5keV, mean')
#plt.fill_between(R,WDM1p5hiresmean-WDM1p5hireserr,WDM1p5hiresmean+WDM1p5hireserr,color='red',alpha=0.3)
#ax4.plot([min(R),max(R)],[1,1],'k-',linewidth=0.2)
#ax4.plot(R,WDM1p5hires/CDMhires,'r-')
#for ax in fig2.get_axes():
#	ax.set_xlim([R[0],1000])
#ax3.set_ylim([1e-22,3e-16])
#ax4.set_ylim([0.7,1.5])
#ax4.set_xlabel(r'R [ckpc]',fontsize=18)
#ax3.set_ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
#ax4.set_ylabel(r'SB/SB$_\textrm{CDM}$', fontsize=18)
#ax3.legend(loc=0,prop={'size':14})
#plt.savefig(fout2, bbox_inches='tight')

#zoom1p5, fb vs fb_boost, CDM 
#fout3=root+"/SBprofile_Halo10_zoom1p5_hydro_fire2_lvl7_FB_z4_t"+str(thresh)+".pdf"
#fig3, (ax5, ax6) = plt.subplots(2,1,sharex=True,figsize=(6,7),height_ratios=[3, 1])
#for ax in fig3.get_axes():
#        ax.set_xscale('log')
#ax5.set_yscale('log')
#ax5.set_title(r'CDM, $m_b \, = \, 7.18 \cdot 10^5 \; M_\odot$',fontsize=18)
#ax5.text(300,5e-18,r'$z=4$',fontsize=16)
#ax5.plot(R, CDM, 'b-',label=r'$v_{ej}$ = 4377 km/s')
#ax5.fill_between(R, CDM16, CDM84, color='blue',alpha=0.2)
#plt.plot(R, CDMmean, 'b--', label=r'$v_{ej}$ = 4377 km/s, mean')
#plt.fill_between(R, CDMmean-CDMerr, CDMmean+CDMerr,color='blue',alpha=0.3)
#ax5.plot(R, CDMfblo, 'm-',label=r'$v_{ej}$ = 3095 km/s') 
#ax5.fill_between(R, CDMfblo16, CDMfblo84, color='magenta',alpha=0.2)
#plt.plot(R, CDMfblomean, 'm--', label=r'$v_{ej}$ = 3095 km/s, mean')
#plt.fill_between(R, CDMfblomean-CDMfbloerr, CDMfblomean+CDMfbloerr,color='magenta',alpha=0.3)
#ax6.plot([min(R),max(R)],[1,1],'k-',linewidth=0.2)
#ax6.plot(R,CDMfblo/CDM,'m-')
#ax6.plot(R,WDM1p5/CDM,'r--',label=r'WDM 1.5 keV')
#for ax in fig3.get_axes():
#	ax.set_xlim([R[0],1000])
#ax5.set_ylim([1e-22,3e-16])
#ax6.set_ylim([0.6,1.75])
#ax6.set_xlabel(r'R [ckpc]',fontsize=18)
#ax5.set_ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
#ax6.set_ylabel(r'SB_{lowFB}/SB$_\textrm{HiFB}$', fontsize=18)
#ax5.legend(loc=0,prop={'size':14})
#ax6.legend(loc=0,prop={'size':14})
#plt.savefig(fout3, bbox_inches='tight')

#ratios
#zoom1p5
#err3p5 = abs(WDM3p5mean/CDMmean)*sqrt((WDM3p5err/WDM3p5mean)**2.+(CDMerr/CDMmean)**2.)
#err1p5 = abs(WDM1p5mean/CDMmean)*sqrt((WDM1p5err/WDM1p5mean)**2.+(CDMerr/CDMmean)**2.)
#errFB  = abs(CDMfblomean/CDMmean)*sqrt((CDMfbloerr/CDMfblomean)**2.+(CDMerr/CDMmean)**2.)
#fout4=root+"/SBprofile_ratios_Halo10_zoom1p5_hydro_fire2_lvl7_z4_t"+str(thresh)+".pdf"
#plt.figure(4)
#plt.xscale('log')
#plt.plot([min(R),max(R)],[1,1],'k-')
#plt.title(r'$m_b \, = \, 7.18 \cdot 10^5 \; M_\odot$',fontsize=18)
#plt.plot(R,WDM3p5/CDM,'g-',label=r'WDM 3.5 keV')
#plt.plot(R,WDM1p5/CDM,'r-',label=r'WDM 1.5 keV')
#plt.plot(R,WDM3p5mean/CDMmean,'g-',label=r'WDM 3.5 keV')
#plt.fill_between(R,WDM3p5mean/CDMmean-err3p5,WDM3p5mean/CDMmean+err3p5,color='green',alpha=0.3)
#plt.plot(R,WDM1p5mean/CDMmean,'r-',label=r'WDM 1.5 keV')
#plt.fill_between(R,WDM1p5mean/CDMmean-err1p5,WDM1p5mean/CDMmean+err1p5,color='red',alpha=0.3)
#plt.plot(R,CDMfblomean/CDMmean,'m-',label=r'$v_{ej}$ = 3095 km/s')
#plt.fill_between(R,CDMfblomean/CDMmean-errFB,CDMfblomean/CDMmean+errFB,color='magenta',alpha=0.3)
#plt.ylim([0.7,1.3])
#plt.xlim([R[0],1000])
#plt.xlabel(r'R [ckpc]',fontsize=20)
#plt.ylabel(r'SB/SB$_\textrm{CDM}$', fontsize=20)
#plt.legend(loc=0,prop={'size':14})
#plt.savefig(fout4, bbox_inches='tight')

#zoom2
#err1p5hires = abs(WDM1p5hiresmean/CDMhiresmean)*sqrt((WDM1p5hireserr/WDM1p5hiresmean)**2.+(CDMhireserr/CDMhiresmean)**2.)
#fout5=root+"/SBprofile_ratios_Halo10_zoom2_hydro_fire2_lvl7_z4_t"+str(thresh)+".pdf"
#plt.figure(5)
#plt.xscale('log')
#plt.plot([min(R),max(R)],[1,1],'k-')
#plt.title(r'$m_b \, = \, 8.97 \cdot 10^4 \; M_\odot$',fontsize=18)
#plt.plot(R,WDM1p5hires/CDMhires,'r-',label=r'WDM 1.5 keV')
#plt.plot(R,WDM1p5hiresmean/CDMhiresmean,'r-',label=r'WDM 1.5 keV')
#plt.fill_between(R,WDM1p5hiresmean/CDMhiresmean-err1p5hires,WDM1p5hiresmean/CDMhiresmean+err1p5hires,color='red',alpha=0.3)
#plt.ylim([0.7,1.3])
#plt.xlim([R[0],1000])
#plt.xlabel(r'R [ckpc]',fontsize=20)
#plt.ylabel(r'SB/SB$_\textrm{CDM}$', fontsize=20)
#plt.legend(loc=0,prop={'size':14})
#plt.savefig(fout5, bbox_inches='tight')

#zoom1p5/zoom2p5
#errCDMz = abs(CDMhiresmean/CDMmean)*sqrt((CDMhireserr/CDMhiresmean)**2.+(CDMerr/CDMmean)**2.)
#err1p5z = abs(WDM1p5hiresmean/WDM1p5mean)*sqrt((WDM1p5hireserr/WDM1p5hiresmean)**2.+(WDM1p5err/WDM1p5mean)**2.)
#fout6=root+"/SBprofile_ratios_Halo10_zooms_hydro_fire2_lvl7_z4_t"+str(thresh)+".pdf"
#plt.figure(6)
#plt.xscale('log')
#plt.text(20,3,r'$z=4$',fontsize=16)
#plt.plot([min(R),max(R)],[1,1],'k-')
#plt.plot(R,CDMhires/CDM,'b-',label=r'CDM')
#plt.plot(R,WDM1p5hires/WDM1p5,'r-',label=r'WDM 1.5 keV')
#plt.plot(R,CDMhiresmean/CDMmean,'b-',label=r'CDM')
#plt.fill_between(R,CDMhiresmean/CDMmean-errCDMz,CDMhiresmean/CDMmean+errCDMz,color='blue',alpha=0.3)
#plt.plot(R,WDM1p5hiresmean/WDM1p5mean,'r-',label=r'WDM 1.5 keV')
#plt.fill_between(R,WDM1p5hiresmean/WDM1p5mean-err1p5z,WDM1p5hiresmean/WDM1p5mean+err1p5z,color='red',alpha=0.3)
#plt.ylim([0.0,4.0])
#plt.xlim([R[0],1000])
#plt.xlabel(r'R [ckpc]',fontsize=18)
#plt.ylabel(r'SB/SB$_\textrm{low}$', fontsize=18)
#plt.legend(loc=0,prop={'size':14})
#plt.savefig(fout6, bbox_inches='tight')

#plt.show()
