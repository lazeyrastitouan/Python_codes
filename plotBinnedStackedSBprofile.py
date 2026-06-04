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

# input parameters
parser  = argparse.ArgumentParser(description='Stack 3 LOS SB profiles and plot + ratios. Args: < SB threshold (float) >')
parser.add_argument('threshold',nargs=1,type=str, help='SB threshold of profiles to stack')
args   = parser.parse_args()
thresh = args.threshold[0]

root   = '/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/Halo10_lvl7/sftresh_cool_low_fb_boost/z4/'

# ---

Mtab = arange(9.0,12.5,0.5).tolist()
dM   = 0.5

fout = root+"/SBprofile_Halo10_zoom1p5_hydro_fire2_lvl7_CDM_WDM_z4_t"+str(thresh)+"_Mbins.pdf"
nbin = len(Mtab)
rows = int(nbin/3)+(nbin%3)
fig, axes = plt.subplots(rows,3,sharex=True,sharey=True,figsize=(16,14),height_ratios=[1,1,1])
fig.subplots_adjust(hspace=0.3)
for ax in fig.get_axes():
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim([1e-22,3e-16])
axes[2,0].set_xlabel(r'R [ckpc]',fontsize=24,labelpad=16)
axes[1,0].set_ylabel(r'SB/$(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=24, labelpad=16)

for m in Mtab:
	im = Mtab.index(m)
	i = int(im/3)
	j = im%3
	milieu = round((m+m+dM)/2.,2)

	#LOAD
	CDMX         = load(root+'/CDM/zoom1p5/mockCubeX_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	CDMY         = load(root+'/CDM/zoom1p5/mockCubeY_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	CDMZ         = load(root+'/CDM/zoom1p5/mockCubeZ_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	WDM1p5X      = load(root+'/WDM_1p5/zoom1p5/mockCubeX_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	WDM1p5Y      = load(root+'/WDM_1p5/zoom1p5/mockCubeY_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	WDM1p5Z      = load(root+'/WDM_1p5/zoom1p5/mockCubeZ_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	WDM3p5X      = load(root+'/WDM_3p5/zoom1p5/mockCubeX_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	WDM3p5Y      = load(root+'/WDM_3p5/zoom1p5/mockCubeY_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')
	WDM3p5Z      = load(root+'/WDM_3p5/zoom1p5/mockCubeZ_SB_s1_n0_z4.0_Ng512_observed_SBprofile_t'+thresh+'_Mmin_'+str(m)+'_Mmax_'+str(m+0.5)+'.npz')

	# ---

	#COMPUTE MEDIAN OF MEDIANS AND PRECENTILES 
	R       = CDMX['binmean'] #same for all
	#-
	CDM     = nanmedian([CDMX['sbprofile'],CDMY['sbprofile'],CDMZ['sbprofile']],axis=0)
	CDM84   = nanmedian([CDMX['p84'],CDMY['p84'],CDMZ['p84']],axis=0)
	CDM16   = nanmedian([CDMX['p16'],CDMY['p16'],CDMZ['p16']],axis=0)
	#-
	WDM3p5     = nanmedian([WDM3p5X['sbprofile'],WDM3p5Y['sbprofile'],WDM3p5Z['sbprofile']],axis=0)
	WDM3p584   = nanmedian([WDM3p5X['p84'],WDM3p5Y['p84'],WDM3p5Z['p84']],axis=0)
	WDM3p516   = nanmedian([WDM3p5X['p16'],WDM3p5Y['p16'],WDM3p5Z['p16']],axis=0)
	#-
	WDM1p5     = nanmedian([WDM1p5X['sbprofile'],WDM1p5Y['sbprofile'],WDM1p5Z['sbprofile']],axis=0)
	WDM1p584   = nanmedian([WDM1p5X['p84'],WDM1p5Y['p84'],WDM1p5Z['p84']],axis=0)
	WDM1p516   = nanmedian([WDM1p5X['p16'],WDM1p5Y['p16'],WDM1p5Z['p16']],axis=0)

	# ---

	#PLOT
	#zoom1p5, fb_boost, CDM vs WDM
	axes[i,j].set_title(r'$\log(M/M_\odot) \, =$ '+str(milieu),fontsize=16)
	axes[i,j].set_xlim([R[0],1000])
	if i == 0 and j == 2:
		axes[i,j].plot(R, CDM, 'b-',label=r'CDM')
		axes[i,j].plot(R, WDM3p5, 'g-',label=r'WDM 3.5keV')
		axes[i,j].plot(R, WDM1p5, 'r-',label=r'WDM 1.5keV')
	else:
		axes[i,j].plot(R, CDM, 'b-')
		axes[i,j].plot(R, WDM3p5, 'g-')
		axes[i,j].plot(R, WDM1p5, 'r-')
	axes[i,j].fill_between(R, CDM16, CDM84, color='blue',alpha=0.2)
	axes[i,j].fill_between(R,WDM3p516,WDM3p584,color='green',alpha=0.2)
	axes[i,j].fill_between(R,WDM1p516,WDM1p584,color='red',alpha=0.2)
	
	axes[i,j].tick_params(axis='both', which='both', labelsize=14,labelleft=True, labelbottom=True)

	del R, CDMX, CDMY, CDMZ, WDM1p5X, WDM1p5Y, WDM1p5Z, WDM3p5X, WDM3p5Y, WDM3p5Z, CDM, CDM84, CDM16, WDM3p5, WDM3p584, WDM3p516, WDM1p5, WDM1p584, WDM1p516 

axes[0,2].legend(bbox_to_anchor=(0.75, -2.0),ncol=1,frameon=False,prop={'size':16})

#customize subplots layout
pos = axes[2,1].get_position()
for ax in axes.flat[nbin:]:
	ax.remove()
axes[2,0].set_position([pos.x0,pos.y0,pos.width,pos.height])

plt.savefig(fout, bbox_inches='tight')
plt.show()
