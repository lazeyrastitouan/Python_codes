from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
import glob
from decimal import Decimal
import re

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

savemean = True

zm      = '1p5'
red     = 3.5
edges   = [100,200]
sims    = ['DaLya fiducial', 'COLIBRE m6', 'HELLO']
grids   = [ [[4,11]], [[4,15]], [[1,2]] ] #LR[ [[2,3]], [[4,14]], [[1,2]] ] #HR[[[4,11]], [[4,15]], [[2,2]] ]
sg      = shape(grids) #num sims, num Mmbin, low grid - high grid
colors  = ['b','orange','g']
paths   = ['/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom'+zm+'_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/2cMpc/','/data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/2cMpc/','/data/titouan/HELLO/CDM/g*e12/2cMpc/']
Mf      = ['mlow']
smooth  = 1.5

#observations
obs = loadtxt("/data/titouan/Andrea_MAGG_data/SB_values_100_200ckpc_CoveringFraction.txt").T
sb  = obs[0]

for m in Mf:
    im   = Mf.index(m)
    fout = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/CoveringFraction_z'+str(red)+'_dalya1p5_smooth1p5.pdf'
    #fig1 mlow, fig2  mhigh, axs are different annuli if needed
    if len(edges)-1 == 1:
        fig, axs = plt.subplots(1,len(edges)-1)
        axs = [axs]
    else:
        fig, axs = plt.subplots(1,len(edges)-1,sharey=True,sharex=True,figsize=(18,6))
        plt.subplots_adjust(wspace=0.07)
    for i in range(len(edges)-1):
        axs[i].set_xscale('log')
        axs[i].set_yscale('log')
        axs[i].set_ylim(4e-4,1.5)
        #axs[i].text(5e-17,0.02,r''+str(edges[i])+' $\leq r \leq$ '+str(edges[i+1]),fontsize=12)
        axs[i].text(1e-14,0.7,r'$n_\star=10$ cm$^{-3}$',fontsize=18)
        #axs[i].text(5e-17,4e-3,r'FWHM = 27.5 ckpc = 0.8"',fontsize=12)
    axs[0].set_xlabel(r'SB$_{\rm intr}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
    axs[0].set_ylabel(r'covering fraction', fontsize=18)
    #if m == 'mlow':
    #    axs[0].set_title(r'$12.0 \leq \log M \leq 12.5$',fontsize=14)
    #elif m == 'mhigh':
    #    axs[0].set_title(r'$12.5 \leq \log M \leq 13.1$',fontsize=14)
    
    for s in range(sg[0]):
        if s > 0:
            continue
        mfcov, sfcov = empty((0,len(edges)-1)), empty((0,len(edges)-1)) #final shape will be number SB tresh, len(edges)-1
        for ie in range(len(edges)-1):
            fcov  = empty((0,len(sb))) #final shape is num of files, len(sb)
            files = []
            for f in glob.glob(paths[s]+'mockCube*_noise5e-20_*gn*_SB_s'+str(smooth)+'_z'+str(red)+'*_Ng256_observed_cubex_neb_singleobj_NB_nonoise_CoveringFraction_'+str(edges[ie])+'-'+str(edges[ie+1])+'.txt'):
                match = re.search(r"gn(\d+)",f)
                if match:
                    fi = int(match.group(1))
                    if grids[s][im][0] <= fi <= grids[s][im][1]:
                        files.append(f)
            #add H1312 in DaLya
            if sims[s] == 'DaLya': #DaLya HR
                for f in glob.glob('/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom'+zm+'_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/2cMpc/mockCube*_noise5e-20_*gn*_SB_s'+str(smooth)+'_z'+str(red)+'*_Ng256_observed_cubex_neb_singleobj_NB_nonoise_CoveringFraction_'+str(edges[ie])+'-'+str(edges[ie+1])+'.txt'):
                    match = re.search(r"gn(\d+)",f)
                    if match:
                        fi = int(match.group(1))
                        if grids[s][im][0]-1 <= fi <= grids[s][im][1]:
                            files.append(f)
            print(files)
            for f in files:
                fcov = vstack((fcov,genfromtxt(f,usecols=(1,))))
            mfcov = mean(fcov,axis=0)
            sfcov = std(fcov,axis=0)/sqrt(len(files)-1)
            del files

            axs[ie].plot(sb,mfcov,ls='-',marker='.',color=colors[s],label=r''+sims[s])
            axs[ie].errorbar(sb,mfcov,xerr=0.0,yerr=sfcov,ls='-',color=colors[s])

    for i in range(len(edges)-1):
        obs = loadtxt("/data/titouan/Andrea_MAGG_data/SB_values_"+str(edges[i])+"_"+str(edges[i+1])+"ckpc_CoveringFraction.txt").T
        axs[i].plot(obs[0,0:-1],obs[1,0:-1],'k-',label=r'Observations, Travascio+ in prep.')
        axs[i].set_xlim(min(obs[0]),1e-13)
        axs[i].legend(loc='lower left',frameon=False,prop={'size':16})
    plt.savefig(fout, bbox_inches='tight')
    print("Saved figure to "+fout)
    plt.show()
