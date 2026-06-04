from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
import glob
from decimal import Decimal

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

edges  = ['100','200']
red    = 3.5
sims   = ['DaLya', 'COLIBRE m6', 'HELLO']
smooth = ['1','1.5','5']
slbl   = ['1','1p5','5']
grids  = [ [['1','3'],['4','11']], [['1','3'],['4','15']], [['1','1'],['2','2']] ]
Ngrids = [ [15,30] , [9,36] , [3,3] ]
sg     = shape(grids) #num sims, num mbin, low grid - high grid
colors = [['b','orange','g'],['b','orange','g'],['b','orange','g']]#,['cyan','orange','lime']]
lsty   = ['-','-','-']
msuf   = ['mhigh','mlow']
paths  = ['/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/HELLO/CDM/g9.20e12/2cMpc/']

binwidth = 0.1 #in log10
obs      = loadtxt('/data/titouan/MQN/SB_PDF_100_200ckpc_densityFalse.txt').T
edges2   = obs[0]-binwidth/2.
edges2   = append(edges2,obs[0,len(obs[0])-1]+binwidth/2.)
Nneb     = 14 #*27678/19575.

means_100_200 = [1.2674083879554989e-15,2.0109307216616325e-15,1.7827606076985994e-15]

#mass loop
for m in range(sg[1]):
    for sm in range(len(smooth)):
        if sm == 0 or sm == 2:
            continue
        fig, axs = plt.subplots(1,1)
        fout = '/data/titouan/COLIBRE/Hist_SB_Npix_R100-200_dalya_zoom1p5_colibre_m6_'+msuf[m]+'_z'+str(red)+'_n5e-20_nonoise_smooth'+slbl[sm]+'.pdf'
        axs.set_ylabel(r'$\textrm{d} N_{pix} / \textrm{dlog} SB$', fontsize=18)
        #if msuf[m] == 'mlow':
        #    axs.set_title(r'$12.0 \leq \log M \leq 12.5$',fontsize=14)
        #elif msuf[m] == 'mhigh':
        #    axs.set_title(r'$12.5 \leq \log M \leq 13.1$',fontsize=14)
        axs.set_yscale('log')
        axs.set_xscale('log')
        axs.set_ylim(4,2e3) #1e3
        axs.set_xlim(3e-18,1e-13)
        #axs.text(-17.85,125,r''+edges[0]+' $\leq r \leq$ '+edges[1],fontsize=11)
        axs.text(4e-18,1000,r'$n_\star=10$ cm$^{-3}$',fontsize=16)
        axs.set_xlabel(r'$\log$(SB$_{\rm intr}$) [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
        #plot obs
        if sm ==1:
            axs.axvline(1e-16, 0, 1, color='grey', linestyle=':')
            axs.axvline(3.4465102094748147e-15,0, 1, color='black', linestyle='--')
            #axs.text(-17.85,300,r'$n_\star = 10$',fontsize=13)
            #axs.text(-14.65,300,r'FWHM = 27.5 ckpc = 0.8"',fontsize=12)
        axs.stairs(obs[1]/Nneb/binwidth,10**edges2,fill=False,color='black',ls='-',label=r'Observations')
        #plot sims
        for s in range(len(sims)):
            if s > 1:
                continue
            data = load(paths[s]+'Hist_SB_Npixel_n5e-20_s'+smooth[sm]+'_Ng256_observed_cubex_neb_singleobj_NB_nonoise_z'+str(red)+'_grids'+grids[s][m][0]+'-'+grids[s][m][1]+'_R'+edges[0]+'-'+edges[1]+'.npz')
            H    = data['hist']/Ngrids[s][m]/binwidth
            E    = 10.**data['edges']*(1.+red)**4.
            axs.stairs(H,E,fill=False,ls=lsty[sm],color=colors[sm][s],label=r''+sims[s])#+' '+slbl[sm])
            if m==1 and sm ==1:
                axs.axvline(means_100_200[s], 0, 1, color=colors[sm][s], linestyle='--')
                #axs.annotate("",xy=(means_100_200[s],4.5),xytext=(means_100_200[s], 9),arrowprops=dict(arrowstyle="->",linewidth=1,color=colors[sm][s]))
        if sm == 2:
            axs.text(-17.85,300,r'FWHM = 91.8 ckpc = 2.7"',fontsize=12)
            axs.text(-17.85,200,r'$n_\star = 10$',fontsize=13)
        elif sm == 0:
            axs.axvline(-16, 0, 1, color='grey', linestyle=':')
            #axs.text(-14.65,300,r'FWHM = 17.2  ckpc = 0.5"',fontsize=12)

        #axs.legend(loc='upper right',prop={'size':12},frameon=False)
    
        plt.savefig(fout, bbox_inches='tight')
        print("Saved figure to "+fout)
        plt.show()
