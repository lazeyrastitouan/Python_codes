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
sims   = ['HELLO','DaLya fiduc.', 'DaLya high res.', 'COLIBRE m6', 'COLIBRE m5']
smooth = ['1.5']
slbl   = ['1p5']
grids  = [ [['1','1']],[['4','11']],[['3','3']],[['4','15']],[['4','14']] ]
Ngrids = [ [12],[30],[6],[36],[33] ]
sg     = shape(grids) #num sims, num mbin, low grid - high grid
colors = [['g','b','b','r','r']]
lsty   = ['--','-','--','-','--']
lw     = [1.5,1,1.5,1,1.5]
msuf   = ['mlow']
paths  = ['/data/titouan/HELLO/CDM/','/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/','/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/COLIBRE/L0050N1504/']

binwidth = 0.1 #in log10
obs      = loadtxt('/data/titouan/MQN/SB_PDF_100_200ckpc_densityFalse.txt').T
edges2   = obs[0]-binwidth/2.
edges2   = append(edges2,obs[0,len(obs[0])-1]+binwidth/2.)
Nneb     = 14 #*27678/19575.

means_100_200 = [1.3946375907597338e-15,1.2674083879554989e-15,1.1854660651228243e-15,2.0109307216616325e-15,2.9482917005876948e-15]

#mass loop
for m in range(sg[1]):
    for sm in range(len(smooth)):
        fig, axs = plt.subplots(1,1)
        fout = '/data/titouan/COLIBRE/Hist_SB_Npix_R100-200_dalya_zoom1p5_2_colibre_m6_m5_hello_'+msuf[m]+'_z'+str(red)+'_n5e-20_nonoise_smooth'+slbl[sm]+'.pdf'
        axs.set_ylabel(r'$\textrm{d} N_{pix} / \textrm{dlog} SB$', fontsize=18)
        #if msuf[m] == 'mlow':
        #    axs.set_title(r'$12.0 \leq \log M \leq 12.5$',fontsize=14)
        #elif msuf[m] == 'mhigh':
        #    axs.set_title(r'$12.5 \leq \log M \leq 13.1$',fontsize=14)
        axs.set_yscale('log')
        axs.set_xscale('log')
        axs.set_ylim(4,3e3) #1e3
        axs.set_xlim(3e-18,2e-13)
        #axs.text(-17.85,125,r''+edges[0]+' $\leq r \leq$ '+edges[1],fontsize=11)
        axs.text(4e-18,1700,r'$n_\star=10$ cm$^{-3}$',fontsize=16)
        axs.set_xlabel(r'$\log$(SB$_{\rm intr}$) [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
        #plot obs
        #axs.axvline(1e-16, 0, 1, color='grey', linestyle=':')
        #axs.axvline(3.4465102094748147e-15,0, 1, color='black', linestyle='--')
        #axs.text(-17.85,300,r'$n_\star = 10$',fontsize=13)
        #axs.text(-14.65,300,r'FWHM = 27.5 ckpc = 0.8"',fontsize=12)
        axs.stairs(obs[1]/Nneb/binwidth,10**edges2,fill=False,color='black',ls='-',label=r'Observations',linewidth=2)
        axs.annotate("",xy=(3.4465102094748147e-15,4.5),xytext=(3.4465102094748147e-15, 9),arrowprops=dict(arrowstyle="->",linewidth=1,color='black'))
        #plot sims
        for s in range(len(sims)):
            data = load(paths[s]+'Hist_SB_Npixel_n5e-20_s'+smooth[sm]+'_Ng256_observed_cubex_neb_singleobj_NB_nonoise_z'+str(red)+'_grids'+grids[s][m][0]+'-'+grids[s][m][1]+'_R'+edges[0]+'-'+edges[1]+'.npz')
            H    = data['hist']/Ngrids[s][m]/binwidth
            E    = 10.**data['edges']*(1.+red)**4.
            axs.stairs(H,E,fill=False,ls=lsty[s],color=colors[sm][s],linewidth=lw[s],label=r''+sims[s])#+' '+slbl[sm])
            #axs.axvline(means_100_200[s], 0, 1, color=colors[sm][s], linestyle='--')
            axs.annotate("",xy=(means_100_200[s],4.5),xytext=(means_100_200[s], 9),arrowprops=dict(linestyle=lsty[s],arrowstyle="->",linewidth=1,color=colors[sm][s]))
            #axs.text(-14.65,300,r'FWHM = 17.2  ckpc = 0.5"',fontsize=12)

        axs.legend(loc="lower center",bbox_to_anchor=(0.5, 1.0),ncol=3,prop={'size':12},frameon=False)
    
        plt.savefig(fout, bbox_inches='tight')
        print("Saved figure to "+fout)
        plt.show()
