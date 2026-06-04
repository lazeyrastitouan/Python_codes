from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
import glob
from decimal import Decimal

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)

edges  = ['100','200']
red    = 3.5
suf    = ['','metal005_','','dust-5_','metal005_','','metal005_']
#sims   = ['DaLya $n_\star = 10$','DaLya $Z_\star = 0.05 Z_\odot$', 'COLIBRE $n_\star = 10$', 'COLIBRE dust $M_{frac} \, < 10^{-5}$', 'COLIBRE $Z_\star = 0.05Z_\odot$', 'HELLO $n_\star = 10$', 'HELLO $Z_\star = 0.05 Z_\odot$']
sims   = ['DaLya $n_\star = 10$','DaLya', 'COLIBRE $n_\star = 10$', 'COLIBRE dust $M_{frac} \, < 10^{-5}$', 'COLIBRE', 'HELLO $n_\star = 10$', 'HELLO']
grids  = [ [['1','3'],['4','11']], [['1','3'],['4','11']], [['1','3'],['4','15']], [['1','3'],['4','15']], [['1','3'],['4','15']], [['1','1'],['2','2']], [['1','1'],['2','2']] ]
Ngrids = [ [15,30] , [15,30] , [9,36] , [9,36] , [9,36] , [3,3], [3,3] ]
sg     = shape(grids) #num sims, num mbin, low grid - high grid
colors = ['cyan','b','magenta','orange','r','lime','green']
msuf   = ['mhigh','mlow']
paths  = ['/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/','/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/HELLO/CDM/g9.20e12/2cMpc/','/data/titouan/HELLO/CDM/g9.20e12/2cMpc/metal005/']

obs    = loadtxt('/data/titouan/MQN/SB_PDF_100_200ckpc_densityFalse.txt').T
edges2 = obs[0]-0.05
edges2 = append(edges2,obs[0,len(obs[0])-1]+0.05)
Nneb   = 14 #*27678/19575.

#mass loop
for m in range(sg[1]):
    fout = '/data/titouan/COLIBRE/Hist_SB_Npix_R100-200_dalya_zoom1p5_colibre_m6_hello_'+msuf[m]+'_z'+str(red)+'_n5e-20_nonoise_dust_metal.pdf'
    #fig1 mlow, fig2  mhigh
    fig, axs = plt.subplots(1,1)#,sharey=True,sharex=True,figsize=(10,5))
    #plt.subplots_adjust(wspace=0.09)
    axs.set_ylabel(r'$\left\langle N_{pix}\right\rangle$', fontsize=15)
    #if msuf[m] == 'mlow':
    #    axs.set_title(r'$12.0 \leq \log M \leq 12.5$',fontsize=14)
    #elif msuf[m] == 'mhigh':
    #    axs.set_title(r'$12.5 \leq \log M \leq 13.1$',fontsize=14)
    axs.set_yscale('log')
    axs.set_ylim(1,6e2) #1e3
    axs.set_xlim(-18.,-12.3)
    #axs.text(-17.85,125,r''+edges[0]+' $\leq r \leq$ '+edges[1],fontsize=11)
    #axs.text(-17.85,80,r'z = '+str(red),fontsize=11)
    axs.text(-13.5,150,r'$Z_\star = 0.05 Z_\odot$',fontsize=13)
    axs.set_xlabel(r'$\log$(SB$_{\rm intr}$) [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=14)
    #plot obs
    axs.stairs(obs[1]/Nneb,edges2,fill=False,color='black',ls='-.',label=r'MAGG')
    #plot sims
    for s in range(len(suf)):
        data = load(paths[s]+'Hist_SB_Npixel_'+suf[s]+'n5e-20_s1_Ng256_observed_cubex_neb_singleobj_NB_nonoise_z'+str(red)+'_grids'+grids[s][m][0]+'-'+grids[s][m][1]+'_R'+edges[0]+'-'+edges[1]+'.npz')
        H    = data['hist']/Ngrids[s][m]
        E    = log10(10.**data['edges']*(1.+red)**4.)
        if suf[s] != "":
            axs.stairs(H,E,fill=False,ls='--',color=colors[s],label=r''+sims[s])
    axs.legend(loc='upper left',ncol=2,prop={'size':12},frameon=False)
    
    plt.savefig(fout, bbox_inches='tight')
    print("Saved figure to "+fout)
    plt.show()
