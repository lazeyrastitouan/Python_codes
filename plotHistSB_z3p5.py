from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
import glob
from decimal import Decimal

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

edges  = ['100','200']
red    = 3.5
suf    = ['','','dust-7_','dust-5_','metal005_']
sims   = ['DaLya $n_\star = 10$','COLIBRE $n_\star = 10$', 'COLIBRE dust frac $< 10^{-7}$', 'COLIBRE dust frac $< 10^{-5}$', 'COLIBRE $Z<0.05Z_\odot$']
grids  = [ [['1','3'],['4','11']] , [['1','3'],['4','15']] , [['1','3'],['4','15']] , [['1','3'],['4','15']] , [['1','3'],['4','15']] ]
sg     = shape(grids) #num sims, num mbin, low grid - high grid
colors = ['b','r','orange','goldenrod','magenta']
msuf   = ['mhigh','mlow']
paths  = ['/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/COLIBRE/L0050N0752/' ]

obs    = loadtxt('/data/titouan/MQN/SB_PDF_100_200ckpc.txt').T
edges2 = obs[0]-0.07705
edges2 = append(edges2,obs[0,len(obs[0])-1]+0.07705)

#mass loop
for m in range(sg[1]):
    fout = '/data/titouan/COLIBRE/Hist_SB_annuliR_dalya_zoom1p5_colibre_m6_'+msuf[m]+'_z'+str(red)+'_n5e-20_nonoise_dust_metal.pdf'
    #fig1 mlow, fig2  mhigh
    fig, axs = plt.subplots(1,1)#,sharey=True,sharex=True,figsize=(10,5))
    #plt.subplots_adjust(wspace=0.09)
    axs.set_ylabel(r'PDF', fontsize=14)
    if msuf[m] == 'mlow':
        axs.set_title(r'$12.0 \leq \log M \leq 12.5$',fontsize=14)
    elif msuf[m] == 'mhigh':
        axs.set_title(r'$12.5 \leq \log M \leq 13.1$',fontsize=14)
    axs.set_yscale('log')
    axs.set_ylim(0.002,6.0)
    axs.set_xlim(-17.5,-12.5)
    axs.text(-17.35,6e-1,r''+edges[0]+' $\leq r \leq$ '+edges[1],fontsize=12)
    axs.text(-17.35,9e-1,r'z = '+str(red),fontsize=12)
    axs.set_xlabel(r'SB$_{\rm intr}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=12)
    #plot obs
    axs.stairs(obs[1],edges2,fill=False,color='black',ls='-.',label=r'MAGG')
    #plot sims
    for s in range(len(suf)):
        data = load(paths[s]+'Hist_SB_'+suf[s]+'n5e-20_s1_Ng256_observed_cubex_neb_singleobj_NB_nonoise_z'+str(red)+'_grids'+grids[s][m][0]+'-'+grids[s][m][1]+'_R'+edges[0]+'-'+edges[1]+'.npz')
        H    = data['hist']
        E    = log10(10.**data['edges']*(1.+red)**4.)
        axs.stairs(H,E,fill=False,ls='--',color=colors[s],label=r''+sims[s])
    axs.legend(loc='upper left',ncol=2,prop={'size':10},frameon=False)
    
    plt.savefig(fout, bbox_inches='tight')
    print("Saved figure to "+fout)
    plt.show()
