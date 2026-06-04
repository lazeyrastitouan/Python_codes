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

edges   = ['100','200']
red     = [3.0,3.5,4.0]

grids  = [ [ [['1','4'],['5','13']] , [['1','3'],['4','11']], [['1','1'],['2','7']] ] , [ [['1','7'],['8','12']], [['1','3'],['4','15']], [['0','0'],['1','6']] ] , [ [['0','0'],['0','0']], [['1','1'],['2','2']], [['1','1'],['0','0']] ] ]
sg     = shape(grids) #num sim type, num reds, num mbins, low grid - high grid
colors = ['b','r','g']
sims   = ['DaLya','COLIBRE','HELLO']
msuf   = ['mhigh','mlow']
paths  = ['/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/','/data/titouan/COLIBRE/L0050N0752/','/data/titouan/HELLO/CDM/g9.20e12/2cMpc/']
obs    = ['/data/titouan/MQN/Hist_SB_OE_Lya_z3.2_R100-200.npz','/data/titouan/MQN/Hist_SB_MQN03_Lya_MergedSegCube.OE_z3.2_R100-200.npz', '/data/titouan/MQN/SB_PDF_100_200ckpc.txt']

#mass loop
for m in range(sg[2]):
    #fig1 mlow, fig2  mhigh
    fout = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Hist_SB_annuliR_dalya_zoom1p5_colibre_m6_hello_'+msuf[m]+'_compObs.pdf'
    fig, axs = plt.subplots(1,3,sharey=True,sharex=True,figsize=(18,6))
    plt.subplots_adjust(wspace=0.07)
    axs[0].set_ylabel(r'PDF', fontsize=14)
    fig.supxlabel(r'SB$_{\rm intr}$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=14)
    #plot obs outside loops as it is a mess
    if msuf[m] == 'mlow':
        axs[1].set_title(r'$12.0 \leq \log M \leq 12.5$',fontsize=14)
        #MQN03
        obs1 = load(obs[1])
        axs[0].stairs(obs1['hist'],log10(10.**obs1['edges']*(1+3.3)**4.),fill=False,color='grey',ls='--')
        axs[1].stairs(obs1['hist'],log10(10.**obs1['edges']*(1+3.3)**4.),fill=False,color='grey',ls='--',label=r'MQN03')
        #MAGG
        obs2   = loadtxt(obs[2]).T
        edges2 = obs2[0]-0.07705
        edges2 = append(edges2,obs2[0,len(obs2[0])-1]+0.07705)
        axs[1].stairs(obs2[1],edges2,fill=False,color='black',ls='-.',label=r'MAGG')
        axs[2].stairs(obs2[1],edges2,fill=False,color='black',ls='-.')
    elif msuf[m] == 'mhigh':
        axs[1].set_title(r'$12.5 \leq \log M \leq 13.1$',fontsize=14)
        #MQN01
        obs3 = load(obs[0])
        axs[0].stairs(obs3['hist'],log10(10.**obs3['edges']*(1+3.2)**4.),fill=False,color='k',ls=':')
        axs[1].stairs(obs3['hist'],log10(10.**obs3['edges']*(1+3.2)**4.),fill=False,color='k',ls=':',label=r'MQN01')
    #sims
    for i in range(len(red)):
        axs[i].set_yscale('log')
        axs[i].set_ylim(0.002,4.0)
        xmin = log10(10.**(-18.3)*(1.+red[i])**4.)
        xmax = log10(10.**(-15.)*(1.+red[i])**4.)
        axs[i].set_xlim(xmin,xmax)
        axs[i].text(log10(10.**(-15.8)*(1.+red[0])**4.),2.0,r''+edges[0]+' $\leq r \leq$ '+edges[1],fontsize=12)
        axs[i].text(log10(10.**(-15.8)*(1.+red[0])**4.),3.0,r'z='+str(red[i]),fontsize=12)
        #sim loop
        for s in range(sg[0]):
            #for legends
            if i == 2:
                axs[i].plot([1,1],[-1,-1],color=colors[s],ls='-',label=r''+sims[s])
            
            if sims[s] == 'HELLO' and red[i] == 3.0:
                continue
            elif sims[s] == 'HELLO' and msuf[m] == 'mlow' and red[i] == 4.0:
                continue
            elif sims[s] == 'COLIBRE' and msuf[m] == 'mhigh' and red[i] == 4.0:
                continue
            else:
                data = load(paths[s]+'Hist_SB_s1_Ng256_observed_cubex_neb_singleobj_NB_z'+str(red[i])+'_grids'+grids[s][i][m][0]+'-'+grids[s][i][m][1]+'_R'+edges[0]+'-'+edges[1]+'.npz')
                H    = data['hist']
                E    = log10(10.**data['edges']*(1+red[0])**4.)
                axs[i].stairs(H,E,fill=False,ls='--',color=colors[s])
            
    axs[1].legend(loc='lower left',prop={'size':12},frameon=False)
    axs[2].legend(loc='lower left',prop={'size':12},frameon=False)
    plt.savefig(fout, bbox_inches='tight')
    print("Saved figure to "+fout)
    plt.show()
