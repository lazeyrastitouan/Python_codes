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
sims   = [r'No ISM cut', r'$n_\star=10$ cm$^{-3}$', r'Dust $M_\textrm{frac}<10^{-5}$'] #['Density cut','No cut','Metal cut','Dust cut']
grid   = ['256','256','256']
smooth = ['1.5']
slbl   = ['1p5']
grids  = [ [['4','14']],[['4','14']],[['4','14']]]#,[['4','15']],[['4','15']]]
Ngrids = [ [33] , [33], [33]]#, [36], [36]]
sg     = shape(grids) #num sims, num mbin, low grid - high grid
colors = [['m','r','orange']]#,'m','orange']]#,['cyan','orange','lime']]
lsty   = ['-','-','-']#,':',':']
msuf   = ['mlow']
paths  = ['/data/titouan/COLIBRE/L0050N1504/z3p5/colibre_0060/2cMpc/noISMcut/','/data/titouan/COLIBRE/L0050N1504/','/data/titouan/COLIBRE/L0050N1504/z3p5/colibre_0060/2cMpc/dust-5/']#'/data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/2cMpc/noISMcut/','/data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/2cMpc/metal005/','/data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/2cMpc/dust-5/']

binwidth = 0.1 #in log10
obs      = loadtxt('/data/titouan/MQN/SB_PDF_100_200ckpc_densityFalse.txt').T
edges2   = obs[0]-binwidth/2.
edges2   = append(edges2,obs[0,len(obs[0])-1]+binwidth/2.)
Nneb     = 14 #*27678/19575.

means_100_200 = [1.0940971498861521e-14,2.9482917005876948e-15] 

#mass loop
for m in range(sg[1]):
    for sm in range(len(smooth)):
        fig, axs = plt.subplots(1,1)
        fout = '/data/titouan/COLIBRE/Hist_SB_Npix_R100-200_colibre_m5_'+msuf[m]+'_z'+str(red)+'_n5e-20_nonoise_smooth'+slbl[sm]+'_impactISMcut.pdf'
        axs.set_ylabel(r'$\textrm{d} N_{pix} / \textrm{dlog} SB$', fontsize=18)
        #if msuf[m] == 'mlow':
        #    axs.set_title(r'$12.0 \leq \log M \leq 12.5$',fontsize=14)
        #elif msuf[m] == 'mhigh':
        #    axs.set_title(r'$12.5 \leq \log M \leq 13.1$',fontsize=14)
        axs.set_yscale('log')
        axs.set_xscale('log')
        axs.set_ylim(1,3e3) #1e3
        axs.set_xlim(3e-20,3e-12)
        #axs.text(-17.85,125,r''+edges[0]+' $\leq r \leq$ '+edges[1],fontsize=11)
        #axs.text(4e-18,1000,r'$n_\star=10$ cm$^{-3}$',fontsize=16)
        axs.set_xlabel(r'$\log$(SB$_{\rm intr}$) [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]', fontsize=18)
        #plot obs
        #if sm ==0:
            #axs.axvline(1e-16, 0, 1, color='grey', linestyle=':')
            #axs.axvline(3.4465102094748147e-15,0, 1, color='black', linestyle='--')
            #axs.text(-17.85,300,r'$n_\star = 10$',fontsize=13)
            #axs.text(-14.65,300,r'FWHM = 27.5 ckpc = 0.8"',fontsize=12)
        axs.stairs(obs[1]/Nneb/binwidth,10**edges2,fill=False,color='black',ls='-',label=r'Observations')
        #plot sims
        for s in range(len(sims)):
            data = load(paths[s]+'Hist_SB_Npixel_n5e-20_s'+smooth[sm]+'_Ng'+grid[s]+'_observed_cubex_neb_singleobj_NB_nonoise_z'+str(red)+'_grids'+grids[s][m][0]+'-'+grids[s][m][1]+'_R'+edges[0]+'-'+edges[1]+'.npz')
            H    = data['hist']/Ngrids[s][m]/binwidth
            E    = 10.**data['edges']*(1.+red)**4.
            axs.stairs(H,E,fill=False,ls=lsty[s],color=colors[sm][s],label=r''+sims[s])#+' '+slbl[sm])
            #if m==1 and sm ==1:
            #axs.axvline(means_100_200[s], 0, 1, color=colors[sm][s], linestyle='--')

        axs.legend(loc='upper left',prop={'size':12},ncol=1,frameon=False)
        axs.set_title(r"COLIBRE m5", fontsize =16)

        plt.savefig(fout, bbox_inches='tight')
        print("Saved figure to "+fout)
        plt.show()
