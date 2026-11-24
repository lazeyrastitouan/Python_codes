import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
from scipy.special import erf

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

#functions for fits
def gaussian(x,mu,sigma):
        '''
        Gaussian pdf for fit
        x = lnn
        u = mean lnn
        '''
        return 1./(sigma*np.sqrt(2.*np.pi)) * np.exp(-((x-mu)**2./(2.*sigma**2.)))

def skewedgaussian(x,mu,sigma,alpha):
        '''
        Skewed aussian pdf for fit
        x = lnn
        '''
        fofx = 2. * gaussian(x,mu,sigma) * 0.5*( 1. + erf(alpha*(x-mu)/(np.sqrt(2.)*sigma)) )
        return fofx

#Number of halos: HELLO [1,1], Dalya [10,5], COLIBRE [12,3]
saveHist = True
nbin     = 65
edges    = np.linspace(-4.,2.5,nbin)*np.log(10.) #bins in ln
binwidth = edges[1]-edges[0] #in ln
mind     = -4.0*np.log(10.)
maxd     = 2.5*np.log(10.)
maxdh1   = -0.5*np.log(10.)
maxdh2   = 0.0*np.log(10.)

#load median densities and stds
Hfile_std = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/g9.20e12/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['std']
Hfile_med = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/g9.20e12/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['median']
Hfile_ave = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/g9.20e12/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['mean']
Dfile_std = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['std']
Dfile_med = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['median']
Dfile_ave = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['mean']
Dfile_asq = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/medians_means_meanssq_stds_Mlow_12.0_Mhigh_12.5.npz")['mean_sq']
Cfile_std = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/COLIBRE/L0050N0752/z3p5/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['std']
Cfile_med = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/COLIBRE/L0050N0752/z3p5/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['median']
Cfile_ave = np.load("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/COLIBRE/L0050N0752/z3p5/medians_stds_Mlow_12.0_Mhigh_12.5.npz")['mean']

#when there is only 1 particle in a cell, std = 0 and log10 won't work so replace 0 with minimum float32
Hfile_std[Hfile_std==0.0]=1e-37
Dfile_std[Dfile_std==0.0]=1e-37
Cfile_std[Cfile_std==0.0]=1e-37

#compute histograms 
hist_med_hlo  = np.histogram(np.log(Hfile_med),density=True,bins=edges,range=(mind,maxdh1))
hist_std_hlo  = np.histogram(np.log(Hfile_std),density=True,bins=100)
hist_ave_hlo  = np.histogram(np.log(Hfile_ave),density=True,bins=edges,range=(mind,maxdh2))
hist_med_dlya = np.histogram(np.log(Dfile_med),density=True,bins=edges,range=(mind,maxd))
hist_std_dlya = np.histogram(np.log(Dfile_std),density=True,bins=100)
hist_ave_dlya = np.histogram(np.log(Dfile_ave),density=True,bins=edges,range=(mind,maxd))
hist_asq_dlya = np.hostpgram(np.log(Dfile_asq),density=True,bins=2*edges,range=(2*mind,2*maxd))
hist_med_col  = np.histogram(np.log(Cfile_med),density=True,bins=edges,range=(mind,maxd))
hist_std_col  = np.histogram(np.log(Cfile_std),density=True,bins=100)
hist_ave_col  = np.histogram(np.log(Cfile_ave),density=True,bins=edges,range=(mind,maxd))

#fits
nbins        = len(hist_med_dlya[0])
nbinsh1      = len(hist_med_hlo[0])
nbinsh2      = len(hist_ave_hlo[0])
dmock        = np.linspace(-3.5,0.5,400)*np.log(10.)
midbins      = hist_med_dlya[1][0:nbins]+binwidth
midbinsh1    = hist_med_hlo[1][0:nbins]+binwidth
midbinsh2    = hist_ave_hlo[1][0:nbins]+binwidth
poptD, pcovD = curve_fit(skewedgaussian, midbins, hist_med_dlya[0], p0=[-1.,1.,-1.])
HfitD        = skewedgaussian(dmock,poptD[0],poptD[1],poptD[2])

poptC, pcovC = curve_fit(skewedgaussian, midbins, hist_med_col[0], p0=[-1.,1.,-1.])
HfitC        = skewedgaussian(dmock,poptC[0],poptC[1],poptC[2])

poptH, pcovH = curve_fit(skewedgaussian, midbinsh2, hist_med_hlo[0], p0=[-1.,1.,-1.])
HfitH        = skewedgaussian(dmock,poptH[0],poptH[1],poptH[2])

if saveHist:
	np.savez("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/Hist_VolumeWeightedMedians_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5",density=hist_med_hlo[1]/np.log(10.),count=hist_med_hlo[0])
	np.savez("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/COLIBRE/L0050N0752/z3p5/Hist_VolumeWeightedMedians_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5",density=hist_med_col[1]/np.log(10.),count=hist_med_col[0])
	np.savez("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/Hist_VolumeWeightedMedians_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5",density=hist_med_dlya[1]/np.log(10.),count=hist_med_dlya[0])
    np.savez("/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/Hist_VolumeWeightedMeans_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5",density=hist_ave_dlya[1]/np.log(10.),count=hist_ave_dlya[0],density_sq=hist_asq_dlya[1]/np.log(10.),count_sq=hist_asq_dlya[0])
	print("Saved dalya to /Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/Hist_VolumeWeightedMedians_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5.npz")

print(np.average(midbins,weights=hist_med_dlya[0])/np.log(10.))
print(np.average(midbins,weights=hist_med_col[0])/np.log(10.))
print(np.average(midbinsh1,weights=hist_med_hlo[0])/np.log(10.))

#plot
fig, axs = plt.subplots(1,3,sharey=True,figsize=(16,4))
plt.subplots_adjust(wspace=0.07)
titles = ['DaLya','COLIBRE','HELLO']
for ix in range(len(axs)):
	axs[ix].set_yscale('log')
	axs[ix].set_title(r''+titles[ix],fontsize=14)
	axs[ix].set_xlabel(r'$\log (\rho_g)$ [atoms cm$^{-3}$]',fontsize = 12)
	axs[ix].set_xlim(-3.,0.5)
	axs[ix].set_ylim(1e-5,2.5)
axs[2].plot(dmock/np.log(10.),HfitH,color='g',ls=':') 
axs[1].plot(dmock/np.log(10.),HfitC,color='r',ls=':')
axs[0].plot(dmock/np.log(10.),HfitD,color='b',ls=':')
axs[2].stairs(hist_med_hlo[0],hist_med_hlo[1]/np.log(10.),color='g')
axs[0].stairs(hist_med_dlya[0],hist_med_dlya[1]/np.log(10.),color='b')
axs[1].stairs(hist_med_col[0],hist_med_col[1]/np.log(10.),color='r')
axs[2].stairs(hist_ave_hlo[0],hist_ave_hlo[1]/np.log(10.),color='g',ls='--')
axs[0].stairs(hist_ave_dlya[0],hist_ave_dlya[1]/np.log(10.),color='b',ls='--')
axs[1].stairs(hist_ave_col[0],hist_ave_col[1]/np.log(10.),color='r',ls='--')
axs[1].text(-2.8,1.2,r'z=3.5, $\; \; 12.0 \leq \log M \leq 12.5$, $\; \; 100 \leq R \leq 200$',fontsize=11)
axs[0].text(-2.8,1.2,r'$10^3 \leq T \, [\textrm{K}] \leq 10^5$',fontsize=11) 
axs[0].text(-2.8,0.005,r'$\mu$ = '+str(round(poptD[0],1)),fontsize=11)
axs[0].text(-2.8,0.003,r'$\sigma$ = '+str(round(poptD[1],1)),fontsize=11)
axs[0].text(-2.8,0.002,r'$\alpha$ = '+str(round(poptD[2],1)),fontsize=11)
axs[1].text(-2.8,0.005,r'$\mu$ = '+str(round(poptC[0],1)),fontsize=11)
axs[1].text(-2.8,0.003,r'$\sigma$ = '+str(round(poptC[1],1)),fontsize=11)
axs[1].text(-2.8,0.002,r'$\alpha$ = '+str(round(poptC[2],1)),fontsize=11)
axs[2].text(-2.8,0.005,r'$\mu$ = '+str(round(poptH[0],1)),fontsize=11)
axs[2].text(-2.8,0.003,r'$\sigma$ = '+str(round(poptH[1],1)),fontsize=11)
axs[2].text(-2.8,0.002,r'$\alpha$ = '+str(round(poptH[2],1)),fontsize=11)
axs[2].plot([10,11],[30,30],'k:',label=r'Skewed lognormal fit')
axs[2].plot([10,11],[30,30],'k-',label=r'Median')
axs[2].plot([10,11],[30,30],'k--',label=r'Mean')
axs[2].legend(loc='upper right',frameon=False,prop={'size':11})
axs[0].set_ylabel(r'PDF',fontsize = 12)
fout = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Hist_medians_mean_density_Mlow_12.0_Mhigh_12.5_z3p5_dalya_zoom1p5_colibre_m6_hello.pdf"
#plt.savefig(fout, bbox_inches='tight')
print("Saved figure to "+fout)

plt.show()

'''
plt.figure(2)
plt.stairs(hist_std_hlo[0],hist_std_hlo[1],color='g',label=r'HELLO')
plt.stairs(hist_std_dlya[0],hist_std_dlya[1],color='b',label=r'DaLya')
plt.stairs(hist_std_col[0],hist_std_col[1],color='r',label=r'COLIBRE')
plt.xlabel(r'std dev $\rho_g$ [atoms cm$^{-3}$]',fontsize = 12)
plt.ylabel(r'PDF',fontsize=12)
plt.xlim(-10,1.5)
plt.ylim(0,0.35)
plt.text(-9.5,0.28,r'z=3.5',fontsize=11)
plt.text(-9.5,0.32,r'$12.0 \leq \log M \leq 12.5$',fontsize=11)
plt.text(-9.5,0.3,r'$100 \leq R \leq 200$',fontsize=11)
plt.legend(loc='upper right',frameon=False,prop={'size':10})
fout2 = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Hist_std_density_Mlow_12.0_Mhigh_12.5_z3p5_dalya_zoom1p5_colibre_m6_hello.pdf"
plt.savefig(fout2, bbox_inches='tight')
plt.show()
'''
