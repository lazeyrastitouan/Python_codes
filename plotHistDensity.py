import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
from scipy.special import erf

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

#function for median and percentiles
def quantile_from_hist(counts, bin_edges, q):
    """
    Estimate quantile q from a histogram.
    
    counts: array of bin counts or weights
    bin_edges: array of bin edges, length len(counts)+1
    q: quantile between 0 and 1
    """
    counts = np.asarray(counts, dtype=float)
    bin_edges = np.asarray(bin_edges, dtype=float)

    total = counts.sum()
    if total <= 0:
        raise ValueError("Histogram is empty.")

    cumulative = np.cumsum(counts)
    target = q * total

    k = np.searchsorted(cumulative, target)

    left = bin_edges[k]
    right = bin_edges[k + 1]
    bin_width = right - left

    cum_before = cumulative[k - 1] if k > 0 else 0.0
    count_in_bin = counts[k]

    if count_in_bin == 0:
        return left

    # Linear interpolation inside the bin:
    fraction = (target - cum_before) / count_in_bin
    return left + fraction * bin_width

def median_and_width_from_hist(counts, bin_edges):
    q16 = quantile_from_hist(counts, bin_edges, 0.16)
    q50 = quantile_from_hist(counts, bin_edges, 0.50)
    q84 = quantile_from_hist(counts, bin_edges, 0.84)

    central_width = q84 - q16
    sigma_like_width = 0.5 * central_width

    return q50, sigma_like_width, (q16, q84)


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

rootChi = '/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/COLIBRE/L0050N1504/z3p5/'
rootClo = '/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/COLIBRE/L0050N0752/z3p5/'
rootDhi = '/Users/titouanlazeyras/Desktop/BICOCCA/Hydrosims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom2/'
rootDlo = '/Users/titouanlazeyras/Desktop/BICOCCA/Hydrosims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom1p5/'
rootH   = '/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/HELLO/CDM/'
ln10 = np.log(10.)


# -------------------- Load mean -------------------- #
meanDhi, meanDlo, meanChi, meanClo, meanH = [],[],[],[],[] #50-160 then 100-200
annuli = ['50-160','100-200']#,'160-250']
for fi in range(len(annuli)):
	meanDlo.append(np.load(rootDlo+'Hist_VolumeWeightedDensity_2D_'+annuli[fi]+'ckpc_z3p5_M12p3.npz')) #log10(density), count
	meanDhi.append(np.load(rootDhi+'Hist_VolumeWeightedDensity_2D_'+annuli[fi]+'ckpc_z3p5_M12p3.npz')) #log10(density), count
	meanClo.append(np.load(rootClo+'Hist_VolumeWeightedDensity_2D_'+annuli[fi]+'ckpc_z3p5_M12p1.npz'))
	meanChi.append(np.load(rootChi+'Hist_VolumeWeightedDensity_2D_'+annuli[fi]+'ckpc_z3p5_M12p1.npz'))
	meanH.append(np.load(rootH+'Hist_VolumeWeightedDensity_2D_'+annuli[fi]+'ckpc_z3p5_M12p0_00128.npz'))

#Load median
medD = np.load(rootDlo+'Hist_VolumeWeightedMedians_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5.npz') #log10(density), count
#medC = np.load(rootC+'Hist_VolumeWeightedMedians_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5.npz')
#medH = np.load(rootH+'Hist_VolumeWeightedMedians_2D_100-200ckpc_z3p5_Mlow_12.0_Mhigh_12.5.npz')

#load parameters
paramsDlo = np.loadtxt(rootDlo+'mu_sigma_alpha_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin1e+3_Tmax1e+5_z3.5_M12.3.txt',usecols=(2,3,4,8)).T #mu, sigma, alpha, mean
paramsDhi = np.loadtxt(rootDhi+'mu_sigma_alpha_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin1e+3_Tmax1e+5_z3.5_M12.3.txt',usecols=(2,3,4,8)).T #mu, sigma, alpha, mean
paramsClo = np.loadtxt(rootClo+'mu_sigma_alpha_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin1e+3_Tmax1e+5_z3.5_M12.1.txt',usecols=(2,3,4,8)).T
paramsChi = np.loadtxt(rootChi+'mu_sigma_alpha_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin1e+3_Tmax1e+5_z3.5_M12.1.txt',usecols=(2,3,4,8)).T
paramsH = np.loadtxt(rootH+'mu_sigma_alpha_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin1e+3_Tmax1e+5_z3.5_M12.0.txt',usecols=(2,3,4,8)).T
mu      = [paramsDlo[0,0],paramsClo[0,0],paramsDhi[0,0],paramsChi[0,0],paramsH[0,0]]
sigma   = [paramsDlo[1,0],paramsClo[1,0],paramsDhi[1,0],paramsChi[1,0],paramsH[1,0]]
alpha   = [paramsDlo[2,0],paramsClo[2,0],paramsDhi[2,0],paramsChi[2,0],paramsH[2,0]]

#Lya/HeII data
lndmock  = np.linspace(-8,6,1400)
mudat    = [-5.3,-6.4,-5.7]  #[-5.3,-6.6,-6.0]
sigmadat = [2.6,2.9,2.75]#[3.0,3.3,3.17]     #min  3.0 3.17 
alphadat = [1.1,2.5,1.9]#[1.0,3.1,2.1]           #hires 1.9 2.4 2.1
#hdata   = skewedgaussian(lndmock,-1,2.2,2.5)
#hdata2  = skewedgaussian(lndmock,-4.5,3.0,2.5)


# -------------------- plot density pdfs -------------------- #
# DALYA fiducial only #
fig0, ax0 = plt.subplots(1,1)
colorD=['cyan','b']
ax0.set_yscale('log')
#mean
for i in range(len(annuli)):
        ax0.stairs(meanDlo[i]['count'],meanDlo[i]['density'],color=colorD[i],label=r'r='+annuli[i])
#medians
ax0.stairs(medD['count'],medD['density'],color='b',ls=':',label=r'Medians, r=100-200')
#obs
ax0.plot(lndmock/ln10,skewedgaussian(lndmock,-4.6,2.6,0),color='grey',ls='-',marker='',alpha=0.5,label=r'Observations, r=50-160')
ax0.errorbar([0,0.75,1.5],skewedgaussian(np.array([0,0.75,1.5])*ln10,-4.6,2.6,0),xerr=[0.3,0.3,0.3],yerr=0.0,xlolims=True,ls='',fmt='grey',alpha=0.5)
#fit
#fit = skewedgaussian(lndmock,mu[0],sigma[0],alpha[0])
#ax0.plot(lndmock/ln10,fit,color=colorD[0],ls='--',label=r'Skewed lognormal')

#for legend
#axs[0].plot([-6,-4],[10,10],'k:',label=r'Medians, R=100-200')
#axs[1].plot([-6,-4],[10,10],'k--',label=r'Skewed lognormal fit')

ax0.set_xlim(-3.2,2.0)
ax0.set_ylim(1e-5,1)

#median and sigma from histogram
medianDlo, widthDlo, intervalDlo = median_and_width_from_hist(meanDlo[0]['count'], meanDlo[0]['density']*ln10)
ax0.text(0.5,10**-0.3,r'\underline{Sim. 50 - 160}:',fontsize=13)
ax0.text(0.5,10**-0.65,r'$\log_{10}n_0$ = '+str(round(medianDlo/ln10,1)),fontsize=13)
ax0.text(0.5,10**-0.95,r'$\sigma$ = '+str(round(sigma[0],1)),fontsize=13)


#ax0.text(1,10**-0.3,r'\underline{50 - 160}:',fontsize=13)
#ax0.text(1,10**-0.6,r'$\mu$ = '+str(round(mu[0],1)),fontsize=13)
#ax0.text(1,10**-0.85,r'$\sigma$ = '+str(round(sigma[0],1)),fontsize=13)
#ax0.text(1,10**-1.1,r'$\alpha$ = '+str(round(alpha[0],1)),fontsize=13)


ax0.set_xlabel(r'$\log_{10} n_c$ [cm$^{-3}$]',fontsize=16)
ax0.set_ylabel(r'PDF',fontsize=16)
ax0.legend(loc='lower left',frameon=False,ncols = 1, prop={'size':14})
fout0="/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Hist_density_part_median_dalya1p5_z3p5_M12-12p5_50-200ckpc.pdf"
plt.savefig(fout0,bbox_inches='tight')
print("Saved figure to "+fout0)
#plt.show()

#---
# ALL SIMULATIONS #
fig, axs = plt.subplots(1,3,sharey=True,figsize=(17,4))
plt.subplots_adjust(wspace=0.07)
for ax in axs:
	ax.set_yscale('log')

colors = [['b','b'],['red','red'],['green','green']] #[['cyan','blue'],['orange','red'],['lime','green']]
titles= ['DaLya','COLIBRE','HELLO']

#plot mean
#for i in range(len(annuli)):
for i in [0]:
	axs[0].stairs(meanDlo[i]['count'],meanDlo[i]['density'],color=colors[0][i],label=r'Fiducial')
	axs[0].stairs(meanDhi[i]['count'],meanDhi[i]['density'],ls='--',color=colors[0][i],label=r'High res.')
	axs[1].stairs(meanClo[i]['count'],meanClo[i]['density'],color=colors[1][i],label=r'm6')	
	axs[1].stairs(meanChi[i]['count'],meanChi[i]['density'],ls='--',color=colors[1][i],label=r'm5')
	axs[2].stairs(meanH[i]['count'],meanH[i]['density'],ls='--',color=colors[2][i])
	#get medians and 1sigma error
	medianDhi, widthDhi, intervalDhi = median_and_width_from_hist(meanDhi[i]['count'], meanDhi[i]['density']*ln10)
	medianClo, widthClo, intervalClo = median_and_width_from_hist(meanClo[i]['count'], meanClo[i]['density']*ln10)
	medianChi, widthChi, intervalChi = median_and_width_from_hist(meanChi[i]['count'], meanChi[i]['density']*ln10)
	medianH, widthH, intervalH       = median_and_width_from_hist(meanH[i]['count'], meanH[i]['density']*ln10)

#plot medians
#axs[0].stairs(np.log10(medD['count']),medD['density']*ln10,color='r',ls=':')#,label=r'Medians, R=100-200')
#axs[1].stairs(np.log10(medC['count']),medC['density']*ln10,color='r',ls=':',label=r'Medians, R=100-200')
#axs[2].stairs(np.log10(medH['count']),medH['density']*ln10,color='r',ls=':')#,label=r'Medians')
#for legend
#axs[0].plot([-6,-4],[10,10],'k:',label=r'Medians, R=100-200')
#axs[1].plot([-6,-4],[10,10],'k--',label=r'Skewed lognormal fit')
for ix in range(len(axs)):
	ax = axs[ix]
	#create fits
	#fit = skewedgaussian(lndmock,mu[ix],sigma[ix],alpha[ix])
	if ix == 2:
		ax.plot(lndmock/ln10,skewedgaussian(lndmock,-4.6,2.6,0),color='grey',ls='-',marker='',alpha=0.5,label=r'Observations')
		#ax.plot(lndmock/ln10,skewedgaussian(lndmock,mudat[ix],sigmadat[ix],alphadat[ix]),color='grey',ls='-',marker='',alpha=0.5,label=r'Observations')
		#ax.plot(lndmock/ln10,fit,color=colors[ix][0],ls=':',label=r'Skewed lognormal')
		#ax.text(0.3, -0.1,r'$z=3.5$', fontsize=13)
	else:
		ax.plot(lndmock/ln10,skewedgaussian(lndmock,-4.6,2.6,0),color='grey',ls='-',marker='',alpha=0.5)
		#ax.plot(lndmock/ln10,fit,color=colors[ix][0],ls=':')
		#ax.plot(lndmock/ln10,skewedgaussian(lndmock,mudat[ix],sigmadat[ix],alphadat[ix]),color='grey',ls='-',marker='',alpha=0.5)
	#top = ax.secondary_xaxis("top",functions=(lambda x: x/ln10, lambda x: x*ln10))
	#top.set_xlabel(r'$\log_{10} n_c$ [cm$^{-3}$]',labelpad=15,fontsize=16)
	ax.errorbar([0,0.75,1.5],skewedgaussian(np.array([0,0.75,1.5])*ln10,-4.6,2.6,0),xerr=[0.3,0.3,0.3],yerr=0.0,xlolims=True,ls='',fmt='grey',alpha=0.5)
	ax.legend(loc='upper right',frameon=False,ncols = 1, prop={'size':13})
	ax.set_xlabel(r'$\log_{10} n_c$ [cm$^{-3}$]',fontsize=16)
	ax.set_title(r''+titles[ix],fontsize=18)
	ax.set_xlim(-3.2,2.0)
	ax.set_ylim(1e-5,2)
	ax.text(-3.1,0.8,r'50 $\leq$ r [ckpc] $\leq$ 160',fontsize=13)

#text
axs[0].text(-3.1,10**-3.8,r'\underline{Fiduc.}:',fontsize=13)
axs[1].text(-3.1,10**-3.8,r'\underline{m6}:',fontsize=13)
axs[0].text(-1,10**-3.8,r'\underline{High res.}:',fontsize=13)
axs[1].text(-1,10**-3.8,r'\underline{m5}:',fontsize=13)
axs[0].text(-3.1,5e-5,r'$\log_{10}n_0$ = '+str(round(medianDlo/ln10,1)),fontsize=13)
axs[0].text(-3.1,2e-5,r'$\sigma$ = '+str(round(sigma[0],1)),fontsize=13)
axs[0].text(-1,5e-5,r'$\log_{10}n_0$ = '+str(round(medianDhi/ln10,1)),fontsize=13)
axs[0].text(-1,2e-5,r'$\sigma$ = '+str(round(sigma[2],1)),fontsize=13)
axs[1].text(-3.1,5e-5,r'$\log_{10}n_0$ = '+str(round(medianClo/ln10,1)),fontsize=13)
axs[1].text(-3.1,2e-5,r'$\sigma$ = '+str(round(widthClo*4/3,1)),fontsize=13)
axs[1].text(-1,5e-5,r'$\log_{10}n_0$ = '+str(round(medianChi/ln10,1)),fontsize=13)
axs[1].text(-1,2e-5,r'$\sigma$ = '+str(round(widthChi*4/3,1)),fontsize=13)
axs[2].text(-3.1,5e-5,r'$\log_{10}n_0$ = '+str(round(medianH/ln10,1)),fontsize=13)
axs[2].text(-3.1,2e-5,r'$\sigma$ = '+str(round(sigma[4],1)),fontsize=13)

axs[0].set_ylabel(r'PDF',fontsize=16)
fout="/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Hist_density_part_dalya_colibre_hello_z3p5_M12-12p5_50-200ckpc_hires.pdf"
#plt.savefig(fout,bbox_inches='tight')
#print("Saved figure to "+fout)
#plt.show()

#--------- 
#PLOT WITH PDFS
fig2, axs2 = plt.subplots(2,1,sharey=True,figsize=(5,9))
plt.subplots_adjust(hspace=0.09)
for ax in axs2:
        ax.set_yscale('log')

colors  = ['b','r','b','r','g']
mapping = {0:0, 1:0, 2:1, 3:1, 4:1}

#plot mean
axs2[0].stairs(meanDlo[i]['count'],meanDlo[i]['density'],color='b',label=r'DaLya fiduc.')
axs2[1].stairs(meanDhi[i]['count'],meanDhi[i]['density'],color='b',label=r'DaLya HR')
axs2[0].stairs(meanClo[i]['count'],meanClo[i]['density'],color='r',label=r'COLIBRE m6')
axs2[1].stairs(meanChi[i]['count'],meanChi[i]['density'],color='r',label=r'COLIBRE m5')
axs2[1].stairs(meanH[i]['count'],meanH[i]['density'],color='g',label='HELLO')

#create fits
for ix in range(len(mu)):
	fit = skewedgaussian(lndmock,mu[ix],sigma[ix],alpha[ix])
	axs2[mapping[ix]].plot(lndmock/ln10,fit,color=colors[ix],ls=':')
#for legend
axs2[0].plot([-2,-1],[1e-6,1e-6],'k:',label=r'Skewed lognormal')
for ax in axs2:
	ax.legend(loc='upper right',frameon=False,ncols = 1, prop={'size':13})
	ax.set_ylabel(r'PDF',fontsize=16)
	ax.set_xlim(-3.2,2.0)
	ax.set_ylim(1e-5,2)
	ax.text(-3.1,2e-5,r'50 $\leq$ r [ckpc] $\leq$ 160',fontsize=13)
axs2[1].set_xlabel(r'$\log_{10} n_c$ [cm$^{-3}$]',fontsize=16)
#ax.set_title(r''+titles[ix],fontsize=18)
fout3="/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Hist_density_part_dalya_colibre_hello_z3p5_M12-12p5_50-200ckpc_hires_wFit.pdf"
#plt.savefig(fout3,bbox_inches='tight')
#print("Saved figure to "+fout3)
plt.show()
