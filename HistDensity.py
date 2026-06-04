from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import sys
from decimal import Decimal 
from scipy.optimize import curve_fit


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

def gaussian(x,mu,sigma):
	'''
	Gaussian pdf for fit
	x = lnn
	u = mean lnn
	'''
	return 1./(sigma*sqrt(2.*pi)) * exp(-((x-mu)**2./(2.*sigma**2.)))

plot_mu_sigma = True
fRv           = [0.0,0.2,0.5,0.7,1.0]
sigma12p3     = empty(len(fRv))
mu12p3        = empty(len(fRv))
sigma12p7     = empty(len(fRv))
mu12p7        = empty(len(fRv))

Tmin      = 1e3
Tmax      = 1e5
red       = 3
nbin      = 40
edges     = linspace(-3.1,2.1,nbin)*log(10.) #bins in ln
midbins   = array([(edges[i]+edges[i+1])/2. for i in range(len(edges)-1)]) #bins cener in l #bins cener in ln 
rootdalya = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/gizmo/Halo10_lvl7/sftresh_cool_low_fb_boost/z"+str(red)+"/CDM/zoom1p5/"
rooteagle = "/Users/titouanlazeyras/Desktop/BICOCCA/HydroSims_CDM_WDM/Eagle/RefL0100N1504/z"+str(red)+"/"
fouteagle = rooteagle+"/Density_Histogram_InOutAll_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_Mbins.pdf"
foutdalya = rootdalya+"/Density_Histogram_InOutAll_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_Mbins.pdf"
foutcomp  = rootdalya+"/Density_Histogram_In_eagle_dalya_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_Mbins.pdf"

# LOAD DaLya
#in
dInD12p3    = load(rootdalya+"/densityS_temperature_in_Hcat_12_12p5.npz")
dInD12p7    = load(rootdalya+"/densityS_temperature_in_Hcat_12p5_13p1.npz")
TinD12p3    = dInD12p3['temperature'] #K
TinD12p7    = dInD12p7['temperature'] #K
maskinD12p3 = where((TinD12p3 >= Tmin) & (TinD12p3 <= Tmax))
maskinD12p7 = where((TinD12p7 >= Tmin) & (TinD12p7 <= Tmax))
DinD12p3    = dInD12p3['density'][maskinD12p3] #ccm^-3
DinD12p7    = dInD12p7['density'][maskinD12p7] #ccm^-3
#out
#dOutD       = load(rootdalya+"/density_temperature_out.npz")
#TOutD       = dOutD['temperature'] #K
#maskoutD    = where((TOutD >= Tmin) & (TOutD <= Tmax))
#DoutD       = dOutD['density'][maskoutD] #ccm^-3
#all
#dAllD       = load(rootdalya+"/density_temperature_all.npz")
#TAllD       = dAllD['temperature'] #K
#maskallD    = where((TAllD >= Tmin) & (TAllD <= Tmax))
#DallD       = dAllD['density'][maskallD] #ccm^-3

#LOAD EAGLE
#in
#dInE12p3    = load(rooteagle+"/densityS_temperature_in_Hcat_12p25_12p5.npz")
#dInE12p8    = load(rooteagle+"/densityS_temperature_in_Hcat_12p75_13p25.npz")
#TinE12p3    = dInE12p3['temperature'] #K
#TinE12p8    = dInE12p8['temperature'] #K
#maskinE12p3 = where((TinE12p3 >= Tmin) & (TinE12p3 <= Tmax))
#maskinE12p8 = where((TinE12p8 >= Tmin) & (TinE12p8 <= Tmax))
#DinE12p3    = dInE12p3['density'][maskinE12p3] #ccm^-3
#DinE12p8    = dInE12p8['density'][maskinE12p8] #ccm^-3
#out
#dOutE       = load(rooteagle+"/densityS_temperature_out_Hcat.npz")
#TOutE       = dOutE['temperature'] #K
#maskoutE    = where((TOutE >= Tmin) & (TOutE <= Tmax))
#DoutE       = dOutE['density'][maskoutE] #ccm^-3
#all
#dAllE       = load(rooteagle+"/densityS_temperature_all_Hcat.npz")
#TAllE       = dAllE['temperature'] #K
#maskallE    = where((TAllE >= Tmin) & (TAllE <= Tmax))
#DallE       = dAllE['density'][maskallE] #ccm^-3

# HISTOGRAMS
HinD12p3 = array(histogram(log(DinD12p3),bins=edges,density=True)[0])
HinD12p7 = array(histogram(log(DinD12p7),bins=edges,density=True)[0])
#HoutD    = array(histogram(log10(DoutD),bins=edges,density=True)[0])
#HallD    = array(histogram(log10(DallD),bins=edges,density=True)[0])
#HinE12p3 = array(histogram(log10(DinE12p3),bins=edges,density=True)[0])
#HinE12p8 = array(histogram(log10(DinE12p8),bins=edges,density=True)[0])
#HoutE    = array(histogram(log10(DoutE),bins=edges,density=True)[0]) 
#HallE    = array(histogram(log10(DallE),bins=edges,density=True)[0])

# Fit lognormal to get sigma
print('Full range')
popt, pcov = curve_fit(gaussian, midbins, HinD12p3, p0=[0.,1.])
#print('12.3 log10(mu) '+str(popt[0]/log(10.)))
#print('12.3 sigma '+str(popt[1]))
mu12p3[4]    = popt[0]/log(10.)
sigma12p3[4] = popt[1]
dmock = linspace(-3.5,2.5,600)*log(10.)
HinD12p3fit = gaussian(dmock,popt[0],popt[1])
del popt, pcov

popt, pcov = curve_fit(gaussian, midbins, HinD12p7, p0=[0.,1.])
#print('12.7 log10(mu) '+str(popt[0]/log(10.)))
#print('12.7 sigma '+str(popt[1]))
mu12p7[4]    = popt[0]/log(10.)
sigma12p7[4] = popt[1]
HinD12p7fit = gaussian(dmock,popt[0],popt[1])
del popt, pcov

#figures 
# all, out and in (for each mass bin) for each simulation set
'''
print("Plotting for DaLya")
plt.figure(1)
plt.title(r'$10^3$K $\le$ T $\le$ $10^5$K',fontsize=12)
plt.text(-1.7,0.01,r'DaLya',fontsize=14)
plt.xlim(min(midbins),max(midbins))
plt.ylim(2e-5,2)
plt.yscale('log')
plt.xlabel(r'$\rho_\textrm{gas}$ [cm$^{-3}$]',fontsize=14)
plt.ylabel(r'PDF $\rho_\textrm{gas}$',fontsize=14)
plt.plot([-1,-1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.plot([0,0],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.plot([1,1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.stairs(HinD12p3,edges,color='blue',fill=False,label=r'Inside $R_\textrm{vir}$, $\log M = 12.3$')
plt.stairs(HinD12p7,edges,color='darkblue',fill=False,label=r'Inside $R_\textrm{vir}$, $\log M = 12.7$')
plt.stairs(HoutD,edges,color='cyan',fill=False,label=r'Outside $R_\textrm{vir}$')
plt.stairs(HallD,edges,color='lightblue',fill=False,label=r'All particles')
plt.legend(loc=0,prop={'size':14})
plt.savefig(foutdalya, bbox_inches='tight')
print("Saved figure to "+foutdalya)

print("Plotting for Eagle")
plt.figure(2)
plt.title(r'$10^3$K $\le$ T $\le$ $10^5$K',fontsize=12)
plt.text(-1.7,0.01,r'Eagle',fontsize=14)
plt.xlim(min(midbins),max(midbins))
plt.ylim(2e-5,2)
plt.yscale('log')
plt.xlabel(r'$\rho_\textrm{gas}$ [cm$^{-3}$]',fontsize=14)
plt.ylabel(r'PDF $\rho_\textrm{gas}$',fontsize=14)
plt.plot([-1,-1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.plot([0,0],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.plot([1,1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.stairs(HinE12p3,edges,color='red',fill=False,label=r'Inside $R_\textrm{vir}$, $\log M = 12.3$')
plt.stairs(HinE12p8,edges,color='darkred',fill=False,label=r'Inside $R_\textrm{vir}$, $\log M = 12.85$')
plt.stairs(HoutE,edges,color='orange',fill=False,label=r'Outside $R_\textrm{vir}$')
plt.stairs(HallE,edges,color='gold',fill=False,label=r'All particles')
plt.legend(loc=0,prop={'size':14})
plt.savefig(fouteagle, bbox_inches='tight')
print("Saved figure to "+fouteagle)

# comparison of eagle and dalya for "in" for 2 mass bins
print("Comparing Eagle and DaLya in and out")
plt.figure(3)
plt.title(r'$10^3$K $\le$ T $\le$ $10^5$K',fontsize=12)
plt.xlim(min(midbins),max(midbins))
plt.ylim(2e-5,2)
plt.yscale('log')
plt.xlabel(r'$\rho_\textrm{gas}$ [cm$^{-3}$]',fontsize=14)
plt.ylabel(r'PDF $\rho_\textrm{gas}$',fontsize=14)
plt.plot([-1,-1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.plot([0,0],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.plot([1,1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
plt.stairs(HinD12p3,edges,color='blue',fill=False,label=r'DaLya, in $R_\textrm{vir}$, $\log M = 12.3$')
plt.stairs(HinD12p7,edges,color='darkblue',fill=False,label=r'DaLya, in $R_\textrm{vir}$, $\log M = 12.7$')
plt.stairs(HinE12p3,edges,color='red',fill=False,label=r'Eagle, in $R_\textrm{vir}$, $\log M = 12.3$')
plt.stairs(HinE12p8,edges,color='darkred',fill=False,label=r'Eagle, in $R_\textrm{vir}$, $\log M = 12.85$')
plt.stairs(HoutD,edges,color='cyan',fill=False,label=r'DaLya, outside $R_\textrm{vir}$')
plt.stairs(HoutE,edges,color='orange',fill=False,label=r'Eagle, outside $R_\textrm{vir}$')
plt.legend(loc=0,prop={'size':12})
plt.savefig(foutcomp, bbox_inches='tight')
print("Saved figure to "+foutcomp)
'''

print("Dalya evolution of in pdf with Rvir fraction")
colorsd   = ['lightblue','cyan','blue','darkblue']
foutdalyaRvirM12p3 = rootdalya+"/Density_Histogram_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_12p3.pdf"
plt.figure(4)
plt.title(r'$10^3$K $\le$ T $\le$ $10^5$K',fontsize=12)
plt.text(-2.9,0.27,r'DaLya',fontsize=14)
plt.text(-2.9,0.25,r'$\log M = 12.3$',fontsize=14)
plt.xlim(-3,2)
plt.ylim(0,0.38)
plt.xlabel(r'$\log (\rho_\textrm{gas})$ [cm$^{-3}$]',fontsize=14)
plt.ylabel(r'PDF $\rho_\textrm{gas}$',fontsize=14)
plt.plot([-1,-1],[1e-5,4],color='grey',ls='--', linewidth=0.7)
plt.plot([0,0],[1e-5,4],color='grey',ls='--', linewidth=0.7)
plt.plot([1,1],[1e-5,4],color='grey',ls='--', linewidth=0.7)
plt.plot(dmock/log(10.),HinD12p3fit,'k:',label=r'Lognormal fit')
plt.stairs(HinD12p3,edges/log(10.),color='black',fill=False,label=r'Full range')
for ir in range(4):
	#print(str(fRv[ir])+"-"+str(fRv[ir+1]))
	data = load(rootdalya+"/densityS_temperature_in_"+str(fRv[ir])+"-"+str(fRv[ir+1])+"Rvir_Hcat_12_12p5.npz")
	T    = data['temperature'] #K
	mask = where((T >= Tmin) & (T <= Tmax))
	D    = data['density'][mask] #ccm^-3
	H    = array(histogram(log(D),bins=edges,density=True)[0])
	
	#Fit lognormal to get sigma
	popt, pcov = curve_fit(gaussian, midbins, H, p0=[0.,1.])
	#print('12.3 log10(mu) '+str(popt[0]/log(10.)))
	#print('12.3 sigma '+str(popt[1]))
	Hfit = gaussian(dmock,popt[0],popt[1])

	mu12p3[ir]    = popt[0]/log(10.)
	sigma12p3[ir] = popt[1]

	plt.stairs(H,edges/log(10.),color=colorsd[ir],fill=False,label=r''+str(fRv[ir])+r' - '+str(fRv[ir+1])+r' $R_\textrm{vir}$')
	plt.plot(dmock/log(10.),Hfit,ls=':',color=colorsd[ir],marker='')

	del data,T,D,H, Hfit, popt, pcov
plt.legend(loc=2,ncols=3,columnspacing=1.5,frameon=False,prop={'size':12})
#plt.savefig(foutdalyaRvirM12p3, bbox_inches='tight')
#print("Saved figure to "+foutdalyaRvirM12p3)

foutdalyaRvirM12p7 = rootdalya+"/Density_Histogram_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_12p7.pdf"
plt.figure(5)
plt.title(r'$10^3$K $\le$ T $\le$ $10^5$K',fontsize=12)
plt.text(-2.9,0.27,r'DaLya',fontsize=14)
plt.text(-2.9,0.25,r'$\log M = 12.7$',fontsize=14)
plt.xlim(-3,2)
plt.ylim(0.,0.38)
#plt.yscale('log')
plt.xlabel(r'$\log (\rho_\textrm{gas})$ [cm$^{-3}$]',fontsize=14)
plt.ylabel(r'PDF $\rho_\textrm{gas}$',fontsize=14)
plt.plot([-1,-1],[1e-5,4],color='grey',ls='--', linewidth=0.7)
plt.plot([0,0],[1e-5,4],color='grey',ls='--', linewidth=0.7)
plt.plot([1,1],[1e-5,4],color='grey',ls='--', linewidth=0.7)
plt.plot(dmock/log(10.),HinD12p7fit,'k:',label=r'Lognormal fit')
plt.stairs(HinD12p7,edges/log(10.),color='black',fill=False,label=r'Full range')
for ir in range(4):
	#print(str(fRv[ir])+"-"+str(fRv[ir+1]))
	data = load(rootdalya+"/densityS_temperature_in_"+str(fRv[ir])+"-"+str(fRv[ir+1])+"Rvir_Hcat_12p5_13p1.npz")
	T    = data['temperature'] #K
	mask = where((T >= Tmin) & (T <= Tmax))
	D    = data['density'][mask] #ccm^-3
	H    = array(histogram(log(D),bins=edges,density=True)[0])
	# lognormal fit to get sigma
	popt, pcov = curve_fit(gaussian, midbins, H, p0=[0.,1.])
	#print('12.7 log10(mu) '+str(popt[0]/log(10.)))
	#print('12.7 sigma '+str(popt[1]))
	Hfit = gaussian(dmock,popt[0],popt[1])
	
	mu12p7[ir]    = popt[0]/log(10.)
	sigma12p7[ir] = popt[1]

	plt.plot(dmock/log(10.),Hfit,ls=':',color=colorsd[ir],marker='')
	plt.stairs(H,edges/log(10.),color=colorsd[ir],fill=False,label=r''+str(fRv[ir])+r'$ - $'+str(fRv[ir+1])+r' $R_\textrm{vir}$')
	
	del data,T,D,H, Hfit, popt, pcov
plt.legend(loc=2,ncols=3,columnspacing=1.5,frameon=False,prop={'size':12})
#plt.savefig(foutdalyaRvirM12p7, bbox_inches='tight')
#print("Saved figure to "+foutdalyaRvirM12p7)

if plot_mu_sigma:
	R12p3 = 98.0
	R12p7 = 130.0 
	fRvmid = array([0.1,0.35,0.6,0.85,1.0])
	plt.figure(1)
	plt.plot(fRvmid,mu12p3,'b.',ls='--',label=r'$\log M = 12.3$')
	plt.plot(fRvmid,mu12p7,'r.',ls='--',label=r'$\log M = 12.7$')
	plt.text(100.,-0.3,r'$z=3$',fontsize=14)
	plt.legend(loc=0,frameon=False,prop={'size':12})
	plt.xlabel(r'r [$R_\textrm{vir}$]',fontsize=14)
	plt.ylabel(r'$\log (\mu)$',fontsize=14)
	plt.savefig(rootdalya+"/Density_mu_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf", bbox_inches='tight')
	print("Saved figure to "+rootdalya+"/Density_mu_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf")

	plt.figure(2)
	plt.plot(fRvmid,sigma12p3,'b.',ls='--',label=r'$\log M = 12.3$')
	plt.plot(fRvmid,sigma12p7,'r.',ls='--',label=r'$\log M = 12.7$')
	plt.text(15.,1.7,r'$z=3$',fontsize=14)
	plt.legend(loc=0,frameon=False,prop={'size':12})
	plt.xlabel(r'r [$R_\textrm{vir}$]',fontsize=14)
	plt.ylabel(r'$\log (\sigma)$',fontsize=14)
	plt.savefig(rootdalya+"/Density_sigma_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf", bbox_inches='tight')	
	print("Saved figure to "+rootdalya+"/Density_sigma_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf")

#print("Eagle evolution of in pdf with Rvir fraction")
#colorse=['goldenrod','orange','red','darkred']
#fouteagleRvirM12p3 = rooteagle+"/Density_Histogram_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_12p3.pdf"
#plt.figure(6)
#plt.title(r'$10^3$K $\le$ T $\le$ $10^5$K',fontsize=12)
#plt.text(-1.7,0.01,r'Eagle',fontsize=14)
#plt.text(-1.7,0.005,r'$\log M = 12.3$',fontsize=14)
#plt.xlim(min(midbins),max(midbins))
#plt.ylim(2e-4,2)
#plt.yscale('log')
#plt.xlabel(r'$\rho_\textrm{gas}$ [cm$^{-3}$]',fontsize=14)
#plt.ylabel(r'PDF $\rho_\textrm{gas}$',fontsize=14)
#plt.plot([-1,-1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
#plt.plot([0,0],[1e-5,2],color='grey',ls='--', linewidth=0.7)
#plt.plot([1,1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
#for ir in range(4):
#        data = load(rooteagle+"/densityS_temperature_in_"+str(fRv[ir])+"-"+str(fRv[ir+1])+"Rvir_Hcat_12p25_12p5.npz")
#        T    = data['temperature'] #K
#        mask = where((T >= Tmin) & (T <= Tmax))
#        D    = data['density'][mask] #ccm^-3
#        H    = array(histogram(log10(D),bins=edges,density=True)[0])
#        plt.stairs(H,edges,color=colorse[ir],fill=False,label=r''+str(fRv[ir])+r'$ - $'+str(fRv[ir+1])+r' $R_\textrm{vir}$')
#        del data,T,D,H
#plt.stairs(HinE12p3,edges,color='darkred',fill=False,label=r'$< R_\textrm{vir}$')
#plt.legend(loc=0,prop={'size':12})
#plt.savefig(fouteagleRvirM12p3, bbox_inches='tight')
#print("Saved figure to "+fouteagleRvirM12p3)

#fouteagleRvirM12p8 = rooteagle+"/Density_Histogram_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_12p8.pdf"
#plt.figure(7)
#plt.title(r'$10^3$K $\le$ T $\le$ $10^5$K',fontsize=12)
#plt.text(-1.7,0.01,r'Eagle',fontsize=14)
#plt.text(-1.7,0.005,r'$\log M = 12.8$',fontsize=14)
#plt.xlim(min(midbins),max(midbins))
#plt.ylim(2e-4,2)
#plt.yscale('log')
#plt.xlabel(r'$\rho_\textrm{gas}$ [cm$^{-3}$]',fontsize=14)
#plt.ylabel(r'PDF $\rho_\textrm{gas}$',fontsize=14)
#plt.plot([-1,-1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
#plt.plot([0,0],[1e-5,2],color='grey',ls='--', linewidth=0.7)
#plt.plot([1,1],[1e-5,2],color='grey',ls='--', linewidth=0.7)
#for ir in range(4):
#        data = load(rooteagle+"/densityS_temperature_in_"+str(fRv[ir])+"-"+str(fRv[ir+1])+"Rvir_Hcat_12p75_13p25.npz")
#        T    = data['temperature'] #K
#        mask = where((T >= Tmin) & (T <= Tmax))
#        D    = data['density'][mask] #ccm^-3
#        H    = array(histogram(log10(D),bins=edges,density=True)[0])
#        plt.stairs(H,edges,color=colorse[ir],fill=False,label=r''+str(fRv[ir])+r'$ - $'+str(fRv[ir+1])+r' $R_\textrm{vir}$')
#        del data,T,D,H
#plt.stairs(HinE12p8,edges,color='darkred',fill=False,label=r'$< R_\textrm{vir}$')
#plt.legend(loc=0,prop={'size':12})
#plt.savefig(fouteagleRvirM12p8, bbox_inches='tight')
#print("Saved figure to "+fouteagleRvirM12p8)

plt.show()
