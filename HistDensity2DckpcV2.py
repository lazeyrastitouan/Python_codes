from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import sys
from decimal import Decimal 
from scipy.optimize import curve_fit
from scipy.special import erf

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

def skewedgaussian(x,mu,sigma,alpha):
	'''
	Skewed aussian pdf for fit
	x = lnn
	'''
	fofx = 2. * gaussian(x,mu,sigma) * 0.5*( 1. + erf(alpha*(x-mu)/(sqrt(2.)*sigma)) )	
	return fofx

def mlp(x,mu,sigma,alpha):
    '''
    Modified lognormal power law (Basu, Gil, auddy 2015)
    '''
    return 0.5*alpha * exp(-alpha*(x-mu)+0.5*alpha**2.*sigma**2.) * ( 1. - erf((alpha*sigma-(x-mu)/sigma)/sqrt(2.)) )

# What to do
parser = argparse.ArgumentParser(description='Plot histogram of density, fit lognormal, and plot sigma and mu. Args: < path to density field file (str) >, < density field file suffix (str) >, < simulation type (str) >, < redshift (float) >, < mass (float) >')
parser.add_argument('filepath',nargs=1,type=str, help='Path to density field (str)')
parser.add_argument('suffix',nargs=1,type=str, help='Density field file suffix (the root is always density_volume_Tmin1e+3_Tmax1e+5_in_X-Xckpc_L), (str)')
parser.add_argument('simtype',nargs=1,type=str,help='Simulation type (dalya/colibre/hello) (str)')
parser.add_argument('redshift',nargs=1,type=float,help='Redshift (float)')
parser.add_argument('mass',nargs=1,type=float,help='log 10 of mass (float)')
args     = parser.parse_args()
filepath = args.filepath[0]
suffix   = args.suffix[0]
sim      = args.simtype[0]
red      = args.redshift[0]
mass     = args.mass[0]

plot_mu_sigma = False
fits          = [True,False] #skewedfit, mlp. We compute both fits in any case but plot only Trues 
nfits         = len(fits)
skewedfit     = sum(fits)
lsty          = [':','--']
models        = ['Skewed lognormal fit','MLP fit']
p0            =[[-1.,1.,1.],[-1.,1.,0.0]]
fRv           = [50,160,100,200] #[50,160,250,350,700] #ckpc
#if sim == 'hello' and suffix == 'M12p0_00128':
#	fRv   = [50,160,250,350]
if skewedfit:
	print("Fitting skewed lognormals")
	mu    = empty((nfits,len(fRv)-1)) 
	sigma = empty((nfits,len(fRv)-1))	
	alpha = empty((nfits,len(fRv)-1))
	funcs = [skewedgaussian,mlp]
else:
	sigma = empty(len(fRv)-1)
	mu    = empty(len(fRv)-1)
dmock         = linspace(-4.0,2.5,650)*log(10.)
Tmin          = 1e3
Tmax          = 1e5
nbin          = 65
edges         = linspace(-4.,2.5,nbin)*log(10.) #bins in ln as skewed PDF for ln n
midbins       = array([(edges[i]+edges[i+1])/2. for i in range(len(edges)-1)]) #bins cener in l #bins cener in ln 
#dfileroot     = 'densityS_Tmin1e+3_Tmax1e+5_in_'
dfileroot     = 'density_volume_Tmin1e+3_Tmax1e+5_in_'
pathout       = filepath
if sim == 'dalya':
	colors = ['darkblue','blue','cyan'] #['lightblue','cyan','blue','darkblue']
	# let's add H1312. The suffix and mass of H10 coincide with those of H10+H1312 combined. We just need to get the correct suffix to load H1312
	if 'zoom1p5' in filepath:
		zm='1p5'
		print('zoom'+zm)
		if red == 3.0:
			pathout   = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/DensityCGM/z3/zoom'+zm+'/'
			#only higher mass bin at z= 3 for H1312
			if mass == 12.7:
				suffix2   = 'M12p8'
				filepath2 = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/Halo1312_lvl7/sftresh_cool_low_fb_boost/z3/CDM/zoom'+zm+'/'
		elif red == 3.5:
			filepath2 = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/Halo1312_lvl7/sftresh_cool_low_fb_boost/z3p5/CDM/zoom'+zm+'/'
			pathout   = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom'+zm+'/'
			if mass == 12.6:
				suffix2 = 'M12p7'
			elif mass == 12.3:
				suffix2 = 'M12p4'
	elif 'zoom2' in filepath:
		zm='2'
		print('zoom'+zm)
		if red == 3.0:
			pathout   = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/DensityCGM/z3/zoom'+zm+'/'
			if mass == 12.6:
				suffix2   = 'M12p9'
				filepath2 = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/Halo1312_lvl7/sftresh_cool_low_fb_boost/z3/CDM/zoom'+zm+'/'
		elif red == 3.5:
			filepath2 = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/Halo1312_lvl7/sftresh_cool_low_fb_boost/z3p5/CDM/zoom'+zm+'/'	
			pathout   = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/gizmo/DensityCGM/z3p5/zoom'+zm+'/'
			if mass == 12.3:
				suffix2 = 'M12p4'
			
elif sim == 'colibre':
	colors = ['darkred','red','orange'] #['goldenrod','orange','red','darkred']
elif sim == 'hello':
	colors     = ['darkgreen','green','lime'] #['greenyellow','lime','mediumseagreen','darkgreen']
	pathout    = '/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/HELLO/CDM/'
	filepaths2 = ['/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/HELLO/CDM/g2.71e12/','/Users/titouanlazeyras/Desktop/Bicocca/Hydrosims_CDM_WDM/HELLO/CDM/g2.32e12/']

#save mu and sigma
if sim == 'dalya' and zm == '2' and mass == 12.6:
	if skewedfit: #at least one skewed fit
		f = open(pathout+'/mu_sigma_alpha_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin'+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_z"+str(red)+"_M12.8.txt",'w')
		fout = pathout+"/2D_VolumeWeighted_Density_Histogram_In_annulusckpc_skewedfit_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_12.8.pdf"
	else:
		f = open(pathout+'/mu_sigma_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin'+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_z"+str(red)+"_M12.8.txt",'w')
		fout = pathout+"/2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_12.8.pdf"
else:
	if skewedfit: #at least one skewed fit
		f = open(pathout+'/mu_sigma_alpha_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin'+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_z"+str(red)+"_M"+str(mass)+".txt",'w')
		fout = pathout+"/2D_VolumeWeighted_Density_Histogram_In_annulusckpc_skewedfit_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_"+str(mass)+".pdf"
	else:
		f = open(pathout+'/mu_sigma_2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin'+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_z"+str(red)+"_M"+str(mass)+".txt",'w')
		fout = pathout+"/2D_VolumeWeighted_Density_Histogram_In_annulusckpc_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+"_"+str(mass)+".pdf"

#figures	 
print("Evolution of pdf of projected gas with Rvir fraction")
plt.figure(1)
if mass == 12.0 or mass == 12.1 or mass == 12.3:
	plt.text(-1.85,0.2,r'$10^3 \le$ $T$ [K] $\le 10^5$, $\; 12.0 \leq \log M \leq 12.5$',fontsize=12)
else:
	plt.text(-1.85,0.2,r'$10^3 \le$ $T$ [K] $\le 10^5$, $\; 12.5 \leq \log M \leq 13.0$',fontsize=12)
plt.text(0.3,0.55,r'z = 3.5',fontsize=12)
if sim == 'dalya':
	plt.title(r'DaLya',fontsize=14)
	meanmedian = -1.83
if sim == 'hello':
	plt.title(r'HELLO',fontsize=14)
	meanmedian = -2.06
if sim == 'colibre':
	plt.title(r'COLIBRE',fontsize=14)
	meanmedian = -1.94
#plt.xlim(-3.6,1.0)
#plt.ylim(-5.0,1.0)
plt.xlabel(r'$\log (n)$ [cm$^{-3}$]',fontsize=14)
plt.ylabel(r'$\log(\textrm{PDF}(\rho_\textrm{gas}))$',fontsize=14)
for ir in range(len(fRv)-1):
	#trick to avoid plotting impossible r-bins
	if fRv[ir] >= fRv[ir+1]:
		continue
	print(str(fRv[ir])+" - "+str(fRv[ir+1])+" ckpc")
	#H10 or not dalya
	print("Loading "+filepath+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix+'.npz')
	DX = load(filepath+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix+'.npz')['density'] #ccm^-3
	DY = load(filepath+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Y_'+suffix+'.npz')['density']
	DZ = load(filepath+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Z_'+suffix+'.npz')['density']
	VX = load(filepath+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix+'.npz')['volume'] #pkpc^3
	VY = load(filepath+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Y_'+suffix+'.npz')['volume']
	VZ = load(filepath+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Z_'+suffix+'.npz')['volume']
	# dalya, add H1312
	if sim =='dalya' and not (mass == 12.3 and red == 3.0):
		print('Adding Halo1312')
		print("Loading "+filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix2+'.npz')
		DX1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix2+'.npz')['density'] #ccm^-3
		DY1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Y_'+suffix2+'.npz')['density'] #ccm^-3
		DZ1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Z_'+suffix2+'.npz')['density'] #ccm^-3
		DX  = concatenate((DX,DX1))
		DY  = concatenate((DY,DY1))
		DZ  = concatenate((DZ,DZ1))
		VX1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix2+'.npz')['volume'] #pkpc^3
		VY1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Y_'+suffix2+'.npz')['volume']
		VZ1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Z_'+suffix2+'.npz')['volume']
		VX  = concatenate((VX,VX1))
		VY  = concatenate((VY,VY1))
		VZ  = concatenate((VZ,VZ1))	
	# hello, add g2.32 and 2.71
	if sim == 'hello':
		for filepath2 in filepaths2:
			DX1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix+'.npz')['density'] #ccm^-3
			DY1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Y_'+suffix+'.npz')['density'] #ccm^-3
			DZ1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Z_'+suffix+'.npz')['density'] #ccm^-3
			DX  = concatenate((DX,DX1))
			DY  = concatenate((DY,DY1))
			DZ  = concatenate((DZ,DZ1))
			VX1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_X_'+suffix+'.npz')['volume'] #pkpc^3
			VY1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Y_'+suffix+'.npz')['volume']
			VZ1 = load(filepath2+"/"+dfileroot+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_Z_'+suffix+'.npz')['volume']
			VX  = concatenate((VX,VX1))
			VY  = concatenate((VY,VY1))
			VZ  = concatenate((VZ,VZ1))
			del DX1, DY1, DZ1, VX1, VY1, VZ1

	D = concatenate((DX,DY,DZ))
	V = concatenate((VX,VY,VZ)) 
	H = array(histogram(log(D),bins=edges,density=True,weights=V)[0])
	
	mean = average(midbins,weights=H)
	
	savez(pathout+"Hist_VolumeWeightedDensity_2D_"+str(fRv[ir])+"-"+str(fRv[ir+1])+'ckpc_z3p5_'+suffix,density=edges/log(10),count=H)	

	print("Mean density in cm^-3 = "+str(average(D,weights=V)))	
	
	plt.stairs(log10(H),edges/log(10.),color=colors[ir],fill=False,label=r''+str(fRv[ir])+r' - '+str(fRv[ir+1])+r' ckpc')	
	#if fRv[ir] == 100:
	#	plt.plot([mean/log(10.),mean/log(10.)],[-5,1],color=colors[ir],ls='--')
	#	plt.plot([meanmedian,meanmedian],[-5,1],color=colors[ir],ls='-')

	new_line=str(fRv[ir])+"\t"+str(fRv[ir+1])+"\t"

	#Fit (skewed) lognormal to get sigma
	if skewedfit:
		for i in range(nfits):
			func = funcs[i]
			popt, pcov  = curve_fit(func, midbins, H, p0=p0[i])
			Hfit        = func(dmock,popt[0],popt[1],popt[2])
			mu[i,ir]    = popt[0]
			sigma[i,ir] = popt[1]
			alpha[i,ir] = popt[2]
			print('mu '+str(popt[0]))
			print('sigma '+str(popt[1]))
			print("alpha "+str(popt[2]))
			new_line+=str(mu[i,ir])+"\t"+str(sigma[i,ir])+"\t"+str(alpha[i,ir])+"\t"			

			if fits[i]:
				plt.plot(dmock/log(10.),log10(Hfit),ls=lsty[i],color=colors[ir],marker='')
				if ir == len(fRv)-2:
					plt.plot([-1,0],[-10,-10],color='k',ls=lsty[i],label=r''+models[i])#for legend			
				if fRv[ir] == 160:
					plt.text(-3.3,-3.6,r'\underline{50 - 160}:',fontsize=12)
					plt.text(-3.3,-4.0,r'$\mu$ = '+str(round(mu[i,ir],1)),fontsize=12)
					plt.text(-3.3,-4.4,r'$\sigma$ = '+str(round(sigma[i,ir],1)),fontsize=12) 
					plt.text(-3.3,-4.8,r'$\alpha$ = '+str(round(alpha[i,ir],1)),fontsize=12) 				

			del popt, pcov, Hfit
	else:
		popt, pcov = curve_fit(gaussian, midbins, H, p0=[0.,1.])
		Hfit       = gaussian(dmock,popt[0],popt[1])
		mu[ir]    = popt[0]
		sigma[ir] = popt[1]
		print('mu '+str(popt[0]))
		print('sigma '+str(popt[1]))
		new_line+=str(mu[ir])+"\t"+str(sigma[ir])+"\t"
	
		plt.plot(dmock/log(10.),log10(Hfit),ls=':',color=colors[ir],marker='')

		del popt, pcov, Hfit
	
	new_line+=str(average(D,weights=V))+"\t"+str(average(log(D),weights=V))+"\n"
	f.write(new_line)

	del DX,DY,DZ,D,H

print("Path out "+pathout)

plt.legend(loc=2,ncols=2,columnspacing=1.5,frameon=False,prop={'size':12})
#plt.savefig(fout, bbox_inches='tight')
#print("Saved figure to "+fout)

f.close()

if plot_mu_sigma:
	fRvmid = array([105,205,300,525]) #0.1
	plt.figure(2)
	plt.plot(fRvmid,mu,'b.',ls='--',label=r'$\log M = $'+str(mass))
	#plt.text(0.7,-0.4,r'$z=$'+str(red),fontsize=14)
	plt.legend(loc=0,frameon=False,prop={'size':12})
	plt.xlabel(r'r [$R_\textrm{vir}$]',fontsize=14)
	plt.ylabel(r'$\log (\mu)$',fontsize=14)
	#plt.savefig(rootdalya+"/2D_Density_mu_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf", bbox_inches='tight')
	#print("Saved figure to "+rootdalya+"/2D_Density_mu_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf")

	plt.figure(3)
	plt.plot(fRvmid,sigma,'b.',ls='--',label=r'$\log M = $'+str(mass))
	#plt.text(0.15,1.95,r'$z=$'+str(red),fontsize=14)
	plt.legend(loc=0,frameon=False,prop={'size':12})
	plt.xlabel(r'r [$R_\textrm{vir}$]',fontsize=14)
	plt.ylabel(r'$\sigma$',fontsize=14)
	#plt.savefig(rootdalya+"/2D_Density_sigma_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf", bbox_inches='tight')
	#print("Saved figure to "+rootdalya+"/2D_Density_sigma_In_annulusRvir_Tmin"+'{:.0e}'.format(Decimal(Tmin))+"_Tmax"+'{:.0e}'.format(Decimal(Tmax))+".pdf")

	if skewedfit:
		plt.figure(4)
		plt.plot(fRvmid,alpha,'b.',ls='--',label=r'$\log M = $'+str(mass))
		#plt.text(0.15,1.95,r'$z=$'+str(red),fontsize=14)
		plt.legend(loc=0,frameon=False,prop={'size':12})
		plt.xlabel(r'r [$R_\textrm{vir}$]',fontsize=14)
		plt.ylabel(r'$\alpha$',fontsize=14)
plt.show()
