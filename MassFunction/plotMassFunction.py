#plot Tinker 08 and data points + lower panel ratio

from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)

######################################################################################
#cosmo
h = 0.7
rho_crit = 27.75e10 #h^2 Msun/Mpc^3
Omega_m = 0.27 #0.3 #a modifier suivant cosmo
rho_m = rho_crit*Omega_m #h^2 Msun/Mpc^3

#load Tinker for continuous curve
tink=loadtxt("Sigma_TinkerMassFunction_OfM") #M^2/rho_m*dndM
Mtink=tink[:,0]
dndlnMtink=tink[:,2]/Mtink*rho_m

#load data
data=loadtxt("500h_512/d0.00/200_snap_z0.000_average_v15") #M^2/rho_m*dndM
M = 10.0**data[:,0]
Mlog=log10(M) 
dndlnM=data[:,1]/M*rho_m
sdndlnM=data[:,2]/M*rho_m

######################################################################################
#compute tinker for M of data
#definition des parametres (Tinker et al (2008)) + cosmo fiducielle
A = 0.186
a = 1.47
b = 2.57
c = 1.19
R = 8.0 #Mpc/h

#charger le PS et creer la masse
PStxt = loadtxt("Pk_lcdm_transfer_out.dat_s8_0.7913_ns_0.9500") #a modifier suivant cosmo
k = PStxt[:,0] #h/Mpc
PS = PStxt[:,1] #(Mpc/h)^3
		
R_M = (3.0*M/(4.0*pi*rho_m))**(1.0/3.0) #Mpc/h - filtre tophat

Sk = shape(k)[0]
SM = shape(M)[0]

#calcul de la mass function pour tout R
#filtre tophat pour tout R
Cc = 4.0*pi*rho_m/3.0
W_M = zeros((SM, Sk))
for i in range(SM) :
	for l in range(Sk) :
		W_M[i,l] = 3.0/(k[l]*R_M[i])**3.0*(sin(k[l]*R_M[i])-k[l]*R_M[i]*cos(k[l]*R_M[i]))

#calcul de sigma pour tous les R_M
s_m_sq = zeros(SM)
d_s_m = zeros(SM)
s_m = zeros(SM)
f_s_m = zeros(SM)

for i in range(SM) : 
	s = 0
	while s < Sk-1 :
		d_s_m[i] = (k[s+1]-k[s])*1.0/(2.0*(pi)**2.0)*(W_M[i,s]**2.0*PS[s]*k[s]**2.0+W_M[i,s+1]**2.0*PS[s+1]*k[s+1]**2.0)/2.0
		s_m_sq[i] = s_m_sq[i] + d_s_m[i]
		s=s+1
#sigma
s_m = sqrt(s_m_sq)
f_s_m = A*((s_m/b)**(-a)+1.0)*exp(-c/s_m_sq) 
	
#dsigma/dM
#refaire une integrale pour d_sigma^2/d_M
d_dSdM = zeros(SM)
dSdM = zeros(SM)
 
for i in range(shape(M)[0]) :
	j=0
	while j < (shape(k)[0]-1) :
		d_dSdM[i] = (k[j+1]-k[j])/(2.0*(pi)**2.0*M[i]**2.0)*PS[j]*2.0*W_M[i,j]/k[j]*Cc*(3.0*(M[i]*k[j]**3.0/Cc)**(1.0/3.0)*cos((M[i]*k[j]**3.0/Cc)**(1.0/3.0))+((M[i]*k[j]**3.0/Cc)**(2.0/3.0)-3.0)*sin((M[i]*k[j]**3.0/Cc)**(1.0/3.0)))
		dSdM[i] = dSdM[i] + d_dSdM[i]
		j=j+1

#mass function
dndM = zeros(SM, dtype = 'd')
dndM = - f_s_m/(2.0*s_m_sq)*rho_m/M*dSdM

toplot = zeros(SM, dtype = 'd')
toplot = M*dndM #M**2.0/rho_m*dndM

#get ratios
ratio_toplot=(dndlnM-toplot)/toplot
s_ratio_toplot=sdndlnM/toplot

################################################################################################################
#plot
f=plt.figure()
gs = gridspec.GridSpec(2, 1,height_ratios=[4,1])
b=plt.subplot(gs[0])
b.plot(Mtink,dndlnMtink,'k-',label=r'$\text{Tinker+ (2008)}$')
b.plot(M,dndlnM,'r.',markersize=12, label=r'\text{Separate Universe}')
b.errorbar(M,dndlnM, xerr=0.0, yerr=sdndlnM,fmt='r.',markersize=12)
b.set_ylabel(r'${\rm d}n/{\rm dln}M$',fontsize=26,labelpad=35)
b.set_xscale('log')
b.set_yscale('log')
#b.set_xlim([12.5,14.8])
#b.set_ylim([0.8,5.0])
b.legend(loc=0, prop={'size':16})

#ratios	
ratio=plt.subplot(gs[1], sharex=b)
ratio.plot([10**12,10**16],[0,0],'k-')
ratio.plot(M,ratio_toplot,'r.')
ratio.errorbar(M,ratio_toplot,xerr=0.0,yerr=s_ratio_toplot,fmt='r.')
ratio.set_xlabel(r'$M \left[\log h^{-1} M_\odot\right]$',fontsize=22, labelpad=15)
ratio.set_ylabel(r'$\frac{data-Tinker}{Tinker}$', fontsize=26)
ratio.set_xscale('log')
#ratio.set_xlim([12.5,14.8])
#ratio.set_ylim([-0.08,0.025])
#ratio.set_yticks(arange(-0.09,0.025,0.04))

f.subplots_adjust(hspace=0.05)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.show()
