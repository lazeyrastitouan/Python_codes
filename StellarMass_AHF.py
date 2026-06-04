from numpy import*
import sys, os, h5py, argparse
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
#matplotlib.use('TkAgg')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

#Stellar mass to halo mass functions in litterature
#Tacchella 2018
def tacchella18(z,Mh):
    Mc    = 1.6e11
    eps_0 = 0.025
    eps_z = eps_0/(1.+z)**0.6
    b_0   = 1.0
    g_0   = 0.3

    return 2.*Mh*eps_z/((Mh/Mc)**(-b_0)+(Mh/Mc)**g_0)

#Moster 2009
def moster09(z,Mh):
    M1_0  = 10.**11.88
    mMo_0 = 0.0282
    b_0   = 1.06
    b_1   = 0.17
    g_0   = 0.556
    g_1   = -0.26
    nu    = -0.72
    mu    = 0.019

    M1z  = 10**(log10(M1_0)*(1.+z)**mu)
    mMoz = mMo_0*(z+1.)**nu
    gz   = g_0*(z+1.)**g_1
    bz   = b_1*z+b_0

    return 2.*Mhm*mMoz/((Mhm/M1z)**(-bz)+(Mhm/M1z)**gz)

#From Moster et al. 2013
def moster(zm,hmass=0):
    Mh=logspace(8,14,1000)
    if hmass.any()>0: Mh=hmass
    lm1  =11.590  + (1.195) *zm/(1+zm)
    N    = 0.0351 - (0.0247)*zm/(1+zm)
    beta = 1.376  - (0.826) *zm/(1+zm)
    gamma= 0.608  + (0.329) *zm/(1+zm)
    m1= 10**lm1
    shm = Mh * 2*N *((Mh/m1)**(-beta)+(Mh/m1)**gamma)**(-1)

    return array([Mh,shm])

#From Moster 2018
def moster18(zm,hmass=0):
    Mh=logspace(8,14,1000)
    #hmass=1e12
    #zm=0.1
    if hmass.any()>0: Mh=hmass

    #My fit
    lm1  = 11.739 + 0.410 * zm/(1+zm)
    N    =  0.160 + 0.010 * zm
    beta =  1.840 - 0.600 * zm/(1+zm)
    gamma=  0.59
        #lm1  = 12.10 # + (1.195) *zm/(1+zm)
        #N    = 0.24 # - (0.0247)*zm/(1+zm)
        #beta = 1.3  #- (0.826) *zm/(1+zm)
        #gamma= 0.60 # + (0.329) *zm/(1+zm)
        #print(lm1,N,beta,gamma)

    m1= 10**lm1
    shm = Mh * 2*N *((Mh/m1)**(-beta)+(Mh/m1)**gamma)**(-1) * 0.0484/0.308
    #print(shm/Mh/0.0484*0.308)
    Msigma = where(zm<1.0,10.85-zm,10.75-0.05*(zm-1))
    alpha  = 2.40775*zm**-0.69197 - 1.8125 #Power-law Fit
    sigma0 = where(zm<1,0.16-0.04*zm,where(zm<4,0.12-0.02*(zm-1),0.06-0.01*(zm-4))) #Piece-wise fit
    sigma=sigma0+log10((Mh/10**Msigma)**(-0.3)+1)
    return array([Mh,shm,10**(log10(shm)-sigma),10**(log10(shm)+sigma),10**(log10(shm)-2*sigma),10**(log10(shm)+2*sigma)])

#From Behroozi et al. 2013
def _fx(x,alpha,delta,gamma):
    return -log10(10**(alpha*x)+1)+delta*(log10(1+exp(x)))**gamma/(1+exp(10**(-x)))

def behroozi(redshift,hmass=0):
    """ 
    Computes Behroozi stellar to halo mass function at z=redshift.
    If hmass=0, then a log-spaced interval between 1e8 and 1e14 solar masses is used
    If hmass!=0, the interval can be defined by the user
    """
    #redshift=str(round(float(redshift),0))
    #redshift=min(redshift,8)
    #Computes the scaling parameters at z
    a=1/(1+redshift)
    nu_a=exp(-4*a**2)
    eps_log=-1.777+(-0.006*(a-1)+0.0*redshift)*nu_a-0.119*(a-1)
    M1_log =11.514+(-1.793*(a-1)+(-0.25*redshift))*nu_a
    alpha=-1.412+(0.731*(a-1))*nu_a
    delta=3.508+(2.608*(a-1)+(-0.043*redshift))*nu_a
    gamma=0.316+(1.319*(a-1)+0.279*redshift)*nu_a
    M1=10**M1_log
    eps=10**eps_log
    Mh=logspace(8,14,1000)
    if(hmass.any()>0): Mh=hmass
    Mstar_log = log10(eps*M1)+_fx(log10(Mh)-M1_log,alpha,delta,gamma)-_fx(0,alpha,delta,gamma)

    return array([Mh,10**Mstar_log])

# -------------------------------------------------------------------------------------------------#
#Start by loading: #[0]ID, [1]Host, [2]Nsub, [3]Mh, [4]fhires, [5]Mstar
root = "/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/AHF_007_allparts/"
z    = 3.0
title = 'z = '+str(z)+', CDM'
print('Loading...')
cat = loadtxt(root+"/200_snap_z"+str(z)+"00_halos_sorted",usecols=(0,1,2,3,37,64),dtype={'names':('ID','Host','Nsub','Mh','fh','Ms'), 'formats':(int64,int64,int32,float,float,float)})
h   = 0.6777

#subhalos
with open(root+"/200_snap_z"+str(z)+"00_substructure") as FILE:
	lines = list(filter(None,(line.rstrip() for line in FILE)))
nhalos = int32(len(lines)/2)
host, nsubs, = [], []
sub          = [[] for x in range(nhalos)]
for i in range(nhalos):
	host.append(int64(lines[2*i].split('        ')[0]))
	nsub = int32(lines[2*i].split('        ')[1])
	nsubs.append(nsub)
	for j in range(nsub):
		sub[i].append(int64(lines[2*i+1].split('\t')[j+1]))

print('Selecting and computing central properties...')
Mhp, Msp = [],[]
#select
for i in range(len(cat['ID'])):
	if cat['Ms'][i] != 0.0 and cat['Host'][i] == 0 and cat['fh'][i] > 0.95:
        #this whole loop can be replaced by e.g. filter=cat['Ms']!=0.0; cat[filter]; 
		Mhp.append(cat['Mh'][i])
		if cat['Nsub'][i] == 0:
			Msp.append(cat['Ms'][i])
		else:	
			#go find its subhalos in the previously formed subhalos catalog
			subid = host.index(cat['ID'][i])
			subs  = sub[subid]
			#remove their stellar mass
			for si in subs:
				subidincat = cat['ID'].tolist().index(si)
				cat['Ms'][i] = cat['Ms'][i]-cat['Ms'][subidincat]
			Msp.append(cat['Ms'][i])		


Mhp=array(Mhp)
Msp=array(Msp)

fb     = 0.0483/0.307
Mhm    = 10.**linspace(log10(1e7),log10(1e14),100)
Msm09  = moster09(z,Mhm)
Mst18  = tacchella18(z,Mhm)
Msm12  = moster(z,Mhm)[1]
Msb13  = behroozi(z,Mhm)[1]
Msm18  = moster18(z,Mhm)[1]

#plot
plt.plot([8.,14.],[8+log10(fb),14.+log10(fb)],'k-',label=r'Cosmic mean')
#plt.plot(log10(Mhm),log10(Msm09),'b-',label=r'Moster+ 2009')
plt.plot(log10(Mhm),log10(Msm12),'g-',label=r'Moster+ 2013')
plt.plot(log10(Mhm),log10(Msb13),'r-',label=r'Behroozi+ 2013')
plt.plot(log10(Mhm),log10(Msm18),'y-',label=r'Moster+ 2018')
if z == 4.0:
    plt.plot(log10(Mhm),log10(Mst18),'c-',label=r'Tacchella+ 2018')
plt.plot(log10(Mhp/h),log10(Msp/h),'m.',label=r'DaLya')
#plt.title(r''+title,fontsize=16) 
plt.text(11.1,11.4,r'z='+str(z),fontsize=18)
plt.legend(loc=0,prop={'size':14}, frameon=False)
plt.xlabel(r"$\log M_h \; [M_\odot]$",fontsize=20, labelpad=14)
plt.ylabel(r"$\log M_\star \; [M_\odot]$",fontsize=20,labelpad=15)
plt.xlim([11.,13.0])
plt.ylim([8.5,11.75])
plt.tight_layout()
plt.savefig(root+'/Mh_Mstar_centrals_z'+str(z)+'_zoom1p5_hydro_fire2.pdf', bbox_inches='tight')
print("Saved figure to "+root+'/Mh_Mstar_centrals_z'+str(z)+'_zoom1p5_hydro_fire2.pdf')
plt.show()
