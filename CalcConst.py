from numpy import*

def calcconst(z, L, grid):
	"""
	Function to calculate SB constant to convert from emissivity based on redshft and resolution
	Arg; redshift (z), 1D boxsize (L) in comoving Mpc, Ngrid
	"""
	L   = 1./(1.+z)*L*1000.
	res = L/grid #phys kpc
	print("Resolution: "+str(res)+" pkpc")

	# hard-coded constants
	alpha     = 1.269e-13 #1.55e-13 [cm**3/s]
	rho_av    = 6.76991e-31 #average baryon density at z=0 [g/cm**-3]
	mH        = 1.673557e-24 #[g]
	d_corr    = (rho_av*(1.+z)**3.)/mH #because P2C outputs overdensity [cm**-3] 
	kpcTOcm   = 3.086e21 #[-]
	# MpcTOcm = 3.086e24
	ELya      = 10.2*1.60218e-12 #Lya energy in [erg]
	Y         = 0.24 #[-]

	const = d_corr**2. * (1.+z)**(-4.) * ((1.-0.5*Y)/(1.-Y)) * ELya/(4.*pi) * alpha * res * kpcTOcm/((2.0626e5)**2.)  
	print(const) 
	
	return const
