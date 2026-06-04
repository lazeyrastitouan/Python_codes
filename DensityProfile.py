from numpy import*
import eagleTools as eT
import matplotlib.pyplot as plt
import matplotlib
import argparse
from astropy.io import fits
import h5py
from decimal import Decimal

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# input parameters
parser  = argparse.ArgumentParser(description='Compute density profile. Args: < P2C file (str) >, < --ave_grid_mass [Msol] >, < --nbins (int) >, < --Tmin (K, float) >, < --Tmax (K, float) >, < --boxsize (ckpc, float) >, < --center (grid units, 3 ints)>')
parser.add_argument('p2cfile',nargs=1,type=str, help='P2C file (hdf5 format, str)')
parser.add_argument('--ave_grid_mass',nargs=1,type=float,help='P2C average mass in cell, default 2.93134E+03 (Msol, float)')
parser.add_argument('--binnum', nargs=1,type=int, help='Optional, mumber of bins over the full interval, will be equally spaced in log. Default 12 (int)')
parser.add_argument('--Tmin',nargs=1, type=float,help = 'Optional, minimum temperature threshold of gas (default 1K, float)')
parser.add_argument('--Tmax',nargs=1, type=float,help = 'Optional, maximum temperature threshold of gas (default 10^9K, float)')
parser.add_argument('--boxsize',nargs=1,type=float, help='Optional, size of box in comoving kpc, default 2000.0 (float)')
parser.add_argument('--center',nargs=3,type=int, help='Coordinate of pixel in grid units around which to compute the profile (optional, default is grid center, 3 ints)')
args    = parser.parse_args()
p2cfile = args.p2cfile[0]
binnum  = 12
if args.binnum is not None:
	binnum = args.binnum[0]
Lbox    = 2000.0
if args.boxsize is not None:	
	Lbox    = args.boxsize[0]
Tmax  = 1e9
Tmin  = 1.
if args.Tmin is not None:
	Tmin = args.Tmin[0]
if args.Tmax is not None:
    Tmax = args.Tmax[0]
ave_grid_mass = 2.93134E+03
if args.ave_grid_mass is not None:
    ave_grid_mass = args.ave_grid_mass[0]

plot    = False
z       = 3.0

# open the P2C file
#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"

f = h5py.File(p2cfile, "r")

#grid size in 1D
Ngrid     = array(f['level_0/boxes'])[0][3]+1
gridtokpc = Lbox/Ngrid
Msuntog   = 1.989e33
kpctocm   = 3.086e21

# navigate to the relevant dataset (1D array to re-order)
dset = f[CPATH]
# load density and temperature
# NB: python does not include the last index in [i:j]
# NB: the reshaping has to use Fortran-like indexing as the chombo-file was originally produced using Fortran-like reshaping
density = dset[0 : Ngrid**3].reshape((Ngrid, Ngrid, Ngrid), order="F") #this is the overdensity rho [cgs]/rho_b(z=0) [cgs]
temp    = dset[4 * Ngrid**3 : 5 * Ngrid**3].reshape((Ngrid, Ngrid, Ngrid), order="F")

center  = [Ngrid/2,Ngrid/2,Ngrid/2]
if args.center is not None:
    center = args.center[0:3]

radDist = eT.calc3DRadialDistance(Ngrid, Ngrid, Ngrid, center[0], center[1], center[2])

#select cells with temperature with Tmin and Tmax
dens    = where((temp>=Tmin) & (temp<=Tmax),density,nan)

#translate overdensity to hydrogen density in cgs (assuming full hydrogen since it is what is done in P2C)
mH      = 1.6735575e-24
denscgs = dens*ave_grid_mass/(gridtokpc*kpctocm)**3.*Msuntog*(1.+z)**3./mH #atom/pcm^3
#denssq     = dens**2.

#Radial profile
bins, dprofile, p16, p84 = eT.radialProfile_log(dens,binnum,radDist,0)
binscgs, dcgsprofile, p16cgs, p84cgs = eT.radialProfile_log(denscgs,binnum,radDist,0)
#bins, sbmean, sberr = eT.radialProfileMean_log(sbmap,binnum,radDist,thresh)
#binssq, dsqprofile, p16sq, p84sq = eT.radialProfile_log(denssq,binnum,radDist,0)

bins    = bins*gridtokpc
binmean = ((bins + roll(bins, 1))/2.0)[1:]

#dsqrtprofile = sqrt(dsqprofile)
#p16sqrt      = sqrt(p16sq)
#p84sqrt      = sqrt(p84sq)

#print(dsqrtprofile/dprofile)

fout    = 'DensityProfile_Tmin_'+str('{:.0e}'.format(Decimal(Tmin)))+'_Tmax_'+str('{:.0e}'.format(Decimal(Tmax)))+'_'+p2cfile[:-5]
print('Writing output to '+fout+'.npz')
savez(fout,bins=bins,binmean=binmean,overdensity=dprofile,p16=p16,p84=p84,density=dcgsprofile,p16dens=p16cgs,p84dens=p84cgs)#,sbmean=sbmean,sberr=sberr)

if plot:
    plt.loglog()
    plt.plot(binmean, dprofile,'b.', label=r'Median')
    plt.fill_between(binmean,p16,p84,alpha=0.3,color='blue')
    #plt.plot(binmean, dsqrtprofile,'g.', label=r'$\sqrt{<n^2>}$')
    #plt.fill_between(binmean,p16sqrt,p84sqrt,alpha=0.3,color='green')
    #plt.plot(binmean,sbmean,'r--',label=r'Mean')
    #plt.fill_between(binmean,sbmean-sberr,sbmean+sberr,alpha=0.3)
    plt.legend(loc=0)
    plt.show()
