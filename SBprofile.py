from numpy import*
import eagleTools as eT
import matplotlib.pyplot as plt
import matplotlib
import argparse
from astropy.io import fits

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# input parameters
parser  = argparse.ArgumentParser(description='Compute SB profiles. Args: < SB file (str) >, < --nbins (grid, int) >, < --boxsize (ckpc, float) >, < --center (grid units, 2 ints)>')
parser.add_argument('sbfile',nargs=1,type=str, help='SB file (fits format, str)')
parser.add_argument('--binnum', nargs=1,type=int, help='Optional, mumber of bins over the full interval, will be equally spaced in log. Default 12 (int)')
parser.add_argument('--boxsize',nargs=1,type=float, help='Optional, size of box in comoving kpc, default 4000.0 (float)')
parser.add_argument('--center',nargs=2,type=int, help='Coordinate of pixel in grid units around which to compute the profile (optional, default is grid center, 2 ints)')
args    = parser.parse_args()
sbfile  = args.sbfile[0]
binnum  = 12
if args.binnum is not None:
	binnum = args.binnum[0]
Lbox    = 2000.0
if args.boxsize is not None:	
	Lbox    = args.boxsize[0]

plot    = False

dhu     = fits.open(sbfile)
sbmap   = dhu[0].data
Ngrid   = shape(sbmap)
center  = [Ngrid[0]/2,Ngrid[1]/2]
if args.center is not None:
    center = args.center[0:2]

radDist = eT.calc2DRadialDistance(Ngrid[0], Ngrid[1], center[0], center[1])

#bins, sbprofile, p16, p84 = eT.radialProfile_log(sbmap,binnum,radDist) #r bins in log
bins, sbmean, sberr = eT.radialProfileMean_log(sbmap,binnum,radDist) #r bins in log	

bins     = bins*Lbox/Ngrid[0]
binmean  = ((bins + roll(bins, 1))/2.0)[1:]


fout    = sbfile[:-5]+'_meanSBprofile'
print('Writing output to '+fout+'.npz')
#savez(fout,bins=bins,binmean=binmean,sbprofile=sbprofile,p16=p16,p84=p84)#,sbmean=sbmean,sberr=sberr)
savez(fout,bins=bins,binmean=binmean,sbprofile=sbmean,sberr=sberr)

if plot:
    plt.loglog()
    #plt.plot(binmean, sbprofile,'b.', label=r'Median')
    #plt.fill_between(binmean,p16,p84,alpha=0.3)
    plt.plot(binmean,sbmean,'r--',label=r'Mean')
    plt.fill_between(binmean,sbmean-sberr,sbmean+sberr,alpha=0.3)
    plt.show()
