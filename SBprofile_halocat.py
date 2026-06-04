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
parser  = argparse.ArgumentParser(description='Compute SB profiles. Args: < SB file (str) >, < log10(Mlow) (float) >, < halo catalog (str) >, < los (str) >, < --nbins (grid, int) >, < --boxsize (ckpc, float) >, < --threshold (float)> ') #< center (grid units, 2 ints) >
parser.add_argument('sbfile',nargs=1,type=str, help='SB file (fits format, str)')
#parser.add_argument('center',nargs=2,type=int, help='Coordinate of pixel in grid units around which to compute the profile (2 ints)')
parser.add_argument('mlow',nargs=1,type=float,help='Mass bin lower edge in log')
parser.add_argument('halocat',nargs=1,type=str,help='Catalog of halos with center position in pixel (output of getPixelPosHalos.py), (str)')
parser.add_argument('los',nargs=1,type=str,help='Line of sight of the MockOb (X/Y/Z) (str)')
parser.add_argument('--binnum', nargs=1,type=int, help='Optional, mumber of bins over the full interval, will be equally spaced in log. Default 12 (int)')
parser.add_argument('--boxsize',nargs=1,type=float, help='Optional, size of box in comoving kpc, default 4000.0 (float)')
parser.add_argument('--threshold',nargs=1,type=float, help='Optional, SB threshold, default 0 (float)')
args    = parser.parse_args()
sbfile  = args.sbfile[0]
mlow    = args.mlow[0]
mhi     = mlow + 0.5
halocat = args.halocat[0]
los     = args.los[0] 
#center  = args.center[0:2]
binnum  = 12
if args.binnum is not None:
	binnum = args.binnum[0]
Lbox    = 4000.0
if args.boxsize is not None:	
	Lbox    = args.boxsize[0]
thresh  = 0
if args.threshold is not None:
	thresh = args.threshold[0]

plot    = False

dhu     = fits.open(sbfile)
sbmap   = dhu[0].data
Ngrid   = shape(sbmap)

if los == 'X':
	centers = loadtxt(halocat,usecols=(2,3))
elif los == 'Y':
	centers = loadtxt(halocat,usecols=(1,3))
elif los == 'Z':
	centers = loadtxt(halocat,usecols=(1,2))
else:
	exit('Invalid LOS!')

sbprofiles, p84s, p16s = [],[],[]
for center in centers:
	radDist = eT.calc2DRadialDistance(Ngrid[0], Ngrid[1], center[0], center[1])

	bins, sbprofile, p16, p84 = eT.radialProfile_log(sbmap,binnum,radDist,thresh)
	print(len(bins))
	print(bins)
	sbprofiles.append(sbprofile)
	p16s.append(p16)
	p84s.append(p84)
	#bins, sbmean, sberr = eT.radialProfileMean_log(sbmap,binnum,radDist,thresh)	

bins     = bins*Lbox/Ngrid[0]
binmean  = ((bins + roll(bins, 1))/2.0)[1:]
sbmedian = nanmedian(sbprofiles,axis=0)
sb16     = nanmedian(p16s,axis=0)
sb84     = nanmedian(p84s,axis=0)


#fout    = sbfile[:-5]+'_SBprofile_t'+str(thresh)+'_Mmin_'+str(mlow)+'_Mmax_'+str(mhi)
#print('Writing output to '+fout+'.npz')
#savez(fout,bins=bins,binmean=binmean,sbprofile=sbmedian,p16=sb16,p84=sb84)#,sbmean=sbmean,sberr=sberr)

if plot:
    plt.loglog()
    plt.plot(binmean, sbprofile,'b-', label=r'Median')
    plt.fill_between(binmean,p16,p84,alpha=0.3)
    plt.plot(binmean,sbmean,'r--',label=r'Mean')
    plt.fill_between(binmean,sbmean-sberr,sbmean+sberr,alpha=0.3)
    plt.show()
