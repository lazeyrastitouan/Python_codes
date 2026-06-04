from numpy import*
import argparse
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
from scipy import ndimage
from pylab import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

### Get mean SB and std dev in squares of nxn pixels in SB maps and output new grids. Output in disperse or HistSB.py readable .fits files. MockOb output must have been transformed to fits by MockOb2fits.py or plotMockOb.py) at various resolution ### 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

#input parameters
parser = argparse.ArgumentParser(description='Get mean SB and std dev in squares of nxn pixels in SB maps and output new grids. Output in disperse or HistSB.py readable .fits files. MockOb output must have been transformed to fits by MockOb2fits.py or plotMockOb.py) at various resolution. Args: < Collapsed MockOb file (str) > < 1D side length of regions to compute mean and variance in pixels (int) >')
parser.add_argument('mockob',nargs=1,type=str,help='Collpased MockOb (SB) .fits file as outputted by MockOb2fits.py or plotMockOb.py (absolute path)')
parser.add_argument('nxn',nargs=1,type=int,help='1D side length of regions to compute mean and variance in pixels')
args   = parser.parse_args()
nxn    = args.nxn[0]
SBfile = args.mockob[0]
foutm  = SBfile[:-5]+'_mean_'+str(nxn)+'x'+str(nxn)+'.fits'
foutv  = SBfile[:-5]+'_stddev_'+str(nxn)+'x'+str(nxn)+'.fits' #we will compute std dev so that it scales with (1+z)^4 as SB and mean SB
plot   = True

hdu    = fits.open(SBfile)
SB     = hdu[0].data

#mean and variance in nxn pixels squares
print("Computing mean...")
meanSB = ndimage.generic_filter(SB, nanmean, size=nxn, mode='constant', cval=NaN)
print("Computing variance...")
stdSB  = ndimage.generic_filter(SB, nanstd, size=nxn, mode='constant', cval=NaN)

#write
datamean = fits.PrimaryHDU(meanSB)
datastd  = fits.PrimaryHDU(stdSB)
datamean.writeto(foutm,overwrite=True)
datastd.writeto(foutv,overwrite=True)

#verification plot
if plot:
    print('Creating the 2d figure...')
    Lbox  = 4.0
    max_q = 5e-15 
    min_q = 5e-21
    
    #create the figure
    fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 1.085]})
    #fig.tight_layout(pad=5.0)
    ax1 = ax[0] #fig.add_subplot(121)
    ax2 = ax[1] #fig.add_subplot(122)

    ax1.set_xlim([0.0,Lbox])
    ax1.set_ylim([0.0,Lbox])
    ax2.set_xlim([0.0,Lbox])
    ax2.set_ylim([0.0,Lbox])

    ax1.set_title(r'Mean SB, '+str(nxn)+'x'+str(nxn)+' pixels', fontsize=16)
    ax1.set_xlabel(r'[cMpc]',fontsize=16)  #x-axis label
    ax1.set_ylabel(r'[cMpc]',fontsize=16)  #y-axis label
    ax2.set_title(r'Std dev SB, '+str(nxn)+'x'+str(nxn)+' pixels', fontsize=16)
    ax2.set_xlabel(r'[cMpc]',fontsize=16)  #x-axis label
    #ax2.set_ylabel(r'$[cMpc]$',fontsize=18)  #y-axis label

    cax1 = ax1.imshow(meanSB,cmap=get_cmap('jet'),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))
    cax2 = ax2.imshow(stdSB,cmap=get_cmap('jet'),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

    #cbar1 = fig.colorbar(cax1)
    cbar = fig.colorbar(cax2,fraction=0.046, pad=0.04)
    cbar.set_label(r'[erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]',fontsize=14)

    figout = SBfile[:-5]+'_mean_stddev_'+str(nxn)+'x'+str(nxn)+'.pdf'
    savefig(figout, bbox_inches='tight')

    #plt.show()
    close(fig)
