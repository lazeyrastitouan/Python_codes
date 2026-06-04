from numpy import*
import h5py, sys, os
import argparse
from astropy.io import fits
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import LogNorm
from disperse_utils import ppairs

# input parameters
parser = argparse.ArgumentParser(description='Routine to plot MockOb SB grids. Args: < MockOb file (str) >, < redshift (float) >, < dark matter type (str) >')
parser.add_argument('mockob',nargs=1,type=str, help='MockOb .fits file')
parser.add_argument('z',nargs=1,type=float,help='Redshift')
parser.add_argument('DM',nargs=1,type=str,help='Dark matter type')
args   = parser.parse_args()
f1     = args.mockob[0]
red    = args.z[0]
DM     = args.DM[0]

#some hard-coded settings
specres   = 1.25 #Angstrom
h         = 0.6777 #0.671
Lbox      = 4.0 #2.0/h #cMpc
scale     = 'log' #'linear' or 'log'
cmap      = 'jet' #'nipy_spectral','Spectral_r','Set3','Paired'

# what to do
verbose = True
obs     = False

#Let's plot the persistence vs field value
fname_ppairs = f1+'.ppairs.NDnet.a.NDnet'
trim, cut = ppairs.read_ppairs_ascii(fname_ppairs, sigma_x=5., sigma_y=5.)
print('trim : ',trim)
print('cut: ', cut)

"""
print("Reading data...")
hdul      = fits.open(f1)
data      = hdul[0].data #velocity first
if obs:
    dataplot  = data.sum(axis=0)*specres*1e-20 #MockOb outputs in MUSE units, observed (so dimmed) SB
else:#
    #intrinsic SB
    dataplot  = data.sum(axis=0)*specres*1e-20*(1.+red)**4.
min_q     = 1e-13 #abs(amin(dataplot))
max_q     = 1e-19 #amax(dataplot)
sideL     = shape(dataplot)[1]
gridtompc = Lbox/sideL

if verbose:
	print(shape(dataplot))
	print("maximum SB: {}".format(amax(dataplot)))
	print("minimum SB: {}".format(amin(dataplot)))

#if obs:
#    save(f1[:-5]+"_s1_n1em23_SB_observed.npy",dataplot)
#else:
#    save(f1[:-5]+"_s1_n1em23_SB_intrinsic.npy",dataplot)

print("Done")

############### 2D IMAGE ###############
print('Creating the 2d figure...')
fig = figure(1)    #create the figure
ax1 = fig.add_subplot(111)

ax1.set_xlim([0.0,Lbox])
ax1.set_ylim([0.0,Lbox])

ax1.set_title(r''+DM+', $z = $'+str(red), fontsize=20)#+', depth = '+str(round(width*gridtompc,2))+' Mpc',fontsize=20)
ax1.set_xlabel(r'$[cMpc]$',fontsize=18)  #x-axis label
ax1.set_ylabel(r'$[cMpc]$',fontsize=18)  #y-axis label

if scale=='linear':
	cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='sinc',vmin=min_q,vmax=max_q)
else:
	cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

cbar = fig.colorbar(cax)
if obs:
    cbar.set_label(r'${\textrm SB}/(1+z)^4$ [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]',fontsize=18)
else:
    cbar.set_label(r'SB [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]',fontsize=18)

#plt.show()
if obs:
    fout = f1[:-5]+'_s1_n1em23_SB_observed.pdf'
else:
    fout = f1[:-5]+'_s1_n1em23_SB_intrinsic.pdf'
savefig(fout, bbox_inches='tight')
close(fig)
"""
