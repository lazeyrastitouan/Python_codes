from numpy import*
import h5py, sys, os
import argparse
from astropy.io import fits
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import LogNorm
from astropy.io import fits

# input parameters
parser = argparse.ArgumentParser(description='Routine to plot MockOb SB grids. Args: < MockOb file (str) >, < redshift (float) >, < dark matter type (str) >,< optional -- second MockOb file (str)>')
parser.add_argument('mockob',nargs=1,type=str, help='MockOb .fits file')
parser.add_argument('z',nargs=1,type=float,help='Redshift')
parser.add_argument('DM',nargs=1,type=str,help='Dark matter type')
parser.add_argument('--mockob2',nargs=1,type=str, help='Second MockOb .fits file (optional)')
args   = parser.parse_args()
f1     = args.mockob[0]
red    = args.z[0]
DM     = args.DM[0]
diff = False
if args.mockob2 is not None:
    diff = True
    f2  = args.mockob2[0]

#some hard-coded settings
specres   = 1.25 #Angstrom
h         = 0.6777 #0.681 0.6777 #0.671
Lbox      = 2.0 #2.0/h #cMpc
scale     = 'log' #'linear' or 'log'
cmap      = 'jet' #'nipy_spectral','Spectral_r','Set3','Paired'
obs       = True
savefits  = False
savefigg  = True

# what to do
verbose  = True

print("Reading data...")
hdul      = fits.open(f1)
data      = hdul[0].data #velocity first
min_q     = 1e-19 #abs(amin(dataplot))
max_q     = 1e-13 #amax(dataplot)
if obs:
    dataplot = data.sum(axis=0)*specres*1e-20 #MockOb outputs in MUSE units, observed (so dimmed) SB
    min_q    = min_q/(1.+red)**4.
    max_q    = max_q/(1.+red)**4.
else:
    #intrinsic SB
    dataplot = data.sum(axis=0)*specres*1e-20*(1.+red)**4.
sideL     = shape(dataplot)[1]
gridtompc = Lbox/sideL

if verbose:
	print(shape(dataplot))
	print("maximum SB: {}".format(amax(dataplot)))
	print("minimum SB: {}".format(amin(dataplot)))

if savefits:
    #output disperse readable fits
    d = fits.PrimaryHDU(dataplot)
    if obs:
        d.writeto(f1[:-5]+"_s5_n0_SB_observed.fits",overwrite=True)
    else:
        d.writeto(f1[:-5]+"_s5_n0_SB_intrinsic.fits",overwrite=True)

#second file
if diff:
    print("Reading second file")
    hdul2     = fits.open(f2)
    data2     = hdul2[0].data
    if obs:
        dataplot2 = data2.sum(axis=0)*specres*1e-20
    else:
        dataplot2 = data2.sum(axis=0)*specres*1e-20*(1.+red)**4.
    if sideL != shape(dataplot2)[0]:
        sys.exit("The grids of the two files have different dimensions, can't compute difference!")

print("Done")

############### 2D IMAGE ###############
print('Creating the 2d figure...')
fig = figure(1)    #create the figure
ax1 = fig.add_subplot(111)

ax1.set_xlim([0.0,Lbox])
ax1.set_ylim([0.0,Lbox])

ax1.set_title(r''+DM+', $z = $'+str(red), fontsize=20)#+', depth = '+str(round(width*gridtompc,2))+' Mpc',fontsize=20)
ax1.set_xlabel(r'[cMpc]',fontsize=18)  #x-axis label
ax1.set_ylabel(r'[cMpc]',fontsize=18)  #y-axis label

if scale=='linear':
	cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='sinc',vmin=min_q,vmax=max_q)
else:
	cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

cbar = fig.colorbar(cax)
if obs:
    cbar.set_label(r'SB [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$ (1+z)$^{-4}$]',fontsize=18)
else:
    cbar.set_label(r'SB [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]',fontsize=18)

#plt.show()
if obs:
    fout = f1[:-5]+'_s5_n0_SB_observed.pdf'
else:
    fout = f1[:-5]+'_s5_n0_SB_intrinsic.pdf'
if savefigg:
    savefig(fout, bbox_inches='tight')
plt.show()
close(fig)

############### DIFFERENCE ###############
if diff:
	print("Not ready yet")
	#d1md2 = (dataplot-dataplot2)/dataplot2
	#print 'Creating the 2d figure of difference...'
	
	#fig = figure(2)
	#ax2 = fig.add_subplot(111)

	#ax2.set_xlim([0.0,Lbox])
	#ax2.set_ylim([0.0,Lbox])

	#ax2.set_title(r'$z = $'+str(red)+', depth = '+str(round(width*gridtompc,2))+' Mpc',fontsize=20)

	#if plane == 'XY':
	#	ax2.set_xlabel(r'$X\, [Mpc]$',fontsize=18)
	#	ax2.set_ylabel(r'$Y\,  [Mpc]$',fontsize=18)
	#elif plane=='XZ':
	#	ax2.set_xlabel(r'$X\, [Mpc]$',fontsize=18)
	#	ax2.set_ylabel(r'$Z\, [Mpc]$',fontsize=18)
	#else:
	#	ax2.set_xlabel(r'$Y\,  [Mpc]$',fontsize=18)
	#	ax2.set_ylabel(r'$Z\, [Mpc]$',fontsize=18)

	#min_q = 1e-3 #amin(d1md2)
	#max_q = amax(d1md2)

	#if scale=='linear':
	#	cax = ax2.imshow(abs(d1md2),cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)
	#else:
	#	cax = ax2.imshow(abs(d1md2),cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

	#cbar = fig.colorbar(cax)
	#cbar.set_label('ratio '+ toplot,fontsize=20)

	#fout = toplot+'_diff_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+'_'+suffix+'.pdf'
	#savefig(fout, bbox_inches='tight')
	#close(fig)

