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


# input parameters
parser = argparse.ArgumentParser(description='Routine to plot CubEx grids. Args: < path to fits file (str) >, < plane to plot (str) >')
parser.add_argument('path',nargs=1,type=str, help='Path to fits file')
parser.add_argument('plane',nargs=1,type=str,help='Plane to plot')
args   = parser.parse_args()
path   = args.path[0]
plane  = args.plane[0]

if plane != 'XY' and plane != 'YZ' and plane != 'XZ':
        sys.exit("Selected plane not valid. Options are XY, YZ, or XZ")

#some hard-coded settings
red       = 3
h         = 0.6777
Lbox      = 2. #Mpc
width     = 256
binary    = True
if binary:
	cmap = 'cividis'
else:
	cmap = 'nipy_spectral' #'cividis' #'nipy_spectral','Spectral_r','Set3','Paired'
print("Reading data...")

dhu       = fits.open(path)
data      = dhu[0].data

if binary:
	dataplot  = where(data[0,:,:]>0,1.0,0.0)
else:
	dataplot  = data[0,:,:]
sideL     = shape(dataplot)[0]
gridtompc = Lbox/sideL

print("Done")

############### 2D IMAGE ###############
print('Creating the 2d figure...')
fig = figure(1)    #create the figure
ax1 = fig.add_subplot(111)

ax1.set_xlim([0.0,Lbox])
ax1.set_ylim([0.0,Lbox])

ax1.set_title(r'$z = $'+str(red),fontsize=14)
	
if plane == 'XY':
	ax1.set_xlabel(r'$X\, [Mpc]$',fontsize=18)  #x-axis label
	ax1.set_ylabel(r'$Y\,  [Mpc]$',fontsize=18)  #y-axis label
elif plane == 'XZ':
	ax1.set_xlabel(r'$X\, [Mpc]$',fontsize=18)  #x-axis label
	ax1.set_ylabel(r'$Z\,  [Mpc]$',fontsize=18)  #y-axis label
else:
	ax1.set_xlabel(r'$Y\,  [Mpc]$',fontsize=18)
	ax1.set_ylabel(r'$Z\,  [Mpc]$',fontsize=18)	

if binary:
	cax  = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none')
else:
	cax  = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=1,vmax=2)

if not binary:
	cbar = fig.colorbar(cax)
	cbar.set_label('Object number',fontsize=20)
if binary:
	fout = path[:-5]+'_binary.pdf'
else:
	fout = path[:-5]+'.pdf'
print("Saving figure to "+fout)
savefig(fout, bbox_inches='tight')
close(fig)
