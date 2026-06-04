from numpy import*
import h5py, sys, os
import argparse

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import LogNorm


# input parameters
parser = argparse.ArgumentParser(description='Routine to plot SB grids. Args: < path to SB file .npy (str) >, < optional -- path to second .npy file (str)>, < plane to plot (str) >, < SB threshold (float) > ')
parser.add_argument('path',nargs=1,type=str, help='Path to .npy file')
parser.add_argument('--path2',nargs=1,type=str, help='Path to second .npy file optional)')
parser.add_argument('plane',nargs=1,type=str,help='Plane to plot')
parser.add_argument('cut',nargs=1,type=float,help='Threshold of SB')
args   = parser.parse_args()
path   = args.path[0]
plane  = args.plane[0]
cut    = args.cut[0]
if args.path2 is not None:
	path2  = args.path2[0]

if plane != 'XY' and plane != 'YZ' and plane != 'XZ':
        sys.exit("Selected plane not valid. Options are XY, YZ, or XZ")

#some hard-coded settings
red       = 3
h         = 0.677
Lbox      = 4./h #Mpc
scale     = 'log' #'linear' or 'log'
cmap      = 'jet' #'nipy_spectral','Spectral_r','Set3','Paired'
width     = 1024

# what to do
verbose  = False
diff     = False
binary   = False
if diff == True and args.path2 is None:
	print("Please provide a second file for difference!")
	exit()

print "Reading data..."

data      = load(path)
if binary:
	dataplot = where(data>=cut,1.0,0.0)
	scale = 'linear'
else:
	dataplot = where(data>=cut,data,cut)
sideL     = shape(dataplot)[0]
gridtompc = Lbox/sideL

if verbose:
	print shape(dataplot)
	print("maximum SB: {}".format(amax(dataplot)))
	print("minimum SB: {}".format(amin(dataplot)))

#second file
if diff:
	print("Reading second file")
	dataplot2 = load(path2)
	if sideL != shape(dataplot2)[0]:
		sys.exit("The grids of the two files have different dimensions, can't compute difference!")

print("Done")

############### 2D IMAGE ###############
print 'Creating the 2d figure...'
fig = figure(1)    #create the figure
ax1 = fig.add_subplot(111)

ax1.set_xlim([0.0,Lbox])
ax1.set_ylim([0.0,Lbox])

ax1.set_title(r'WDM, $z = $'+str(red), fontsize=20)#+', depth = '+str(round(width*gridtompc,2))+' Mpc',fontsize=20)
	
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
	
	fout = path[:-4]+'_cut'+str(cut)+'_binary.pdf'

else:
	min_q = amin(cut)
	max_q = 1e-15 #amax(dataplot)

	if scale=='linear':
		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='sinc',vmin=min_q,vmax=max_q)
	else:
		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

	cbar = fig.colorbar(cax)
	cbar.set_label(r'SB [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]',fontsize=18)

	fout = path[:-4]+'_cut'+str(cut)+'.pdf'

savefig(fout, bbox_inches='tight')
close(fig)

############### DIFFERENCE ###############
if diff:
	print "Note ready yet"
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

