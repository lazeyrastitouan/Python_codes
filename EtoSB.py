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
from scipy import ndimage

# input parameters
parser = argparse.ArgumentParser(description='Routine to go from P2C emissivity to mock SB.')
parser.add_argument('path',nargs=1,type=str, help='Path to P2C output (hdf5)')
parser.add_argument('suffix',nargs=1,type=str, help='Suffix of outputfile') 
parser.add_argument('plane',nargs=1,type=str, help='Plane to plot')
parser.add_argument('width',nargs=1,type=int, help='Width of slice (grid cells)')
parser.add_argument('smooth',nargs=1,type=int, help='Width of smoothing kernel (grid cells)')
args   = parser.parse_args()
path   = args.path[0]
suffix = args.suffix[0]
plane  = args.plane[0]
width  = args.width[0]
smooth = args.smooth[0]

#some hard-coded settings
red       = 3
h         = 0.677
Lbox      = 4./h #Mpc
scale     = 'log' #'linear' or 'log'
cmap      = 'jet' #'viridis', 'nipy_spectral','Spectral_r','Set3','Paired'
if red == 3:
	const = 1.04519029e-25
elif red == 3.5:
	const = 1.17767656e-25
elif red == 7:
	const = 2.09328431e-25

# what to do
verbose = False
plot    = True #whether to save the data in a .npy file


print "Plot: "+str(plot)

print "Reading emissivity..."

#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"
Num   = 7 #emissivity

# open the P2C file(s)
f = h5py.File(path, "r")

#grid size in 1D
sideL     = array(f['level_0/boxes'])[0][3]+1
gridtompc = Lbox/sideL
if width > sideL:
	sys.exit("Width is greater than 1D number of cells!")

# navigate to the relevant dataset (1D array to re-order)
dset = f[CPATH]

# check whether nvar and sLength from ReadParam make sense for the chombo-file being read
if dset.shape[0] / nvar != sideL ** 3:
	print("chombofile.shape: {}".format(dset.shape[0]))
	print("sideLength: {}".format(sideL))
	print("nvar: {}".format(nvar))

	sys.exit("chombofile.shape/nvar != sideLength**3 --> wrong sideLength and/or wrong nvar have been given")

# load the emissivity and emissivity weighted velocities into the array "data"
# NB: python does not include the last index in [i:j]
# NB: the reshaping has to use Fortran-like indexing as the chombo-file was originally produced using Fortran-like reshaping
data  = dset[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F") #still in fortran order yxz

f.close()

if verbose:
	print shape(data)
	print("maximum"+toplot+": {}".format(amax(data)))
	print("minimum"+toplot+": {}".format(amin(data)))

print("Done")

#compute the projected quantity on the width in the good plane
print("Converting to SB and smoothing...")
dataplot=empty((sideL,sideL))
if verbose:
	print "Cells summed: "+str((int(sideL/2)-int(width/2)))+":"+str((int(sideL/2)+int(width/2)))
if plane == 'XY':
	for i in xrange(sideL):
		for j in xrange(sideL):
			dataplot[i,j]=sum(data[j,i,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2))])
elif plane == 'YZ':
	for i in xrange(sideL):
		for j in xrange(sideL):
			dataplot[i,j]=sum(data[(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),j,i])
elif plane == 'XZ':
	for i in xrange(sideL):
		for j in xrange(sideL):
			dataplot[i,j]=sum(data[j,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),i])
else:
	sys.exit("Selected plane not valid. Options are XY, YZ, or XZ") 

#Convert to SB
dataplot = dataplot*const #this factor depends on the redshift and spatial resolution

#Smooth to mimick PSF
dataplot = ndimage.gaussian_filter(dataplot, smooth, order=0, mode='wrap', cval=0.0, truncate=4.0)

#Add noise (STILL TO DO)

#Save
print "Saving data to SB_z"+str(red)+"_Halo4_"+plane+"_w"+str(width)+"_s"+str(smooth)+"p_"+suffix+".fits"
hdu = fits.PrimaryHDU(dataplot)
hdu.writeto('SB_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+'_s'+str(smooth)+'p_'+suffix+'.fits')

if verbose:
	print shape(dataplot)

print "Done"

############### 2D IMAGE ###############
if plot:
	print 'Creating the 2d figure...'
	fig = figure(1)    #create the figure
	ax1 = fig.add_subplot(111)

	ax1.set_xlim([0.0,Lbox])
	ax1.set_ylim([0.0,Lbox])

	ax1.set_title(r'$z = $'+str(red)+', depth = '+str(round(width*gridtompc,2))+' Mpc',fontsize=20)
	
	if plane == 'XY':
		ax1.set_xlabel(r'$X\, [Mpc]$',fontsize=18)  #x-axis label
		ax1.set_ylabel(r'$Y\,  [Mpc]$',fontsize=18)  #y-axis label
	elif plane == 'XZ':
		ax1.set_xlabel(r'$X\, [Mpc]$',fontsize=18)  #x-axis label
		ax1.set_ylabel(r'$Z\,  [Mpc]$',fontsize=18)  #y-axis label
	else:
		ax1.set_xlabel(r'$Y\,  [Mpc]$',fontsize=18)
		ax1.set_ylabel(r'$Z\,  [Mpc]$',fontsize=18)	

	min_q = 1e-23 #amin(dataplot)	
	max_q = 1e-15 #amax(dataplot)

	if scale=='linear':
    		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)
	else:
    		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

	cbar = fig.colorbar(cax)
	cbar.set_label(r'SB [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$]',fontsize=18)
	
	plt.show()
	close(fig)
