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
parser.add_argument('smooth',nargs=1,type=float, help='Sigma of gaussian smoothing kernel')
args   = parser.parse_args()
path   = args.path[0]
suffix = args.suffix[0]
plane  = args.plane[0]
width  = args.width[0]
smooth = args.smooth[0]

h5files  = [f for f in os.listdir(path+"/") if os.path.isfile(os.path.join(path+"/", f)) and f[-11:] == 'subbox.hdf5']
h5files  = sorted(h5files)

#some hard-coded settings
red       = 3
h         = 0.677
Lbox      = 4./h #Mpc
scale     = 'log' #'linear' or 'log'
cmap      = 'jet'
if red == 3:
	#const = 1.04519029e-25  #const for 2.88 pkpc spatial resolution
	const = 5.2332108e-26   # const for 1.44 pkpc spatial resolution
elif red == 3.5:
	const = 1.17767656e-25
elif red == 7:
	const = 2.09328431e-25

# what to do
verbose = True
plot    = True 
outform = 'fits' #npy or fits

print "Plot: "+str(plot)

print "Reading emissivity..."

#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"
Num   = 7 #emissivity
datasub = []

# open the P2C file(s)
for name in h5files:
	f = h5py.File(path+name, "r")

	# navigate to the relevant dataset (1D array to re-order)
	dset = f[CPATH]

	#subbox grid size in 1D
	sidesubL  = array(f['level_0/boxes'])[0][3]+1

	# check whether nvar and sLength from ReadParam make sense for the chombo-file being read
	if dset.shape[0] / nvar != sidesubL ** 3:
		print("chombofile.shape: {}".format(dset.shape[0]))
		print("sideLength: {}".format(sideL))
		print("nvar: {}".format(nvar))

		sys.exit("chombofile.shape/nvar != sideLength**3 --> wrong sideLength and/or wrong nvar have been given")

	# load quantity into the array "data"
	# NB: python does not include the last index in [i:j]
	# NB: the reshaping has to use Fortran-like indexing as the chombo-file was originally produced using Fortran-like reshaping
	datasub.append(dset[(Num - 1) * sidesubL ** 3 : (Num * sidesubL ** 3)].reshape((sidesubL, sidesubL, sidesubL), order="F"))
	f.close()

print("Done")

#now we stich, assuming 2x2x2 box, stich in z, y, x directions
# datasub [[000],[002],[020],[022],[200],[202],[220],[222]]
#x:vstack, y:hstack, z:dstack

print("Stitching boxes...")

#stitch at "x=0"
data00 = dstack((datasub[0],datasub[1]))
data01 = dstack((datasub[2],datasub[3]))
data0  = hstack((data00,data01))

del data00, data01

#stitch at "x=1"
data10 = dstack((datasub[4],datasub[5]))
data11 = dstack((datasub[6],datasub[7]))
data1  = hstack((data10,data11))

del data10, data11

#final stitch
data   = vstack((data0,data1))
del datasub, data0, data1

sideL  = shape(data)[0]
gridtompc = Lbox/sideL

if width > sideL:
	sys.exit("Width is greater than 1D number of cells!")

if verbose:
	print shape(data)
	print("maximum SB: {}".format(amax(data)))
	print("minimum SB: {}".format(amin(data)))

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
#w in pixel
#w = 2*int(truncate*sigma + 0.5) + 1 or (w-1)/2=int(truncate*sigma + 0.5)
# so for truncate = 4.0, smooth=0.1, or for truncate=3 smooth=0.15
truncate = 4.0
pix      = int(truncate*smooth + 0.5)
print("Smoothing over "+str(pix)+" pixels on each side")
dataplot = ndimage.gaussian_filter(dataplot, smooth, order=0, mode='wrap', cval=0.0, truncate=truncate)

#Add noise (STILL TO DO)

#Save
if outform == 'fits':
	print "Saving data to SB_z"+str(red)+"_Halo4_"+plane+"_w"+str(width)+"_s"+str(pix)+"p_"+suffix+".fits"
	hdu = fits.PrimaryHDU(dataplot)
	hdu.writeto('SB_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+'_s'+str(pix)+'p_'+suffix+'.fits')
else:
	print "Saving data to SB_z"+str(red)+"_Halo4_"+plane+"_w"+str(width)+"_s"+str(pix)+"p_"+suffix+".npy"
	save('SB_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+"_s"+str(pix)+"p_"+suffix,dataplot)

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

	min_q = 1e-23#amin(dataplot)	
	max_q = 1e-15#amax(dataplot)

	if scale=='linear':
    		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)
	else:
    		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

	cbar = fig.colorbar(cax)
	cbar.set_label(r'SB [erg cm$^{-2}$ s${-1}$ arcsec$^{-2}$]',fontsize=18)
	
	plt.show()
	close(fig)
