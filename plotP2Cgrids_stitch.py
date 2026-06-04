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


def find_index(i):
	switcher={'density':1, 'velocity_x':2, 'velocity_y':3, 'velocity_z':4, 'temperature':5, 'stars':6, 'emissivity':7, 'SB':7, 'density_weighted_velz':8, 'density_weighted_vely':9, 'density_weighted_velz':10}
	
	return switcher.get(i,"Invalid quantity to plot\n Choices are density, velocity_x, velocity_y, velocity_z, temperature, stars, emissivity (or SB), density_weighted_velx, density_weighted_vely, density_weighted_velz")

def colormap(i):
	switcher={'density':'viridis', 'temperature':'plasma','emissivity':'inferno', 'SB':'jet'}

	return switcher.get(i,"Invalid quantity to plot")

# input parameters
parser = argparse.ArgumentParser(description='Routine to plot P2C output grids. Args: < path to hdf5 files (str) >, < suffix of outputfile (str) >, < what to plot (str) >, < plane to plot (str: XY, YZ, XZ)>, < width of slice (int, in cells)>')
parser.add_argument('path',nargs=1,type=str, help='path to P2C output (hdf5)')
parser.add_argument('suffix',nargs=1,type=str, help='Suffix of outputfile') 
parser.add_argument('toplot',nargs=1,type=str, help='Quantity to plot')
parser.add_argument('plane',nargs=1,type=str, help='Plane to plot')
parser.add_argument('width',nargs=1,type=int, help='Width of slice (grid cells)')
args   = parser.parse_args()
path   = args.path[0]
suffix = args.suffix[0]
toplot = args.toplot[0]
plane  = args.plane[0]
width  = args.width[0]
h5files  = [f for f in os.listdir(path+"/") if os.path.isfile(os.path.join(path+"/", f)) and f[-11:] == 'subbox.hdf5']
h5files  = sorted(h5files)

#some hard-coded settings
red       = 3
h         = 0.677
Lbox      = 4./h #Mpc
scale     = 'log' #'linear' or 'log'
cmap      = colormap(toplot) #'viridis', 'nipy_spectral','Spectral_r','Set3','Paired', 'jet'
if toplot == 'emissivity':
	min_q = 1e-2
	max_q = 1e7
elif toplot == 'SB':
	min_q = 1e-23
	max_q = 1e-15
	if red == 3:
	#	const = 1.04519029e-25	# const for 2.88 pkpc spatial resolution
		const = 5.2332108e-26   # const for 1.44 pkpc spatial resolution
	elif red == 3.017:
		 const = 5.49007030827e-26 # Eagle sim from stephanie	
	elif red == 3.5:
		const = 1.17767656e-25  # const for 2.88 pkpc spatial resolution
	elif red == 7:
		const = 2.09328431e-25  # const for 2.88 pkpc spatial resolution
elif toplot == 'density':
	min_q = 1e-1
	max_q = 7e2
elif toplot == 'temperature':
	min_q = 1e4
	max_q = 1e7
else:
	min_q = None
	max_q = None

# what to do
verbose  = False
img2d    = True
hist     = False
savedata = False #whether to save the data in a .npy file

print "2D image: "+str(img2d)
print "Histogram: "+str(hist)
print "Saving data: "+str(savedata)

if savedata == True:
	outform = 'npy' #'fits' or 'npy'
	print "Output data format: "+outform

print "Reading data..."

#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"
Num   = find_index(toplot)
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

#stich at "x=0"
data00 = dstack((datasub[0],datasub[1]))  
data01 = dstack((datasub[2],datasub[3]))
data0  = hstack((data00,data01))

del data00, data01

#stich at "x=1"
data10 = dstack((datasub[4],datasub[5]))
data11 = dstack((datasub[6],datasub[7]))
data1  = hstack((data10,data11))

del data10, data11

#finals stich
data   = vstack((data0,data1))
del datasub, data0, data1

sideL  = shape(data)[0]
gridtompc = Lbox/sideL

if width > sideL:
	sys.exit("Width is greater than 1D number of cells!")

if verbose:
	print shape(data)
	print("maximum "+toplot+": {}".format(amax(data)))
	print("minimum "+toplot+": {}".format(amin(data)))

print("Done")

#compute the projected quantity on the width in the good plane
# transpose by hand since the plotting routine imshow transposes
print("Preparing data to plot...")
dataplot=empty((sideL,sideL))
if verbose:
	print "Cells summed: "+str((int(sideL/2)-int(width/2)))+":"+str((int(sideL/2)+int(width/2)))
if plane == 'XY':
	for i in xrange(sideL):
		for j in xrange(sideL):
			dataplot[i,j]=sum(data[j,i,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2))])/width
elif plane == 'YZ':
	for i in xrange(sideL):
		for j in xrange(sideL):
			dataplot[i,j]=sum(data[(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),j,i])/width
elif plane == 'XZ':
	for i in xrange(sideL):
		for j in xrange(sideL):
			dataplot[i,j]=sum(data[j,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),i])/width
else:
	sys.exit("Selected plane not valid. Options are XY, YZ, or XZ") 

if toplot == 'SB':
	dataplot = dataplot*width*const #this factor depends on the redshift and spatial resolution

if savedata:
	print "Saving data..."
	if outform == 'npy':
		save(toplot+'_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+'_'+suffix,dataplot)
	else:
		hdu = fits.PrimaryHDU(dataplot)
		hdu.writeto(toplot+'_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+'_'+suffix+'.fits')

if verbose:
	print shape(dataplot)

print "Done"

############### 2D IMAGE ###############
if img2d:
	print 'Creating the 2d figure...'
	fig = figure(1)    #create the figure
	ax1 = fig.add_subplot(111)

	#ax1.set_xlim([0.0,Lbox])
	#ax1.set_ylim([0.0,Lbox])
	ax1.set_xlim([1.42,4.49])
	ax1.set_ylim([1.42,4.49])

	ax1.set_title(r'$z = $'+str(red)+', depth = '+str(round(width*gridtompc,2))+' cMpc',fontsize=20)
	
	if plane == 'XY':
		ax1.set_xlabel(r'X [cMpc]',fontsize=18)  #x-axis label
		ax1.set_ylabel(r'Y [cMpc]',fontsize=18)  #y-axis label
	elif plane == 'XZ':
		ax1.set_xlabel(r'X [cMpc]',fontsize=18)  #x-axis label
		ax1.set_ylabel(r'Z [cMpc]',fontsize=18)  #y-axis label
	else:
		ax1.set_xlabel(r'Y [cMpc]',fontsize=18)
		ax1.set_ylabel(r'Z [cMpc]',fontsize=18)	

	if min_q is None:  min_q = amin(dataplot)	
	if max_q is None:  max_q = amax(dataplot)

	if scale=='linear':
    		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)
	else:
    		cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))

	cbar = fig.colorbar(cax)
	if toplot == 'SB':
                cbar.set_label(toplot+r" $[\rm{erg} \rm{s}^{-1} \rm{cm}^{-2} \rm{arcsec}^{-2}]$",fontsize=20)
        elif toplot == 'temperature':
                cbar.set_label(toplot+" [K]",fontsize=20)
        else:
                cbar.set_label(toplot,fontsize=20)

	fout = toplot+'_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+'_'+suffix+'.pdf'
	savefig(fout, bbox_inches='tight')
	close(fig)

############### HIST ###############
if hist:
	print 'Creating the histogram...'
	lg    = True
	nbins = 100

	plt.figure(2)
	if lg:
		plt.hist(log10(dataplot.flatten()),bins=nbins,log=lg)
		plt.xlabel('log['+toplot+']',fontsize=18)
	else:
		plt.hist(dataplot.flatten(),bins=nbins,log=lg)
		plt.xlabel(toplot,fontsize=18)
	plt.show()
