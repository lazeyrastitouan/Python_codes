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
import pyfftw
from fftw import*


#FUNCTIONS
def find_index(i):
	switcher={'density':1, 'velocity_x':2, 'velocity_y':3, 'velocity_z':4, 'temperature':5, 'stars':6, 'emissivity':7, 'SB':7, 'density_weighted_velz':8, 'density_weighted_vely':9, 'density_weighted_velz':10}
	
	return switcher.get(i,"Invalid quantity to plot\n Choices are density, velocity_x, velocity_y, velocity_z, temperature, stars, emissivity (or SB), density_weighted_velx, density_weighted_vely, density_weighted_velz")


def colormap(i):
	switcher={'density':'viridis', 'temperature':'plasma','emissivity':'inferno', 'SB':'jet', 'velocity':'seismic'}

	return switcher.get(i,"Invalid quantity to plot")


#INPUT PARAMETERS
parser = argparse.ArgumentParser(description='Routine to compute curl of a P2C velocity grid. Args: < path to P2C file (str) >, < hdf5 file (str) >')
parser.add_argument('path',nargs=1,type=str, help='path to P2C output') 
parser.add_argument('hdf5file',nargs=1,type=str, help='P2C output (hdf5)')
args    = parser.parse_args()
path    = args.path[0]
p2cfile = args.hdf5file[0]

#some hard-coded settings
red  = 3.017
h    = 0.6777
Lbox = 50.0 #cMpc
cpus = 1

#WHAT TO DO
verbose  = False
img2d    = True
print("2d image: "+str(img2d))
outform = 'fits' #'fits' or 'npy' or 'hdf5 (chombo)
print("Output data format: "+outform)

#LOAD
print("Reading data...")

#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"
# open the P2C file(s)
f = h5py.File(path+'/'+p2cfile, "r")
#grid size in 1D
sideL     = array(f['level_0/boxes'])[0][3]+1
gridtompc = Lbox/sideL
# navigate to the relevant dataset (1D array to re-order)
dset = f[CPATH]
# NB: python does not include the last index in [i:j]
# NB: the reshaping has to use Fortran-like indexing as the chombo-file was originally produced using Fortran-like reshaping
Num = find_index('velocity_x')
vx  = dset[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
Num = find_index('velocity_y')
vy  = dset[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
Num = find_index('velocity_z')
vz  = dset[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
Num = find_index('density')
dta = dset[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
f.close()

if verbose:
	print(shape(vx))
print("Done")

#compute density contrast
mean = mean(dta)
dta  = (dta/mean)-1.

#FFTW
kF  = 2.*pi/Lbox
mid = sideL/2
kN  = 2.*pi/(2.*Lbox/sideL)
vxk = FFT3Dr_d(vx,cpus)
vyk = FFT3Dr_d(vy,cpus)
vzk = FFT3Dr_d(vz,cpus)

#CURL AND DIVERGENCE
curlxk = empty((sideL,sideL,sideL/2+1))
curlyk = empty((sideL,sideL,sideL/2+1))
curlzk = empty((sideL,sideL,sideL/2+1))
divk   = empty((sideL,sideL,sideL/2+1))
for kxx in range(sideL):
	kx = (kF*(kxx-sideL) if (kxx>mid) else kF*kxx)
	for kyy in range(sideL):
		ky = (kF*(kyy-sideL) if (kyy>mid) else kF*kyy)
		for kzz in range(sideL/2+1):
			kz = (kF*(kzz-sideL) if (kzz>mid) else kF*kzz)

			# kz=0 and kz=middle planes are special
			if kz==0 or (kz==mid and sideL%2==0):
				if kx<0: continue
				elif kx==0 or (kx==mid and sideL%2==0):
					if ky<0.0: continue

			curlxk[kxx,kyy,kzz] = I*(ky*vzk - kz*vyk) 
			curlyk[kxx,kyy,kzz] = I*(kz*vxk - kx*vzk)
			curlzk[kxx,kyy,kzz] = I*(kx*vyk - ky*vxk)
			divk[kxx,kyy,kzz]   = I*(kx*vxk + ky*vyk + kz*vzk)

#IFFTW
curlx = IFFTW3Dr_d(curlxk,cpus)
curly = IFFTW3Dr_d(curlyk,cpus)
curlz = IFFTW3Dr_d(curlzk,cpus)
div   = IFFTW3Dr_d(divk,cpus)

#SAVING
print("Saving data...")
if outform == 'npy':
        save('curl_'+p2cfile[:-5],curl)
elif outform == 'fits':
        hdu = fits.PrimaryHDU(curl)
        hdu.writeto('curl_'+p2cfile[:-5]+'.fits')
#else:
print("Done")

#PLOT
if img2d:
	# plane = 'XY' assumed here. Otherwise use plotP2Cgrid.py
	width = sideL	
	cmap1 = colormap('density') #'viridis', 'nipy_spectral','Spectral_r','Set3','Paired'  
	cmap2 = colormap('velocity')
	cmap3 = colormap('velocity')

	#compute the projected quantity on the width in the good plane
	#transpose by hand as imshow transposes later when plotting
	print("Preparing data to plot...")
	dtaplot  = empty((sideL,sideL))
	divplot  = empty((sideL,sideL))
	curlplot = empty((sideL,sideL))
	if verbose:
		print("Cells summed: "+str((int(sideL/2)-int(width/2)))+":"+str((int(sideL/2)+int(width/2))))
	for i in range(sideL):
		for j in range(sideL): #assumes plane = XY
			dtaplot[i,j]  = sum(dta[j,i,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2))])/width
			divplot[i,j]  = sum(div[j,i,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2))])/width
			curlplot[i,j] = sum(curlx[j,i,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2))])/width

	print('Creating the 2d figure...')
	fig = figure(1)    #create the figure
	ax1 = fig.add_subplot(131)
	ax2 = fig.add_subplot(132)
	ax3 = fig.add_subplot(133)

	ax1.set_xlim([0.0,Lbox])
	ax1.set_ylim([0.0,Lbox])
	ax2.set_xlim([0.0,Lbox])
	ax2.set_ylim([0.0,Lbox])
	ax3.set_xlim([0.0,Lbox])
	ax3.set_ylim([0.0,Lbox])	

	ax1.set_title(r'Density',fontsize=20)
	ax2.set_title(r'Divergence',fontsize=20)
	ax3.set_title(r'Curl x component',fontsize=20)
	
	ax1.set_xlabel(r'X [cMpc]',fontsize=18)
	ax1.set_ylabel(r'Y [cMpc]',fontsize=18)
	ax2.set_xlabel(r'X [cMpc]',fontsize=18)  
	ax2.set_ylabel(r'Y [cMpc]',fontsize=18)
	ax3.set_xlabel(r'X [cMpc]',fontsize=18)  
	ax3.set_ylabel(r'Y [cMpc]',fontsize=18)

    	#cax1 = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)
    	#log scale ssumed, change if necessary
	cax1  = ax1.imshow(dtaplot+2.,cmap=get_cmap(cmap2),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=1e-1,vmax=1e3))
	cbar1 = fig.colorbar(cax1)
	#cbar1.set_label('Density',fontsize=20)
	cax2  =	ax2.imshow(divplot,cmap=get_cmap(cmap2),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=amin(divplot),vmax=amax(divplot)))
	cbar2 = fig.colorbar(cax2)
	cax3  = ax3.imshow(curlplot,cmap=get_cmap(cmap3),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=amin(curlplot),vmax=amax(curlplot)))
	cbar3 = fig.colorbar(cax3)

	fout = 'dta_div_curlx_'+p2cfile[:-5]+'_XY_w'+str(width)+'.pdf'
	savefig(fout, bbox_inches='tight')
	close(fig)
