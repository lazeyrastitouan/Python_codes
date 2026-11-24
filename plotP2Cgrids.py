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
from CalcConst import*

def find_index(i):
	switcher={'density':1, 'velocity_x':2, 'velocity_y':3, 'velocity_z':4, 'temperature':5, 'stars':6, 'emissivity':7, 'SB':7, 'density_weighted_velz':8, 'density_weighted_vely':9, 'density_weighted_velz':10}
	
	return switcher.get(i,"Invalid quantity to plot\n Choices are density, velocity_x, velocity_y, velocity_z, temperature, stars, emissivity (or SB), density_weighted_velx, density_weighted_vely, density_weighted_velz")

def colormap(i):
	switcher={'density':'viridis', 'temperature':'plasma','emissivity':'inferno', 'SB':'jet'}

	return switcher.get(i,"Invalid quantity to plot")

tmin   = None
tmax   = None
suffix = ''

# input parameters
parser = argparse.ArgumentParser(description='Routine to plot P2C output grids. Args: < hdf5 file (str) >, < optional -- second hdf5 file (str)>, < suffix of outputfile (str) >, < what to plot (str) >, < plane to plot (str: XY, YZ, XZ)>, < width of slice (int, in cells)>, <optional -- max T (float, K)>, <optional -- min T (float, K)>')
parser.add_argument('hdf5file',nargs=1,type=str, help='P2C output (hdf5)')
parser.add_argument('--hdf5file2',nargs=1,type=str, help='P2C output for difference (hdf5, optional)')
parser.add_argument('--suffix',nargs=1,type=str, help='Suffix of outputfile (str, optional)')
parser.add_argument('toplot',nargs=1,type=str, help='Quantity to plot')
parser.add_argument('plane',nargs=1,type=str, help='Plane to plot')
parser.add_argument('width',nargs=1,type=int, help='Width of slice (grid cells)')
parser.add_argument('--tempmin',nargs=1,type=float, help='Minimum temperature (float, optional)')
parser.add_argument('--tempmax',nargs=1,type=float, help='Maximum temperature (float, optional)')
args   = parser.parse_args()
path   = args.hdf5file[0]
if args.suffix is not None:
	suffix = '_'
	suffix += args.suffix[0]
toplot = args.toplot[0]
plane  = args.plane[0]
width  = args.width[0]
if args.tempmin is not None:
	tmin   = args.tempmin[0]
if args.tempmax is not None:
	tmax   = args.tempmax[0]
diff = False
if args.hdf5file2 is not None:
	diff = True
	path2  = args.hdf5file2[0]

#SB with temperature condition not implemented correctly
if toplot == 'SB':
	if tmin is not None or tmax is not None:
		print("Sorry, temperature condition not implemented for SB :(")
		exit()

#some hard-coded settings
gn        = 3
red       = 3.5 #3.596 
a         = 1./(1.+red)
Lbox      = 4.0 #cMpc
scale     = 'log' #'linear' or 'log'
cmap      = colormap(toplot) #'viridis', 'nipy_spectral','Spectral_r','Set3','Paired'

# what to do
verbose  = True
img2d    = True
plotRH   = False
hist     = False
savedata = False #whether to save the data in a .npy file

print("2D image: "+str(img2d))
print("Plotting halo radius: "+str(plotRH))
print("Histogram: "+str(hist))
print("2D image difference: "+str(diff))
print("Saving data: "+str(savedata))

if savedata == True:
	outform = 'npy' #'fits' or 'npy'
	print("Output data format: "+outform)

print("Reading data...")

#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"
Num   = find_index(toplot)

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
#if dset.shape[0] / nvar != sideL ** 3:
#	print("chombofile.shape: {}".format(dset.shape[0]))
#	print("sideLength: {}".format(sideL))
#	print("nvar: {}".format(nvar))

#	sys.exit("chombofile.shape/nvar != sideLength**3 --> wrong sideLength and/or wrong nvar have been given")

# load the emissivity and emissivity weighted velocities into the array "data"
# NB: python does not include the last index in [i:j]
# NB: the reshaping has to use Fortran-like indexing as the chombo-file was originally produced using Fortran-like reshaping
data  = dset[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")


#if temperature boundaries imposed, load temperature and put to 0 all values not in interval
if tmin is not None or tmax is not None:
	temp  = dset[4 * sideL ** 3 : (5 * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
if tmin is not None:
	data  = where(temp>=tmin,0.0,data)
if tmax is not None:
	data  = where(temp<=tmax,0.0,data)
	
f.close()

if verbose:
	print(shape(data))
	print("maximum "+toplot+": {}".format(amax(data)))
	print("minimum "+toplot+": {}".format(amin(data)))


#second file
if diff:
	print("Reading second file")
	f2 = h5py.File(path2, "r")
	if sideL != array(f2['level_0/boxes'])[0][3]+1:
		sys.exit("The grids of the two files have different dimensions, can't compute difference!")
	dset2 = f2[CPATH]
	data2 = dset2[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F") #still in fortran order yxz
	
	#respect temperature condition
	if tmin is not None or tmax is not None:
		temp2 = dset2[4 * sideL **3 : 5 * sideL **3].reshape((sideL, sideL, sideL), order="F")
	if tmin is not None:
		data2 = where(temp2<tmin,0.0,data2)
	if tmax is not None:
		data2 = where(temp2>tmax,0.0,data2)	

	f2.close()
	
	dataplot2 = empty((sideL,sideL))	

print("Done")

#compute the projected quantity on the width in the good plane
#transpose by hand as imshow transposes later when plotting
print("Preparing data to plot...")
dataplot = empty((sideL,sideL))
Ncell    = width
if verbose:
	print("Cells summed: "+str((int(sideL/2)-int(width/2)))+":"+str((int(sideL/2)+int(width/2))))
if plane == 'XY':
    for i in range(sideL):
        for j in range(sideL):
			#we need to remove the 0 cells from the number of cells for the mean if temperature condition
            if tmin is not None or tmax is not None:
                Ncell = len(data[j,i][data[j,i]!=0.0])
                if Ncell == 0:
                    print("Warning, no cell satisfy temperature condition on location ("+str(j)+","+str(i)+")!")
                    Ncell == 1
            dataplot[i,j]=sum(data[j,i,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2))])/Ncell
    if diff:
        for i in range(sideL):
            for j in range(sideL):
                if tmin is not None or tmax is not None:
                    Ncell = len(data[j,i][data[j,i]!=0.0])
                dataplot2[i,j]=sum(data2[j,i,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2))])/Ncell
elif plane == 'YZ':
	for i in range(sideL):
		for j in range(sideL):
			#we need to remove the 0 cells from the number of cells for the mean if temperature condition
			if tmin is not None or tmax is not None:
				 Ncell = len(data[j,i][data[j,i]!=0.0])  
			dataplot[i,j]=sum(data[(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),j,i])/Ncell
	if diff:
		for i in range(sideL):
			for j in range(sideL):
				if tmin is not None or tmax is not None:
					Ncell = len(data[j,i][data[j,i]!=0.0])
				dataplot2[i,j]=sum(data2[(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),j,i])/Ncell
elif plane == 'XZ':
	for i in range(sideL):
		for j in range(sideL):
			#we need to remove the 0 cells from the number of cells for the mean if temperature condition
			if tmin is not None or tmax is not None:
				Ncell = len(data[j,i][data[j,i]!=0.0])
			dataplot[i,j]=sum(data[j,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),i])/Ncell
	if diff:
		for i in range(sideL):
			for j in range(sideL):
				if tmin is not None or tmax is not None:
					Ncell = len(data[j,i][data[j,i]!=0.0])
				dataplot2[i,j]=sum(data2[j,(int(sideL/2)-int(width/2)):(int(sideL/2)+int(width/2)),i])/Ncell
else:
	sys.exit("Selected plane not valid. Options are XY, YZ, or XZ") 

if toplot == 'emissivity':
    min_q = 1e-2
    max_q = 1e7
elif toplot == 'SB':
    const = calcconst(red, Lbox, sideL)
    dataplot = dataplot*width*const*(1.+red)**(-4.) #this factor depends on the redshift and spatial resolution
    if diff:
        dataplot2 = dataplot2*width*const*(1.+red)**(-4.)
    min_q = 1e-26
    max_q = 1e-18
elif toplot == 'density':
    min_q = 1e-4
    max_q =2e-3
    mH  = 1.6735575e-24 #gram
    mHe = 6.6464731e-24 #gram
    mA  = 0.76*mH+0.24*mHe
    Msuningram    = 1.989e33
    Mpcincm       = 3.086e24
    ave_grid_mass = 2.93139e3 #Msun
    #density is in fact M/ave_grid_mass so we convert to atom/cc
    dataplot = dataplot * ave_grid_mass / (a*gridtompc*Mpcincm)**3. * Msuningram/mA
elif toplot == 'temperature':
        min_q = 5e4
        max_q = 1e7
else:
    min_q = None
    max_q = None


if savedata:
	print("Saving data...")
	if outform == 'npy':
		save(toplot+'_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+suffix,dataplot)
	else:
		hdu = fits.PrimaryHDU(dataplot)
		hdu.writeto(toplot+'_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+suffix+'.fits')

if verbose:
	print(shape(dataplot))

print("Done")

############### 2D IMAGE ###############
if img2d:
    print('Creating the 2d figure...')
    fig = figure(1)    #create the figure
    ax1 = fig.add_subplot(111)
	
    ax1.set_xlim([0.0,Lbox])
    ax1.set_ylim([0.0,Lbox])

    ax1.set_title(r'$z = $'+str(red),fontsize=14) #+', depth = '+str(round(width*gridtompc,2))+' cMpc',fontsize=20)
	
    if plane == 'XY':
        ax1.set_xlabel(r'X [cMpc]',fontsize=14)  #x-axis label
        ax1.set_ylabel(r'Y [cMpc]',fontsize=14)  #y-axis label
    elif plane == 'XZ':
        ax1.set_xlabel(r'X [cMpc]',fontsize=14)  #x-axis label
        ax1.set_ylabel(r'Z [cMpc]',fontsize=14)  #y-axis label
    else:
        ax1.set_xlabel(r'Y [cMpc]',fontsize=14)
        ax1.set_ylabel(r'Z [cMpc]',fontsize=14)	
    
    if min_q is None:  min_q = amin(dataplot)	
    if max_q is None:  max_q = amax(dataplot)
    
    if scale=='linear':
        cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)
    else:
        cax = ax1.imshow(dataplot,cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))
    
    cbar = fig.colorbar(cax)
    if toplot == 'SB':
        cbar.set_label(toplot+r" $[\rm{erg} \rm{s}^{-1} \rm{cm}^{-2} \rm{arcsec}^{-2}]$",fontsize=14)
    elif toplot == 'temperature':
        cbar.set_label(toplot+" [K]",fontsize=14)
    elif toplot == 'density':
        #cbar.set_label(r'$\rho / \bar{\rho}_b(z=0)$',fontsize=20)
        cbar.set_label(r'$\rho$ [atom/cm$^{-3}]$',fontsize=14)
    else:
        cbar.set_label(toplot,fontsize=20)
	
    if plotRH:
        Rh = 0.6968
        c  = [Lbox/2.,Lbox/2.,Lbox/2.]
        if plane == 'XY':
            circle1  = plt.Circle((c[0], c[1]), Rh, color='r', fill=False,alpha=0.5)
            circle07 = plt.Circle((c[0], c[1]), 0.7*Rh, color='r', ls ='--',fill=False,alpha=0.5)
            circle05 = plt.Circle((c[0], c[1]), 0.5*Rh, color='r', ls ='--',fill=False,alpha=0.5)
            circle02 = plt.Circle((c[0], c[1]), 0.2*Rh, color='r', ls ='--',fill=False,alpha=0.5)
        elif plane == 'XZ':
            circle = plt.Circle((c[0], c[2]), Rh, color='r', fill=False)
        else:
            circle = plt.Circle((c[1], c[2]), Rh, color='r', fill=False)
        ax1.add_patch(circle1)
        ax1.add_patch(circle07)
        ax1.add_patch(circle05)
        ax1.add_patch(circle02)
    #fout = toplot+'_z'+str(red)+'_Halo4_'+plane+'_w'+str(width)+suffix+'.pdf'
    #fout = toplot+'_z'+str(red)+'_'+plane+'_w'+str(width)+suffix+'.pdf'
    #print("saving figure to "+fout)
    #savefig(fout, bbox_inches='tight')
    plt.show()
    close(fig)

############### DIFFERENCE ###############
if diff:
    d1md2 = (dataplot-dataplot2)/dataplot2
    print('Creating the 2d figure of difference...')    

    fig = figure(3)		
    ax2 = fig.add_subplot(111)
    
    ax2.set_xlim([0.0,Lbox])
    ax2.set_ylim([0.0,Lbox])
	
    ax2.set_title(r'$z = $'+str(red)+', depth = '+str(round(width*gridtompc,2))+' Mpc',fontsize=20)
	
    if plane == 'XY':
        ax2.set_xlabel(r'$X\, [Mpc]$',fontsize=18)
        ax2.set_ylabel(r'$Y\,  [Mpc]$',fontsize=18) 
    elif plane=='XZ':
        ax2.set_xlabel(r'$X\, [Mpc]$',fontsize=18)
        ax2.set_ylabel(r'$Z\, [Mpc]$',fontsize=18)
    else:
        ax2.set_xlabel(r'$Y\,  [Mpc]$',fontsize=18)
    ax2.set_ylabel(r'$Z\, [Mpc]$',fontsize=18)
    
    min_q = -2. #amin(d1md2)
    max_q = 2. #amax(d1md2)
   
    scale = 'linear'
    if scale=='linear':
        cax = ax2.imshow(d1md2,cmap=get_cmap('bwr'),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)       

    else:
        cax = ax2.imshow(abs(d1md2),cmap=get_cmap(cmap),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))
    
    cbar = fig.colorbar(cax)
    cbar.set_label('difference '+ toplot,fontsize=20)
    
    fout2 = toplot+'_diff_z'+str(red)+'_gn'+str(gn)+'_'+plane+'_w'+str(width)+suffix+'.pdf'
    print("saving figure to "+fout2)
    savefig(fout2, bbox_inches='tight')
    plt.show()
    close(fig)

############### HIST ###############
if hist:
	print('Creating the histogram...')
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
