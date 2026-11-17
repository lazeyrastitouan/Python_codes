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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

tmin   = None
tmax   = None
suffix = ''

# input parameters
parser = argparse.ArgumentParser(description='Routine to plot P2C output grids. Args: < hdf5 file (str) >, < optional -- second hdf5 file (str)>, < suffix of outputfile (str) >, < what to plot (str) >, < plane to plot (str: XY, YZ, XZ)>, < width of slice (int, in cells)>, <optional -- max T (float, K)>, <optional -- min T (float, K)>')
parser.add_argument('hdf5file',nargs=1,type=str, help='P2C output (hdf5)')
parser.add_argument('--suffix',nargs=1,type=str, help='Suffix of outputfile (str, optional)')
parser.add_argument('plane',nargs=1,type=str, help='Plane to plot')
parser.add_argument('--tempmin',nargs=1,type=float, help='Minimum temperature (float, optional)')
parser.add_argument('--tempmax',nargs=1,type=float, help='Maximum temperature (float, optional)')
args   = parser.parse_args()
path   = args.hdf5file[0]
if args.suffix is not None:
	suffix = '_'
	suffix += args.suffix[0]
plane  = args.plane[0]
if args.tempmin is not None:
	tmin   = args.tempmin[0]
if args.tempmax is not None:
	tmax   = args.tempmax[0]

#some hard-coded settings
red           = 3.0 #3.596 
a             = 1./(1.+red)
h             = 0.6777 #0.6710917274583887
Lbox          = 2.0 #cMpc
ave_grid_mass = 2.93136e3
scale         = 'log' #'linear' or 'log'

# what to do
verbose  = True
img2d    = True
savedata = False 

print("2D image: "+str(img2d))
print("Saving data: "+str(savedata))

if savedata == True:
	outform = 'npy' #'fits' or 'npy'
	print("Output data format: "+outform)

print("Reading data...")

#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"

# open the P2C file(s)
f = h5py.File(path, "r")

#grid size in 1D
sideL     = array(f['level_0/boxes'])[0][3]+1
cellinmpc = Lbox/sideL #cMpc

# load mass/ave_grid_mass ("density", 1D array to re-order)
dset = f[CPATH]
data  = dset[0:sideL**3].reshape((sideL, sideL, sideL), order="F")

#if temperature boundaries imposed, load temperature and put to 0 all values not in interval
if tmin is not None or tmax is not None:
	temp  = dset[4 * sideL ** 3 : (5 * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
if tmin is not None:
	data  = where(temp<tmin,0.0,data)
if tmax is not None:
	data  = where(temp>=tmax,0.0,data)

#density is in fact M/ave_grid_mass so we convert to atom/cm^-2
mH         = 1.6735575e-24 #gram
mHe        = 6.6464731e-24 #gram
mA         = 0.76*mH+0.24*mHe
Msuningram = 1.989e33
Mpcincm    = 3.086e24
data       = data * ave_grid_mass / (a*cellinmpc*Mpcincm)**3. * Msuningram/mA

if verbose:
	print(shape(data))

#compute the column density for selected plane
if verbose:
	print("Cells summed: 0:"+str(int(sideL)))
if plane == 'XY':
    Ngas = sum(data,axis=2)
elif plane == 'YZ':
    Ngas = sum(data,axis=0)
elif plane == 'XZ':
    Ngas = sum(data,axis=1)
else:
	sys.exit("Selected plane not valid. Options are XY, YZ, or XZ") 

#Multiply by LoS in pcm (since we summed, Los is just 1 cell, and sum creates the full LoS
Ngas = Ngas*cellinmpc*a*Mpcincm

if savedata:
	print("Saving data...")
	if outform == 'npy':
		save('/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/2cMpc/Ngas_z'+str(red)+'_'+plane+suffix,Ngas)
	else:
		hdu = fits.PrimaryHDU(dataplot)
		hdu.writeto('/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/2cMpc/Ngas_z'+str(red)+'_'+plane+suffix+'.fits')

if verbose:
    print("maximum col density [cm^-2]: {}".format(amax(Ngas)))
    print("minimum col density [cm^-2]: {}".format(amin(Ngas)))
print("Done")

############### 2D IMAGE ###############
if img2d:
    print('Creating the 2d figure...')
    fig = figure(1)    #create the figure
    ax1 = fig.add_subplot(111)
	
    ax1.set_xlim([0.0,Lbox])
    ax1.set_ylim([0.0,Lbox])

    '''
    ax1.set_title(r'$z = $'+str(red),fontsize=14)
	
    if plane == 'XY':
        ax1.set_xlabel(r'X [cMpc]',fontsize=14)  #x-axis label
        ax1.set_ylabel(r'Y [cMpc]',fontsize=14)  #y-axis label
    elif plane == 'XZ':
        ax1.set_xlabel(r'X [cMpc]',fontsize=14)  #x-axis label
        ax1.set_ylabel(r'Z [cMpc]',fontsize=14)  #y-axis label
    else:
        ax1.set_xlabel(r'Y [cMpc]',fontsize=14)
        ax1.set_ylabel(r'Z [cMpc]',fontsize=14)	
    '''
    min_q = 1e19 #amin(Ngas)	
    max_q = 1e22
    
    if scale=='linear':
        cax = ax1.imshow(Ngas.T,cmap=get_cmap('viridis'),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',vmin=min_q,vmax=max_q)
    else:
        cax = ax1.imshow(Ngas.T,cmap=get_cmap('viridis'),origin='lower',extent=[0.0,Lbox,0.0,Lbox],interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))
    cax.axes.get_xaxis().set_ticks([])
    cax.axes.get_yaxis().set_ticks([])
    ### inset colorbar ###
    #cbaxes = inset_axes(ax1, width="5%", height="50%", loc="upper left") 
    #cbar   = fig.colorbar(cax,cax=cbaxes, orientation='vertical',ticks=[1e19,1e20,1e21,1e22])
    #cbar.ax.tick_params(color='white',which='both', labelcolor='white')
    ### horizontal colorbar ###
    cbar = fig.colorbar(cax,orientation='vertical',fraction=0.043, pad=0.04)
    cbar.set_label(r'N [atom/cm$^{-2}]$',fontsize=12)
    fout = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/2cMpc/Ngas_z'+str(red)+'_'+plane+suffix+'.pdf'
    print("Saving figure to "+fout)
    plt.savefig(fout,bbox_inches='tight')
    plt.show()
    close(fig)
