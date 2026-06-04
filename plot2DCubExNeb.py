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

def obs_to_int(x):
    return x*(1+red)**4.

def int_to_obs(x):
    return x/(1+red)**4.

# input parameters
parser = argparse.ArgumentParser(description='Routine to plot CubEx projected 2D nebulae. Args: < Fits file (str) >, < Grid number (int) >')
parser.add_argument('path',nargs=1,type=str, help='Fits file')
parser.add_argument('gn',nargs=1,type=int, help='Grid number')
args   = parser.parse_args()
path   = args.path[0]
gn     = args.gn[0]

#some hard-coded settings
red   = 3.5
Lbox  = 2. #Mpc
width = 256
plane = 'XY'
cmap  = 'jet' #'cividis' #'nipy_spectral',#'Spectral_r',#'Set3',#'Paired' #'jet'
print("Redshift: "+str(red)+"\nPlane: "+plane+"\nLbox: "+str(Lbox)+"\nWidth: "+str(width))
print("Reading data...")

dhu   = fits.open(path)
data  = dhu[0].data
data  = data*(1.+red)**4.

print(shape(data))
sideL     = shape(data)[0]
gridtompc = Lbox/sideL

print("Done")

############### 2D IMAGE ###############
print('Creating the 2d figure...')

min_q     = 5e-20*(1.+red)**4. 
max_q     = 1e-13#/(1.+red)**4.  

fig = figure(1)    #create the figure
ax1 = fig.add_subplot(111)

ax1.set_xlim([0.0,2.0])#Lbox])
ax1.set_ylim([0.0,2.0])#Lbox])

#ax1.set_title(r'$z = $'+str(red),fontsize=14)
	
#ax1.set_xlabel(r'$[cMpc]$',fontsize=18)  #x-axis label
#ax1.set_ylabel(r'$[cMpc]$',fontsize=18)  #y-axis label

#if plane == 'XY':
#	ax1.set_xlabel(r'$X\, [Mpc]$',fontsize=18)  #x-axis label
#	ax1.set_ylabel(r'$Y\,  [Mpc]$',fontsize=18)  #y-axis label
#elif plane == 'XZ':
#	ax1.set_xlabel(r'$X\, [Mpc]$',fontsize=18)  #x-axis label
#	ax1.set_ylabel(r'$Z\,  [Mpc]$',fontsize=18)  #y-axis label
#else:
#	ax1.set_xlabel(r'$Y\,  [Mpc]$',fontsize=18)
#	ax1.set_ylabel(r'$Z\,  [Mpc]$',fontsize=18)	

cax  = ax1.imshow(data,cmap=get_cmap(cmap),origin='lower',extent=[0.0,2.0,0.0,2.0],interpolation='none',norm=LogNorm(vmin=min_q,vmax=max_q))
cax.axes.get_xaxis().set_ticks([])
cax.axes.get_yaxis().set_ticks([])
circle1  = plt.Circle((1, 1), 0.1, color='k', lw=2,fill=False)
circle2  = plt.Circle((1, 1), 0.2, color='k', lw=2,fill=False)
#ax1.add_patch(circle1)
#ax1.add_patch(circle2)
### inset colorbar ###
#cbaxes = inset_axes(ax1, width="3%", height="40%",loc="upper left")  #for Halo10 gn3: loc="lower left", bbox_to_anchor=(0.1,0,1,1), bbox_transform=ax1.transAxes)
#cbar   = fig.colorbar(cax,cax=cbaxes, orientation='vertical')#,ticks=[1e19,1e20,1e21,1e22])
#cbar.ax.tick_params(labelsize=8) #which='both',labelsize=8,right=False,left=True) #color='white',which='both', labelcolor='white')
#cbaxes.yaxis.set_label_position('right')
#cbaxes.yaxis.set_ticks_position('right') #'left' for for Halo10 gn3
#cbar.set_label(r'[erg/s/cm$^2$/arcsec$^2$]',fontsize=8)
### outside colorbar ###
cbar = fig.colorbar(cax, orientation='horizontal',fraction=0.043, pad=0.04)#,pad = 0.3)
cbar.set_label(r'SB$_{\rm{intr}}$ [erg/s/cm$^2$/arcsec$^2$]',fontsize=12)
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.xaxis.set_label_position('bottom')
fout = path[:-5]+'_bottomcolorbar.pdf'
print("Saving figure to "+fout)
savefig(fout, bbox_inches='tight')
plt.show()
close(fig)
