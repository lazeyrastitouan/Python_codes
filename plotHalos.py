from numpy import*
from astropy.io import fits
from pylab import *
from matplotlib.colors import LogNorm

cat    = loadtxt("/data/titouan/Eagle_cosmo_gizmo/100_1024/WDM_3p5keV/base/AHF_002/200_snap_z3.000_halos_sorted",usecols=(5,6,7)).T/0.6777
hdul   = fits.open("/data/titouan/Eagle_cosmo_gizmo/100_1024/WDM_3p5keV/base/AHF_002/200_snap_z3.000_halos_sorted_density_100.0Mpc_Ng256.fits")
dh     = hdul[0].data
Lbox   = 100.
width  = Lbox/256.*1000.
filter = cat[2]<=width
toplot = dh[:,:,0]
mindh  = amin(toplot)
maxdh  = amax(toplot)

#check density field
fig=figure()
ax1=fig.add_subplot(111)
cax=ax1.imshow(toplot.T,cmap=get_cmap('viridis'),origin='lower',extent=[0,Lbox,0,Lbox],interpolation='bicubic',norm=LogNorm(vmin=1e3*mindh,vmax=maxdh))
ax1.plot(cat[0][filter]/1000.,cat[1][filter]/1000,'ko')
cbar = fig.colorbar(cax)
plt.show()

#Histogram
dh1D = dh.flatten()
hist(dh1D,bins=logspace(log10(mindh),log10(maxdh),100))#,stacked=True,density=True)
gca().set_xscale('log')
plt.show()
