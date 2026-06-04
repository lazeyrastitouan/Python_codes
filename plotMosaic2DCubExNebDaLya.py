from numpy import*
import h5py, sys, os
import argparse
from astropy.io import fits
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import *
#from matplotlib.ticker import ScalarFormatter
#import matplotlib.gridspec as gridspec
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import LogNorm
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe

#some hard-coded settings
red   = 3.5
h     = 0.6777
Lbox  = 2. #Mpc
Lfrac = 0.5 #fraction of box plotted
shift = Lbox*Lfrac/2. #shift of (0,0)
width = 256
cmap  = 'jet' #'nipy_spectral' 'cividis' #'nipy_spectral','Spectral_r','Set3','Paired'
print("Redshift: "+str(red)+"\nLbox: "+str(Lbox)+"\nWidth: "+str(width))
print("Fraction of box plotted: "+str(Lfrac)+"\nCorresponding to: "+str(int(width*Lfrac))+" cells")
print("Reading data...")

rootH10   = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/2cMpc/'
rootH1312 = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/2cMpc/'
arcsec    = 0.034 #cMpc
x0, y0    = (Lbox-0.4)*Lfrac+shift, 0.1*Lfrac+shift

rootsD  = array([[rootH1312,rootH10,rootH10],[rootH1312,rootH10,rootH10],[rootH10,rootH10,rootH10]])
losD    = array([['X','Y','Z'],['Z','Y','X'],['Y','X','Z']])
gnD     = array([[3,4,6],[4,7,8],[9,10,11]])
MD      = array([[12.47,12.41,12.39],[12.28,12.23,12.18],[12.18,12.14,12.03]])

min_q   = 5e-20*(1.+red)**4.
max_q   = 1e-13

#rescale 0.2,1,1.8 to 0 at center and +- on left and right
lbls    = [int((0.2*Lfrac+shift-Lbox/2.)*1e3),int((Lbox/2.*Lfrac-shift)*1e3),int(((Lbox-0.2)*Lfrac-shift)*1e3)]  

print('Creating the 2d figure for DaLya...')
fig1, axs1 = plt.subplots(3,3,sharey=True,sharex=True,figsize=(7,7))
plt.subplots_adjust(hspace=0.01,wspace=0.01)
for i in range(shape(axs1)[0]):
    for j in range(shape(axs1)[1]):
        mockCube = rootsD[i,j]+'mockCube'+losD[i,j]+'_noise5e-20_gn'+str(gnD[i,j])+'_SB_s1_z3.5_Ng256_observed_cubex_neb_singleobj_NB_nonoise.fits' 
        dhu   = fits.open(mockCube)
        data  = dhu[0].data
        data  = data*(1.+red)**4.
        sideL     = shape(data)[0]
        gridtompc = Lbox/sideL
        ilow  = int(0+round(shift/gridtompc,0)) #lower index of data compatible with what we want to plot
        ihigh = int(sideL-round(shift/gridtompc,0)) #higher index of data compatible with what we want to plot

        axs1[i,j].set_xlim([0.0+shift,Lbox-shift])
        axs1[i,j].set_ylim([0.0+shift,Lbox-shift])
        axs1[i,j].text(0.05*Lfrac+shift,(Lbox-0.2)*Lfrac+shift,r'$\log M = $'+str(MD[i,j]),color='w',fontsize=11,path_effects=[pe.Stroke(linewidth=2, foreground="black"),pe.Normal()])	
        axs1[i,j].set_xticks([0.2*Lfrac+shift,Lbox/2.,Lbox-0.2*Lfrac-shift],labels=[r''+str(lbls[0]),r''+str(lbls[1]),r''+str(lbls[2])],fontsize=12)#,minor=True) #in units of data, so Lbox/2 in center, but we put the label corresponding to zoom. For the left and right ticks we want 0.25 rescaled by Lfrac from 0 and Lboc, and we have to move by shift
        axs1[i,j].set_yticks([0.2*Lfrac+shift,Lbox/2.,Lbox-0.2*Lfrac-shift],labels=[r''+str(lbls[0]),r''+str(lbls[1]),r''+str(lbls[2])],fontsize=12)#,minor=True)

        cax     = axs1[i,j].imshow(data[ilow:ihigh+1,ilow:ihigh+1],cmap=get_cmap(cmap),origin='lower',extent=[0.0+shift,Lbox-shift,0.0+shift,Lbox-shift],interpolation='none',norm=LogNorm(vmin=min_q,vmax=max_q))
        circle1 = plt.Circle((1, 1), 0.1, color='white', lw=1.5,fill=False)
        circle1.set_path_effects([pe.Stroke(linewidth=2.5, foreground="black"),pe.Normal()])
        circle2 = plt.Circle((1, 1), 0.2, color='white', lw=1.5,fill=False)
        circle2.set_path_effects([pe.Stroke(linewidth=2.5, foreground="black"),pe.Normal()])
        axs1[i,j].add_patch(circle1)
        axs1[i,j].add_patch(circle2)

        #10*Lfrac arcsec bar
        axs1[i,j].plot([x0, x0 + 10*Lfrac*arcsec],[y0, y0],color="white",solid_capstyle="butt")
        axs1[i,j].text(x0 + 10*Lfrac*arcsec / 2,y0 + 0.2*Lfrac,''+str(int(10*Lfrac))+'″',color="white",ha="center",va="top")

fout = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/mosaic_mockCubes_DaLya_n5e-20_z3p5_intr_cubex_neb_singleobj_NB_nonoise_new.pdf'
print("Saving DaLya mosaic to "+fout)
savefig(fout, bbox_inches='tight')
plt.show()
close(fig1)
