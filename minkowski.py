from numpy import*
import h5py, sys, os
import argparse
from astropy.io import fits
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import *
from minkfncts2d import MF2D

# input parameters
parser = argparse.ArgumentParser(description='Routine to compute minkowski functionals of SB maps.')
parser.add_argument('DM',nargs=1,type=str, help='Type of DM (CDM, WDM_3p5)')
parser.add_argument('zoom',nargs=1,type=str,help='Zoom level (1, 1p5, 2, 3)')
parser.add_argument('plane',nargs=1,type=str, help='Plane to plot (XY, YZ, XZ)')
parser.add_argument('width',nargs=1,type=int, help='Width of slice (grid cells)')
parser.add_argument('nvox',nargs=1,type=int, help='CubEx NVox (int)')
parser.add_argument('smooth',nargs=1,type=int, help='Sigma of smoothing kernel (int)')
parser.add_argument('red',nargs=1,type=float, help='Redshift (float)')
args   = parser.parse_args()
dm     = args.DM[0]
zoom   = args.zoom[0]
plane  = args.plane[0]
width  = args.width[0]
nvox   = args.nvox[0]
s      = args.smooth[0]
red    = args.red[0]

plot = False

threshold = 0.0 #in cubex anything above 0 is an object
cut       = [1e-20, 5e-20, 1e-19,5e-19]
if dm == 'WDM':
	root = '/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/'+dm+'_3p5/zoom'+zoom+'_hydro/'
else:
	root = '/cobra/ptmp/titouan/Eagle_cosmo_arepo/200_512/'+dm+'/zoom'+zoom+'_hydro/'
F, U, Chi = [],[],[]
fi = open(root+"minkowski_z"+str(int(red))+'_Halo4_'+plane+'_w'+str(width)+'_s'+str(s)+'p_NVox'+str(nvox)+'.dat','w')
for c in cut:
	fname = 'SB_z'+str(int(red))+'_Halo4_'+plane+'_w'+str(width)+'_s'+str(s)+'p_zoom'+zoom+'_'+dm+'.Objects_Id_cut'+str(c)+'_NVox'+str(nvox)+'.fits'
	hdu   = fits.open(root+fname)
	SB    = hdu[0].data[0,:,:]
	SB    = SB.newbyteorder()
	SB    = SB.astype('float64')
	
	#minkowski functionals
	(f,u,chi) = MF2D(SB, threshold)
	F.append(f)
	U.append(u)
	Chi.append(chi)
	
	#write
	fi.write(str(c)+"\t"+str(f)+"\t"+str(u)+"\t"+str(chi)+"\n")
fi.close()

#check plot
if plot:
	plt.yscale('log')
	plt.plot(cut, F, cut, U, cut, Chi)
	plt.title("2D Minkowski Functions")
	plt.legend(["F (Area)", "U (Boundary)", "$\chi$ (Euler characteristic)"])
	plt.xlabel("Threshold")
	plt.show()
