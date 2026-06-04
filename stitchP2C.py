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


# input parameters
parser = argparse.ArgumentParser(description='Routine to plot P2C output grids. Args: < path to hdf5 files (str) >')
parser.add_argument('path',nargs=1,type=str, help='path to P2C output (hdf5)')
args   = parser.parse_args()
path   = args.path[0]
h5files  = [f for f in os.listdir(path+"/") if os.path.isfile(os.path.join(path+"/", f)) and f[-5:] == '.hdf5']
h5files  = sorted(h5files)

verbose  = True

print("Reading data...")
#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = [0,1,2,3,4,5,6,7,8,9]
ln    = len(nvar)
CPATH = "/level_0/data:datatype=0"
datasub = [[] for x in range(ln)]

# open the P2C file(s) (loaded in order from gn1 to gn8)
for name in h5files:
    f = h5py.File(path+'/'+name, "r")

    # navigate to the relevant dataset (1D array to re-order)
    dset = f[CPATH]

    sidesubL  = array(f['level_0/boxes'])[0][3]+1
	
    for Num in nvar:
        ni = nvar.index(Num)
		# load quantity into the array "data"
		# NB: python does not include the last index in [i:j]
		# NB: the reshaping has to use Fortran-like indexing as the chombo-file was originally produced using Fortran-like reshaping
        datasub[ni].append(dset[Num * sidesubL ** 3 : ((Num + 1) * sidesubL ** 3)].reshape((sidesubL, sidesubL, sidesubL), order="F")) 
    
    f.close()

print("Done")

#now we stich in z, y, x directions
#datasub gn [1,2,3,4,5,6,7,8] or [[+++],[++-],[+-+],[+--],[-++],[-+-],[--+],[---]]
#x:vstack, y:hstack, z:dstack
print("Stitching boxes...")
data = []

for Num in range(ln):
	#stich at "x=+"
	data00 = dstack((datasub[Num][1],datasub[Num][0]))  
	data01 = dstack((datasub[Num][3],datasub[Num][2]))
	data0  = hstack((data01,data00))

	del data00, data01

	#stich at "x=-"
	data10 = dstack((datasub[Num][5],datasub[Num][4]))
	data11 = dstack((datasub[Num][7],datasub[Num][6]))
	data1  = hstack((data11,data10))

	del data10, data11

	#final stich
	data.append(vstack((data1,data0)))
	del data0, data1

data   = array(data)
sideL  = shape(data)[1]

if verbose:
	print("Shape data before reshaping "+str(shape(data)))

	dataplot = empty((1024,1024))
	for i in range(sideL):
		for j in range(sideL):
			dataplot[i,j]=sum(data[0,j,i,(int(1024/2)-int(1024/2)):(int(1024/2)+int(1024/2))])/1024.
	
	print('Creating the 2d figure...')
	fig = figure(1)    #create the figure
	ax1 = fig.add_subplot(111)

	min_q = amin(dataplot)
	max_q = amax(dataplot)

	ax1.set_xlabel('X [Mpc]',fontsize=18)  #x-axis label
	ax1.set_ylabel('Y [Mpc]',fontsize=18)  #y-axis label
	cax = ax1.imshow(dataplot,cmap='viridis',origin='lower',interpolation='none',norm = LogNorm(vmin=min_q,vmax=max_q))
	cbar = fig.colorbar(cax)
	plt.show()

print("Done")
print("Reshaping...")

datars = empty(int(ln*sideL**3))
for i in range(ln):
	datars[i*sideL**3:(i+1)*sideL**3] = reshape(data[i],int(sideL**3),order='F')

if verbose:
	print("Shape data after reshaping "+str(shape(data)))

print("Done")

#rewrite
print("Writing")
#reopen the first file
f = h5py.File(path+'/'+h5files[0],'r')

#outfile
of = h5py.File(path+'/P2C_grid_'+str(int(sideL))+'.hdf5','w')

#groups
groups = ['Chombo_global','level_0']
for g in groups:
	of.create_group(g)

#attrs
for a in f['Chombo_global'].attrs.items():
	of['Chombo_global'].attrs.create(a[0],a[1])
for a in f['level_0'].attrs.items():
	of['level_0'].attrs.create(a[0],a[1]) #change dx and prob_domain

#datasets
for ds in delete(list(f['level_0'].keys()),(4,5)):
	of['level_0'].create_dataset(ds,data=f['level_0/'+ds][0])
of['level_0'].create_dataset('boxes',data=array([(0, 0, 0, sideL-1, sideL-1, sideL-1)],dtype=[('lo_i', '<i4'), ('lo_j', '<i4'), ('lo_k', '<i4'), ('hi_i', '<i4'), ('hi_j', '<i4'), ('hi_k', '<i4')]))
of['level_0'].create_dataset('data:datatype=0',data=datars)

of.close()
f.close()	
