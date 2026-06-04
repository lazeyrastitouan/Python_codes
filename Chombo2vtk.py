from numpy import*
import h5py, sys, os
import argparse
from astropy.io import fits
from pyevtk.hl import gridToVTK,imageToVTK

def data2vtk(xtab,ytab,ztab,delta,rdata_f_name):

    print("x dimension = "+str(shape(xtab))+", y dimension = "+str(shape(ytab))+", z dimension = "+str(shape(ztab)))
    print("shape delta = "+str(shape(delta)))

    gridToVTK(rdata_f_name, xtab, ytab, ztab, pointData = { toplot[0:4] : delta})

def find_index(i):
	switcher={'density':1, 'velocity_x':2, 'velocity_y':3, 'velocity_z':4, 'temperature':5, 'stars':6, 'emissivity':7, 'SB':7, 'density_weighted_velz':8, 'density_weighted_vely':9, 'density_weighted_velz':10}
	
	return switcher.get(i,"Invalid quantity to plot\n Choices are density, velocity_x, velocity_y, velocity_z, temperature, stars, emissivity (or SB), density_weighted_velx, density_weighted_vely, density_weighted_velz")

# input parameters
parser = argparse.ArgumentParser(description='Convert P2C output grids to vtk. Args: < hdf5 file (str) >, < what to convert (str) >')
parser.add_argument('hdf5file',nargs=1,type=str, help='P2C output (hdf5)')
parser.add_argument('toplot',nargs=1,type=str, help='Quantity to convert')
args   = parser.parse_args()
path   = args.hdf5file[0]
toplot = args.toplot[0]

# what to do
verbose  = True
Tfloors  = False

#some hard-coded settings
Lbox = 4.0 #cMpc

if verbose:
    print("Box size = "+str(Lbox)+" cMpc")

print("Reading data...")

#indices for the various attributes: 1=den, 2=velx, 3=vely, 4=velz, 5=temp, 6=stars, 7=emissivity, 8=den**2weighted velx, 9=den**2weighted vely, 10=den**2weighted velz
nvar  = 10
CPATH = "/level_0/data:datatype=0"
Num   = find_index(toplot)

# open the P2C file(s)
f = h5py.File(path, "r")

#grid size in 1D
sideL = array(f['level_0/boxes'])[0][3]+1
xtab  = linspace(0,Lbox,sideL)
ytab  = linspace(0,Lbox,sideL)
ztab  = linspace(0,Lbox,sideL)

# navigate to the relevant dataset (1D array to re-order)
dset = f[CPATH]
# load data
# NB: python does not include the last index in [i:j]
# NB: the reshaping has to use Fortran-like indexing as the chombo-file was originally produced using Fortran-like reshaping
data  = dset[(Num - 1) * sideL ** 3 : (Num * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
if toplot=='emissivity':
    Elya      = 1.6342236e-11 #[erg] 
    ave_b_rho = 6.76991e-31 #[g/cm^-3]
    alpha     = 1.269e-13 #[cm^-3/s]
    Y         = 0.24 #[-]
    mH        = 1.673557e-24 #[g]
    data      *= Elya * ave_b_rho**2. * alpha * (1.-0.5*Y)/(1.-Y)/mH**2. #[erg/s/cm^-3] 

if Tfloors:
    T1   = 4.5
    T2   = 5.5
    temp = dset[(5 - 1) * sideL ** 3 : (5 * sideL ** 3)].reshape((sideL, sideL, sideL), order="F")
    temp = log10(array(temp))

if verbose:
    print("Dividing cold, cool and hot gas: "+str(Tfloors))
    if Tfloors:
        print("Using T1 = "+str(T1)+" and T2 = "+str(T2))

print("Saving data...")
vtkfilename = path[:-5]+'_'+toplot
xtab = linspace(0,Lbox,sideL)
ytab = linspace(0,Lbox,sideL)
ztab = linspace(0,Lbox,sideL)
if Tfloors:
    vtkfilename1=vtkfilename+'_cold'
    vtkfilename2=vtkfilename+'_cool'
    vtkfilename3=vtkfilename+'_hot'
    data2vtk(xtab,ytab,ztab,where(temp<T1,data,1e-10),vtkfilename1)
    data2vtk(xtab,ytab,ztab,where((temp>=T1)&(temp<=T2),data,1e-10),vtkfilename2)
    data2vtk(xtab,ytab,ztab,where(temp>T2,data,1e-10),vtkfilename3)
else:
    data2vtk(xtab,ytab,ztab,data,vtkfilename)

print("Done")
