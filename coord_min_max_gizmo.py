from numpy import*
import h5py
import os
import sys


folder='/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2p5_hydro_fire2_lvl7/filament/snapdir_017/'
h5files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

xmin=135540.;ymin=135540.;zmin=135540.
xmax=0.0;ymax=0.0;zmax=0.0

for i in h5files:
    if i[0:4] == 'snap' and (i[-4:] == 'hdf5' or i[-4:] == 'HDF5'):
        f = h5py.File(folder+i, "r")
    else:
        continue

    x=array(f['/PartType0/Coordinates'])
    xmint=amin(x[:,0])
    ymint=amin(x[:,1])
    zmint=amin(x[:,2])
    xmaxt=amax(x[:,0])
    ymaxt=amax(x[:,1])
    zmaxt=amax(x[:,2])

    if xmint < xmin:
        xmin=xmint
    if ymint < ymin:
        ymin=ymint
    if zmint < zmin:
        zmin=zmint
    if xmaxt > xmax:
        xmax=xmaxt
    if ymaxt > ymax:
        ymax=ymaxt
    if zmaxt > zmax:
        zmax=xmaxt
    print('Min: '+str(xmin)+', '+str(ymin)+', '+str(zmin))
    print('Max: '+str(xmax)+', '+str(ymax)+', '+str(zmax))

lx = xmax-xmin
ly = ymax-ymin
lz = zmax-zmin

Lbox=max(array([lx,ly,lz]))
middle=array([(xmin+xmax)/2.,(ymin+ymax)/2.,(zmin+zmax)/2.])
cutle=(middle-Lbox/2.)/1e3
cutre=(middle+Lbox/2.)/1e3

print('New box size: '+str(max(array([lx,ly,lz]))))

print('For P2C: cutle '+str(cutle)+', cutre '+str(cutre)+', new center '+str(middle))
