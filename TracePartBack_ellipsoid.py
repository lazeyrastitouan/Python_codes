from numpy import*
import sys
import os
import h5py
import argparse
from multiprocessing import Pool
import morph

# What to do
parser = argparse.ArgumentParser(description='Trace particles of a cubic region of a simulation centered on a given position back to initial condition. Args: < path to snapshot (hdf5) (str) >, < initial conditions file (str) >, < xc, yc, zc (floats) >, < zoom size (float) >')
parser.add_argument('snappath',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('icspath',nargs=1,type=str, help='Path to initial conditions (hdf5)')
parser.add_argument('xyz',nargs=3,type=float, help='Center of region xc, yc, zc [kpc/h]')
parser.add_argument('size',nargs=1,type=float, help='New box size [kpc/h]')
args     = parser.parse_args()
snappath = args.snappath[0]
icspath  = args.icspath[0]
xc 	     = args.xyz[0]
yc 	     = args.xyz[1]
zc 	     = args.xyz[2]
size     = args.size[0]
xmin     = xc-size/2.
ymin     = yc-size/2.
zmin     = zc-size/2.
xmax     = xc+size/2.
ymax     = yc+size/2.
zmax     = zc+size/2.
safety   = 0.0 #size/25.
ptypes   = [0,1,2]
ntypes   = len(ptypes)
verbose  = False

h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5']
h5files  = sorted(h5files)
h5files  = sorted(h5files,key=len)
if verbose:
	print(h5files)

print("Keeping particule in cube\n"+str(xmin/1000.)+" < x < "+str(xmax/1000.)+" Mpc/h\n"+str(ymin/1000.)+" < y < "+str(ymax/1000.)+" Mpc/h\n"+str(zmin/1000.)+" < z < "+str(zmax/1000.)+" Mpc/h\n")
print("Size of safety region: "+str(safety/1000.)+" Mpc/h")
 
#1 Find particle ids in a given snapshot region
print("Selecting particles of type "+str(ptypes)+" in volume...")
ids=[[] for x in range(ntypes)]
for i in h5files:
    f = h5py.File(snappath+"/"+i, "r")
    
    for p in ptypes:
        ip   = ptypes.index(p)
        idst = array(f['PartType'+str(p)+'/ParticleIDs'])
        pos  = array(f['PartType'+str(p)+'/Coordinates'])

        numparts = f['Header'].attrs['NumPart_ThisFile'][p]	
        for j in range(numparts):
            if xmin <= pos[j,0] <= xmax and ymin <= pos[j,1] <= ymax and zmin <= pos[j,2] <= zmax:		
                ids[ip].append(idst[j])

npin = sum(len(ids[i]) for i in range(ntypes))
print("Done. Found "+str(npin)+" particles")

#2 Find the corresponding volume occupied by these particles in the ICs
print("Finding corresponding region in ICs...")

h5ics        = [f for f in os.listdir(icspath+"/") if os.path.isfile(os.path.join(icspath+"/", f)) and f[-4:] == 'hdf5']
h5ics        = sorted(h5ics)
h5ics        = sorted(h5ics,key=len)
idsic, posic = [[] for x in range(ntypes)],[[] for x in range(ntypes)]
box          = h5py.File(icspath+'/'+h5ics[0],'r')['Header'].attrs['BoxSize']
for ic in h5ics:
    f = h5py.File(icspath+"/"+ic,'r')

    for p in ptypes:
        ip   = ptypes.index(p)
        idsict = array(f['PartType'+str(p)+'/ParticleIDs'])
        posict = array(f['PartType'+str(p)+'/Coordinates'])

        numparts = f['Header'].attrs['NumPart_ThisFile'][p]
        for jc in range(numparts):
            idsic[ip].append(idsict[jc])
            posic[ip].append(posict[jc])

posiccommon = []
for ip in range(ntypes):
    common = nonzero(in1d(idsic[ip],ids[ip]))[0] #indices of ids in idsic
    posiccommon.append(array([posic[ip][j] for j in common]))  #positions of selected particles
nptot       = sum(len(posiccommon[i]) for i in range(ntypes))
posiccommon = reshape(morph.flatten(posiccommon),(nptot,3))
posiccommon = posiccommon/box

#Write
fname = snappath+"/partpos_ic_center_x_"+str(round(xc/1000.,2))+"_y_"+str(round(yc/1000.,2))+"_z_"+str(round(zc/1000.,2))+"_"+str(size/1000.)+"_hMpc.txt"
savetxt(fname, posiccommon,fmt='%s')#,delimiter="\t")

print("Wrote region to "+fname)
