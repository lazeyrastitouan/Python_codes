from numpy import*
import sys
import os
import h5py
import argparse
from multiprocessing import Pool

# What to do
parser = argparse.ArgumentParser(description='Trace particles of a cubic region of a simulation centered on a given position (an corresponding lower left center) back to initial condition. Args: < path to snapshot (hdf5) (str) >, < initial conditions file (str) >, < xc, yc, zc (floats) >, < new box size (float) >')
parser.add_argument('snappath',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('icspath',nargs=1,type=str, help='Path to initial conditions (hdf5)')
parser.add_argument('xyz',nargs=3,type=float, help='Center of cube xc, yc, zc [kpc/h]')
parser.add_argument('size',nargs=1,type=float, help='New box size [kpc/h]')
args     = parser.parse_args()
snappath = args.snappath[0]
icspath  = args.icspath[0]
xc 	 = args.xyz[0]
yc 	 = args.xyz[1]
zc 	 = args.xyz[2]
size     = args.size[0]
xmin     = xc-size/2.
ymin     = yc-size/2.
zmin     = zc-size/2.
xmax     = xc+size/2.
ymax     = yc+size/2.
zmax     = zc+size/2.
safety   = size/25.
verbose  = False

h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5']
h5files  = sorted(h5files)
h5files  = sorted(h5files,key=len)
if verbose:
	print(h5files)

print("Keeping particule in cube\n"+str(xmin/1000.)+" < x < "+str(xmax/1000.)+" Mpc/h\n"+str(ymin/1000.)+" < y < "+str(ymax/1000.)+" Mpc/h\n"+str(zmin/1000.)+" < z < "+str(zmax/1000.)+" Mpc/h\n")
print("Size of safety region: "+str(safety/1000.)+" Mpc/h")
 
#1 Find particle ids in a given snapshot region
print("Selecting particles in volume (only considering highest refinment ones)...")
ids=[]
for i in h5files:
	f    = h5py.File(snappath+"/"+i, "r")
	idst = array(f['PartType1']['ParticleIDs'])
	pos  = array(f['PartType1']['Coordinates'])
	
	numparts = f['Header'].attrs['NumPart_ThisFile'][1]	
	for j in range(numparts):
		if xmin <= pos[j,0] <= xmax and ymin <= pos[j,1] <= ymax and zmin <= pos[j,2] <= zmax:		
			ids.append(idst[j])

print("Done. Found "+str(len(ids))+" particles")

#2 Find the corresponding volume occupied by these particles in the ICs
print("Finding corresponding region in ICs...")

h5ics        = [f for f in os.listdir(icspath+"/") if os.path.isfile(os.path.join(icspath+"/", f)) and f[-4:] == 'hdf5']
h5ics        = sorted(h5ics)
h5ics        = sorted(h5ics,key=len)
idsic, posic = [],[]
box          = h5py.File(icspath+'/'+h5ics[0],'r')['Header'].attrs['BoxSize']
for ic in h5ics:
	f      = h5py.File(icspath+"/"+ic,'r')
	idsict = array(f['PartType1']['ParticleIDs'])
	posict = array(f['PartType1']['Coordinates'])
	
	numparts = f['Header'].attrs['NumPart_ThisFile'][1]
	for jc in range(numparts):
		idsic.append(idsict[jc])
		posic.append(posict[jc])

common = nonzero(in1d(idsic,ids))[0] #indices of ids in idsic
posk   = array([posic[j] for j in common])  #positions of selected particles
numk   = len(posk)

#3 Make region cubic and add safety boundary
#Find smallest norm and largest norm to get corner and size of cube
norm     = linalg.norm(posk,axis=1)
normmin  = argmin(norm)
normmax  = argmax(norm)
llc      = posk[normmin] #lower left corner
urc      = posk[normmax] #upper right corner
if verbose:
	print("Lower left corner before safety region: "+str(llc/1000.)+"Mpc/h")
	print("Upper right corner before safety region: "+str(urc/1000.)+"Mpc/h")
llc      -= array([safety,safety,safety])
urc      += array([safety,safety,safety])
if verbose:
	print("Lower left corner after safety region: "+str(llc/1000.)+"Mpc/h")
	print("Upper right corner after safety region: "+str(urc/1000.)+"Mpc/h")
	
x      = min(xmin-safety,llc[0]) 
y      = min(ymin-safety,llc[1]) 
z      = min(zmin-safety,llc[2])
xup    = max(xmax+safety,urc[0]) 
yup    = max(ymax+safety,urc[1]) 
zup    = max(zmax+safety,urc[2])
sizeic = max(xup-x,yup-y,zup-z)

print("Lower left corner at ["+str(round(x/1000.,3))+", "+str(round(y/1000.,3))+", "+str(round(z/1000.,3))+"] [Mpc/h]")
print("1D Size of region in ICs: "+str(round(sizeic/1000.,3))+" Mpc/h")
print("Center of region at ["+str(round((x+sizeic/2.)/1000.,3))+", "+str(round((y+sizeic/2.)/1000.,3))+", "+str(round((z+sizeic/2.)/1000.,3))+"] [Mpc/h]")

fname = snappath+"/ic_region_center_x_"+str(round(xc/1000.,2))+"_y_"+str(round(yc/1000.,2))+"_z_"+str(round(zc/1000.,2))+"_"+str(size/1000.)+"_Mpc.txt"
f     = open(fname,"w")
f.write(str(x+sizeic/2.)+"\t"+str(y+sizeic/2.)+"\t"+str(z+sizeic/2.)+"\t"+str(sizeic/2.)+"\n")
f.write(str(x)+"\t"+str(y)+"\t"+str(z)+"\t"+str(sizeic)+"\n")
f.write(str((x+sizeic/2.)/box)+"\t"+str((y+sizeic/2.)/box)+"\t"+str((z+sizeic/2.)/box)+"\t"+str(sizeic/(2.*box))+"\n") # box units needed for music
f.write(str(x/box)+"\t"+str(y/box)+"\t"+str(z/box)+"\t"+str(sizeic/box)+"\n") # box units needed for music
f.close()

print("Wrote region to "+fname)
