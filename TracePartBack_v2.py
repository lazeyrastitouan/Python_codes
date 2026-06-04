from numpy import*
import sys
import os
import h5py
import argparse
from multiprocessing import Pool
import morph

# What to do
parser = argparse.ArgumentParser(description='Trace particles of a cubic region of a simulation centered on a given particle back to initial condition. Args: < path to snapshot (hdf5) (str) >, < initial conditions file (str) >, < xc, yc, zc (floats) >, < new box size (float) >, < shift_x, shift_y, shift_z (floats) > < ID most bound (int64, optional) >')
parser.add_argument('snappath',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('icspath',nargs=1,type=str, help='Path to initial conditions (hdf5)')
parser.add_argument('xyz',nargs=3,type=float, help='Center of cube xc, yc, zc [kpc/h]')
parser.add_argument('size',nargs=1,type=float, help='New box size [kpc/h]')
parser.add_argument('--shifts',nargs=3,type=float, help='In case the region is determined from a zoom-in, the shift in kpc/h that was applied in MUSIC (will be applied backwards)')
parser.add_argument('--ID',nargs=1,type=int64, help='ID of most bound particle (optional)')
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
safety   = size/100.
ptypes   = [1,2]
ntypes   = len(ptypes)
verbose  = False
shiftx   = 0.0
shifty   = 0.0
shiftz   = 0.0
if args.shifts is not None:
    shiftx = args.shifts[0]
    shifty = args.shifts[1]
    shiftz = args.shifts[2]
IDMB     = 0
if args.ID is not None:
    IDMB     = args.ID[0]

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
        idst = array(f['PartType'+str(p)]['ParticleIDs'])
        pos  = array(f['PartType'+str(p)]['Coordinates'])

        numparts = f['Header'].attrs['NumPart_ThisFile'][p]	
        for j in range(numparts):
            if xmin <= pos[j,0] <= xmax and ymin <= pos[j,1] <= ymax and zmin <= pos[j,2] <= zmax:		
                ids[ip].append(idst[j])
    
print("Done. Found "+str(len(ids[0])+len(ids[1]))+" particles")

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
        idsict = array(f['PartType'+str(p)]['ParticleIDs'])
        posict = array(f['PartType'+str(p)]['Coordinates'])

        numparts = f['Header'].attrs['NumPart_ThisFile'][p]
        for jc in range(numparts):
            idsic[ip].append(idsict[jc])
            posic[ip].append(posict[jc])

posk  = []
for ip in range(ntypes):
    common = nonzero(in1d(idsic[ip],ids[ip]))[0] #indices of ids in idsic
    posk.append(array([posic[ip][j] for j in common]))  #positions of selected particles
nptot = sum(len(posk[i]) for i in range(ntypes))
posk  = reshape(morph.flatten(posk),(nptot,3))
numk = len(posk)

if IDMB > 0:
    #find central particle in ICs
    IDMBic = idsic.index(IDMB)
    posMB  = posic[IDMBic]-array([shiftx,shifty,shiftz])

#3 Make region cubic and add safety boundary
#Find smallest norm and largest norm to get corner and size of cube
norm     = linalg.norm(posk,axis=1)
normmin  = argmin(norm)
normmax  = argmax(norm)
llc      = posk[normmin] #lower left corner
urc      = posk[normmax] #upper right corner

#second option
xmin2 = min(posk[:,0])
ymin2 = min(posk[:,1])
zmin2 = min(posk[:,2])
xmax2 = max(posk[:,0])
ymax2 = max(posk[:,1])
zmax2 = max(posk[:,2])
llc2  = array([xmin2,ymin2,zmin2])
urc2  = array([xmax2,ymax2,zmax2])

if verbose:
	print("Lower left corner before safety region: "+str(llc/1000.)+"Mpc/h")
	print("Lower left corner second way before safety region:"+str(llc2/1000.)+"Mpc/h")
	print("Upper right corner before safety region: "+str(urc/1000.)+"Mpc/h")
	print("Upper right corner second way before safety region:"+str(urc2/1000.)+"Mpc/h")
llc      -= array([safety,safety,safety])
urc      += array([safety,safety,safety])
llc2     -= array([safety,safety,safety])
urc2     += array([safety,safety,safety])
if verbose:
	print("Lower left corner after safety region: "+str(llc/1000.)+"Mpc/h")
	print("Lower left corner second way after safety region: "+str(llc2/1000.)+"Mpc/h")
	print("Upper right corner after safety region: "+str(urc/1000.)+"Mpc/h")
	print("Upper right corner second way after safety region: "+str(urc2/1000.)+"Mpc/h")
	
x      = min(xmin-safety,llc[0]) 
y      = min(ymin-safety,llc[1]) 
z      = min(zmin-safety,llc[2])
xup    = max(xmax+safety,urc[0]) 
yup    = max(ymax+safety,urc[1]) 
zup    = max(zmax+safety,urc[2])

#apply shift back
llc  -= array([shiftx,shifty,shiftz])
urc  -= array([shiftx,shifty,shiftz])
llc2 -= array([shiftx,shifty,shiftz])
urc2 -= array([shiftx,shifty,shiftz])
if verbose:
    print("Lower left corner after shift region back: "+str(llc/1000.)+"Mpc/h")
    print("Lower left corner second way after shift region back: "+str(llc2/1000.)+"Mpc/h")
    print("Upper right corner after shift region back: "+str(urc/1000.)+"Mpc/h")
    print("Upper right corner second way after shidt region back: "+str(urc2/1000.)+"Mpc/h")


x    -= shiftx
y    -= shifty
z    -= shiftz
xup  -= shiftx
yup  -= shifty
zup  -= shiftz

sizeic = max(xup-x,yup-y,zup-z)
sizeic2 = max(xmax2-xmin2,ymax2-ymin2,zmax2-zmin2)

print("Lower left corner at ["+str(round(x/1000.,3))+", "+str(round(y/1000.,3))+", "+str(round(z/1000.,3))+"] [Mpc/h]")
print("Lower left corner second way at ["+str(round(llc2[0]/1000.,3))+", "+str(round(llc2[1]/1000.,3))+", "+str(round(llc2[2]/1000.,3))+"] [Mpc/h]")
print("1D Size of region in ICs: "+str(round(sizeic/1000.,3))+" Mpc/h")
print("1D Size of region in ICs second way: "+str(round(sizeic2/1000.,3))+" Mpc/h")
print("Center of region at ["+str(round((x+sizeic/2.)/1000.,3))+", "+str(round((y+sizeic/2.)/1000.,3))+", "+str(round((z+sizeic/2.)/1000.,3))+"] [Mpc/h]")
print("Center of region second way at ["+str(round((llc2[0]+sizeic2/2.)/1000.,3))+", "+str(round((llc2[1]+sizeic2/2.)/1000.,3))+", "+str(round((llc2[2]+sizeic2/2.)/1000.,3))+"] [Mpc/h]")
if IDMB>0:
    print("Center particle at z=3 is now at ["+str(round(posMB[0]/1000.,3))+", "+str(round(posMB[1]/1000.,3))+", "+str(round(posMB[2]/1000.,3))+"] [Mpc/h]")

fname = snappath+"/ic_region_center_x_"+str(round((xc-shiftx)/1000.,2))+"_y_"+str(round((yc-shifty)/1000.,2))+"_z_"+str(round((zc-shiftz)/1000.,2))+"_"+str(size/1000.)+"_Mpc.txt"
f     = open(fname,"w")
f.write("Center from norm and half size\n")
f.write(str(x+sizeic/2.)+"\t"+str(y+sizeic/2.)+"\t"+str(z+sizeic/2.)+"\t"+str(sizeic/2.)+"\n")
f.write("Lower left corner from norm and size\n")
f.write(str(x)+"\t"+str(y)+"\t"+str(z)+"\t"+str(sizeic)+"\n")
f.write("Box units\n")
f.write(str((x+sizeic/2.)/box)+"\t"+str((y+sizeic/2.)/box)+"\t"+str((z+sizeic/2.)/box)+"\t"+str(sizeic/(2.*box))+"\n") # box units needed for music
f.write(str(x/box)+"\t"+str(y/box)+"\t"+str(z/box)+"\t"+str(sizeic/box)+"\n") # box units needed for music
f.write("Center from min of coords and half size\n")
f.write(str(llc2[0]+sizeic2/2.)+"\t"+str(llc2[1]+sizeic2/2.)+"\t"+str(llc2[2]+sizeic2/2.)+"\t"+str(sizeic2/2.)+"\n")
f.write("Lower left corner from min ofcoords and size\n")
f.write(str(llc2[0])+"\t"+str(llc2[1])+"\t"+str(llc2[2])+"\t"+str(sizeic2)+"\n")
f.write("Box units\n")
f.write(str((llc2[0]+sizeic2/2.)/box)+"\t"+str((llc2[1]+sizeic2/2.)/box)+"\t"+str((llc2[2]+sizeic2/2.)/box)+"\t"+str(sizeic2/(2.*box))+"\n") # box units needed for music
f.write(str(llc2[0]/box)+"\t"+str(llc2[1]/box)+"\t"+str(llc2[2]/box)+"\t"+str(sizeic2/box)+"\n") # box units needed for music
if IDMB>0:
    f.write("Central particle in kpc/h and box units\n")
    f.write(str(posMB[0])+"\t"+str(posMB[1])+"\t"+str(posMB[2])+"\n")
    f.write(str(posMB[0]/box)+"\t"+str(posMB[1]/box)+"\t"+str(posMB[2]/box)+"\n")
f.close()

print("Wrote region to "+fname)
