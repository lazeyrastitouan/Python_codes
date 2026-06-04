from numpy import*
import sys
import os
import h5py
import argparse
from multiprocessing import Pool
import morph

# What to do
parser = argparse.ArgumentParser(description='Find closest halos to an input position in AHF output and display line number, M [Msun/h], x, y, z [ckpc/h], vx, vy, vz ckm/s, radius [ckpc/h]. Input: < x,y,z [ckpc/h] (3floats) >, < Path to AHF catalog (str) >')
parser.add_argument('xyz',nargs=3,type=float, help='Input position x, y, z [ckpc/h]')
parser.add_argument('path2cat',nargs=1,type=str, help='Absoulte path to AHF catalog')
args     = parser.parse_args()
path2cat = args.path2cat[0]
xc 	     = args.xyz[0]
yc 	     = args.xyz[1]
zc 	     = args.xyz[2]
xyz      = array([xc,yc,zc])
verbose  = False

cat      = loadtxt(path2cat,usecols=(3,5,6,7,8,9,10,11)).T
xyzhalo  = array([cat[1],cat[2],cat[3]]).T #loadtxt(path2cat,usecols=(5,6,7))

diff     = xyzhalo-xyz
norm     = linalg.norm(diff,axis=1)
normmin  = min(norm)
idmin    = norm.tolist().index(normmin)

path     = path2cat.rsplit('/', 1)[0]
outname  = 'closest_halo_'+str(int(xc))+'_'+str(int(yc))+'_'+str(int(zc))+'.txt'
outfile  = open(path+'/'+outname,'w')
outfile.write(str(int(xc))+"\t"+str(int(yc))+"\t"+str(zc)+"\t"+str(normmin)+"\n")
#write cat in correct order for genP2Cgrid input: x,y,z,R,M,vx,vy,vz
outfile.write(str(cat[1,idmin])+"\t"+str(cat[2,idmin])+"\t"+str(cat[3,idmin])+"\t"+str(cat[7,idmin])+"\t"+str(cat[0,idmin])+"\t"+str(cat[4,idmin])+"\t"+str(cat[5,idmin])+"\t"+str(cat[6,idmin])+"\n")
outfile.close()
print("Wrote output to "+path+'/'+outname)
