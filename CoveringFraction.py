from numpy import*
import h5py, sys, os
import argparse
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import eagleTools as eT
from decimal import Decimal 

# input parameters
parser = argparse.ArgumentParser(description='Routine get covering fraction of CubEx projected 2D nebulae. Args: < Fits file (str) >, < red (float) >')
parser.add_argument('file',nargs=1,type=str, help='Fits file')
parser.add_argument('red',nargs=1,type=float, help='Redshift')
args     = parser.parse_args()
fitsfile = args.file[0]
red      = args.red[0]

tresh = array([1e-19,5e-19,1e-18,5e-18,1e-17,2e-17,5e-17,7e-17,1e-16,2e-16])*red**4.

#some hard-coded settings
Lbox  = 2000. #ckpc
edges = [100.,200.]#,300.]

dhu       = fits.open(fitsfile)
data      = dhu[0].data*(1.+red)**4.
sideL     = shape(data)[0]
center    = sideL/2
gridTokpc = Lbox/sideL

radDistGrid = eT.calc2DRadialDistance(sideL, sideL, center, center)*gridTokpc

for i in range(len(edges)-1):
    fname      = fitsfile[:-5]+'_CoveringFraction_'+str(int(edges[i]))+'-'+str(int(edges[i+1]))+'.txt'
    fout       = open(fname,'w')
    Rfilter    = logical_and(edges[i] <= radDistGrid, radDistGrid < edges[i + 1])
    rSelection = data[Rfilter]
    if len(rSelection) == 0:
        continue
    else:
        Npixel    = len(rSelection)
        for t in tresh:
            SBfilter = rSelection >= t
            Nselect  = len(rSelection[SBfilter])
            fcover   = Nselect/Npixel
            fout.write(str(t)+"\t"+str(fcover)+"\n")
            del fcover
    fout.close()
    print("Wrote results to "+fname)

