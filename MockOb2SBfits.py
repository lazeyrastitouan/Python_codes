from numpy import*
import argparse
from astropy.io import fits

# input parameters
parser = argparse.ArgumentParser(description='Routine to transform MockOb.fits to DisPersE/SBprofile.py readable intrinsic SB grids fits. Args: < MockOb file (str) >, < redshift (float) >, < noise (str) >, < sigma smoothing > (int) > < --cubex (optional) >, < --nonoise_file (optional) >')
parser.add_argument('mockob',nargs=1,type=str, help='MockOb .fits file')
parser.add_argument('redshift',nargs=1,type=float,help='Redshift')
parser.add_argument('sigma',nargs=1,type=float,help='Sigma for the smoothing')
parser.add_argument('--cubex',action='store_true',help="Optional, indicates to use the cubex .Objects_Id.fits file to keep only 1 nebula and nan otherwise (str). Must be accompanied by at .cat file with same name")
parser.add_argument('--nonoise_file',type=str,nargs=1,help='Optional, if you want to get the 2D image of the mockCube w/o noise based on the nebulae extracted on noisy data')
args  = parser.parse_args()
f1    = args.mockob[0]
red   = args.redshift[0]
s     = args.sigma[0]
if args.cubex:
    fcube = f1[:-5]+'.Objects_Id.fits'
    dhuc  = fits.open(fcube)
    cube  = dhuc[0].data
    cen   = 128.

#some hard-coded settings
specres   = 1.25 #Angstrom
verbose   = True
singleobj = False #True #wether to keep only 1 object from cubex or all of them
NB        = False #True #pseudo narrow-band images instead of optimally extracted

print("Reading data...")
print(f1)
hdul = fits.open(f1)
data = hdul[0].data #z direction first

if args.cubex and singleobj:
    #find ID of most extended object
    ID = 1
    xy = array([cen,cen])
    with open(f1[:-5]+'.cat') as cat:
        lines = cat.readlines()
    cat.close()
    l = 0
    while lines[l][0] == '#':
        l+=1
    xyn  = [float(lines[l].split()[2]),float(lines[l].split()[3])]
    cell = int(lines[l].split()[1])
    l_n  = len(lines)
    if l_n > l+1:
        IDS   = [int(lines[l+i].split()[0]) for i in range(l_n-l)]
        cells = [int(lines[l+i].split()[1]) for i in range(l_n-l)]
        xyneb = [[float(lines[l+i].split()[2]),float(lines[l+i].split()[3])] for i in range(l_n-l)]
        xyneb = array((xyneb))
        dist  = linalg.norm(xy-xyneb,axis=1)
        
        #first let's take most extended object as one associated to halo
        obj   = cells.index(max(cells))
        ID    = IDS[obj]
        xyn   = xyneb[obj]
        cell  = cells[obj]
        #check that halo center is in nebula
        if ID not in cube[:,int(xy[0]-1):int(xy[0]+1),int(xy[1]-1):int(xy[1]+1)]:
            #if not, let's take the nebula with the closet "center" to the halo center
            obj   = dist.tolist().index(min(dist))
            ID    = IDS[obj]
            xyn   = xyneb[obj]
            cell  = cells[obj]
            #check that halo center is in nebula
            if ID not in cube[:,int(xy[0]-1):int(xy[0]+1),int(xy[1]-1):int(xy[1]+1)]:
                exit("Neither the most extended object nor the closest one cover halo center, disgarding nebula")

    print("ID of central objec: "+str(ID))
    print("Number of pixels: "+str(cell))
    print("Position: "+str(xyn)+" [pixels]")

    #get nebula (ev. in narrow band), 0 otherwise (so won't count in sum)
    if NB: 
        minz = min(where(cube==ID)[0]) #where(cube==ID)[0] is the z coordinate of nebula
        maxz = max(where(cube==ID)[0])
        print("Narrow band image between z-layers "+str(minz)+" and "+str(maxz))
        if args.nonoise_file is not None:
            print("Reading no noise data")
            hduln0 = fits.open(args.nonoise_file[0])
            data = hduln0[0].data #replace the data by the non noisy one
        data[0:minz]=0.
        data[maxz+1:shape(data)[0]-1]=0.
    else:
        data = where(cube==ID,data,0)

elif args.cubex and not singleobj:
    data = where(cube!=0,data,0)

#get 2D nebula
dataplot  = data.sum(axis=0)*specres*1e-20 #*(1.+red)**4. specres is the size of voxel in z direction so sum *specres is the integral
grid      = shape(dataplot)[0]

if verbose:
	print("Image size in pixel "+str(shape(dataplot)))
	print("maximum SB: {}".format(amax(dataplot)))
	print("minimum SB: {}".format(amin(dataplot)))

#replace zeros by nan to not be taken into account in profile when doing nanmedian
dataplot = where(dataplot==0,nan,dataplot)

d  = fits.PrimaryHDU(dataplot)
if args.cubex and singleobj and NB:
    if args.nonoise_file is not None:
        fo = f1[:-5]+'_SB_s'+str(s)+'_z'+str(red)+'_Ng'+str(grid)+'_observed_cubex_neb_singleobj_NB_nonoise.fits'
    else:
        fo = f1[:-5]+'_SB_s'+str(s)+'_z'+str(red)+'_Ng'+str(grid)+'_observed_cubex_neb_singleobj_NB.fits'
elif args.cubex and singleobj and not NB:
    fo = f1[:-5]+'_SB_s'+str(s)+'_z'+str(red)+'_Ng'+str(grid)+'_observed_cubex_neb_singleobj.fits'
elif args.cubex and not singleobj:
    fo = f1[:-5]+'_SB_s'+str(s)+'_z'+str(red)+'_Ng'+str(grid)+'_observed_cubex_neb.fits'
else:
    fo = f1[:-5]+'_SB_s'+str(s)+'_z'+str(red)+'_Ng'+str(grid)+'_observed.fits'
print("Writing to "+fo)
d.writeto(fo,overwrite=True)

print("Done")
