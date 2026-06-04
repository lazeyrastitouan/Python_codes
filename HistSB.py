from numpy import*
import h5py, sys, os
import argparse
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import eagleTools as eT
from decimal import Decimal 
import glob

# input parameters
parser = argparse.ArgumentParser(description='Routine get pixel histogram of CubEx projected 2D nebulae. Args: < simulation (str) >, < redhsift (float) >, < grid list file (str) >, < --nonoise (optional) > ')
parser.add_argument('sim',nargs=1,type=str, help='simulation type (dalya/colibre/hello)')
parser.add_argument('red',nargs=1,type=float, help='redshift')
parser.add_argument('grids',nargs=1,type=str, help='file with list of grids')
parser.add_argument('--nonoise',action='store_true',help="Optional, substracts the histogram of pure noise")
args   = parser.parse_args()
sim    = args.sim[0]
red    = args.red[0]
gridsf = args.grids[0]

#some hard-coded settings
PDF    = False
Lbox   = 2000. #4000. #ckpc
boxint = '2' #'4'
grid   = '256' #'512'
edges  = [40.,100.,200.,300.]
s      = '1.5'
n      = '5e-20'
suf    = ''#'noISMcut_' #'dust-7_' '' 'dust-5_' 'metal005_' 
print("Noise = "+n)
if PDF:
    print("Histogram of PDF")
else:
    print("Histogram of number of pixels")

#get files
files = []
if sim == 'dalya':
    path ='/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/'
    H     = ['1312','10']
    if red == 4.0:
        snaps = ['0','1'] 
    elif red == 3.5:
        snaps = ['1','4']
    elif red == 3.0:
        snaps = ['2','7']
    for si in range(len(snaps)):
        grids = genfromtxt(path+'Halo'+H[si]+'/sftresh_cool_low_fb_boost/AHF_00'+snaps[si]+'_allparts/'+gridsf,names=True,delimiter=';')['gn']
        if not isinstance(grids, list): #need this for empty lists
            if isinstance(grids,float64):
                grids = [float(grids)]
        print("Loading grids "+str(grids))
        grids = [float(grids)]
        for g in grids:
            newfiles = glob.glob(path+'Halo'+H[si]+'/sftresh_cool_low_fb_boost/snapdir_00'+snaps[si]+'/'+boxint+'cMpc/'+suf[:-1]+'/mockCube*noise'+n+'*_gn'+str(int(g))+'_SB_s'+s+'_z'+str(red)+'_Ng'+grid+'_observed_cubex_neb_singleobj_NB_nonoise.fits') 
            files    = files + newfiles
elif sim == 'colibre':
    colsim = 'L0050N1504'
    if red == 4.0:
        snap = '56'
        path = '/data/titouan/COLIBRE/'+colsim+'/z4/'
    elif red == 3.5:
        snap = '60'
        path = '/data/titouan/COLIBRE/'+colsim+'/z3p5/'
    elif red == 3.0:
        snap = '64'
        path = '/data/titouan/COLIBRE/'+colsim+'/z3/'
    grids = genfromtxt(path+gridsf,names=True,delimiter=';')['gn']
    print("Loading grids "+str(grids))
    for g in grids:
        newfiles = glob.glob(path+'colibre_00'+snap+'/'+boxint+'cMpc/'+suf[:-1]+'/mockCube*noise'+n+'*_gn'+str(int(g))+'_SB_s'+s+'_z'+str(red)+'_Ng'+grid+'_observed_cubex_neb_singleobj_NB_nonoise.fits')
        files    = files + newfiles
    #new path for outfile
    path = '/data/titouan/COLIBRE/'+colsim+'/'
elif sim == 'hello':
    if red == 3.0:
        exit("HELLO has no redshift 3 output!")
    path = '/data/titouan/HELLO/CDM/'
    snaps = ['g9.20e12','g2.71e12','g2.32e12']
    for si in snaps:
        grids = genfromtxt(path+'/'+si+'/'+si+'.'+gridsf,names=True,delimiter=';')['gn']
        if not isinstance(grids, list):
            grids = atleast_1d(grids).tolist()
        print("Loading grids "+str(grids))
        for g in grids:
            newfiles = glob.glob(path+'/'+si+'/'+boxint+'cMpc/'+suf[:-1]+'/mockCube*noise'+n+'*_gn'+str(int(g))+'_SB_s'+s+'_z'+str(red)+'*_Ng'+grid+'_observed_cubex_neb_singleobj_NB_nonoise.fits')
            files    = files + newfiles

print('Files '+str(files))
print("Loading "+str(len(files))+" mockCubes")

#get SBs
SBs   = [[] for x in range(len(edges)-1)] #final shape is (number of pixels, len(edges)-1)
for f in files:
    dhu       = fits.open(f)
    data      = dhu[0].data
    sideL     = shape(data)[0]
    center    = sideL/2
    gridTokpc = Lbox/sideL
    radDistGrid = eT.calc2DRadialDistance(sideL, sideL, center, center)*gridTokpc
    for i in range(len(edges)-1):
        Rfilter = logical_and(edges[i] <= radDistGrid, radDistGrid < edges[i + 1])
        SBs[i]  = SBs[i] + data[Rfilter].tolist()


print(log10(mean(SBs[1])*(1+red)**4.))


SBn = [[] for x in range(len(edges)-1)]
if args.nonoise:
    print("Removing noise from histogram")
    filesn = glob.glob('/data/titouan/pure_noise/mockCube*noise'+n+'_z'+str(red)+'_Ng_256_SB_s'+s+'*observed.fits')
    suf += 'wo_noise_'
    for f in filesn:
        dhu       = fits.open(f)
        data      = dhu[0].data
        radDistGrid = eT.calc2DRadialDistance(sideL, sideL, center, center)*gridTokpc
        for i in range(len(edges)-1):
            Rfilter = logical_and(edges[i] <= radDistGrid, radDistGrid < edges[i + 1])
            SBn[i]  = SBn[i] + data[Rfilter].tolist()

#histogram
if not PDF:
    suf = 'Npixel_'+suf
for i in range(len(edges)-1):
    hist, sb_edges = histogram(log10(array(SBs[i])),range=(-23,-14),density=PDF,bins=70) #abs for negative noise
    if i == 1:
        low = (1+red)**4.*average(10**sb_edges[:-1], weights=hist)
        high = (1+red)**4.*average(10**sb_edges[1:], weights=hist)
        print('The average is in the range '+str(low)+' - '+str(high))
    if args.nonoise:
        histn, sb_edgesn = histogram(log10(abs(array(SBn[i]))),range=(-21,-15),density=PDF,bins=60) #abs for negative noise
        #if i == 1:
        #    plt.stairs(hist, 10**sb_edges, fill=False,label='w/ noise')
        #    plt.stairs(histn, 10**sb_edgesn, fill=False,label='noise')
        hist = hist-histn
        #if i == 1:
        #    plt.stairs(hist, 10**sb_edges, fill=False,label='w/o noise')
        #    plt.legend(loc=0,frameon=False)
        #    plt.show()
    fout = path+'/Hist_SB_'+suf+'n'+n+'_s'+s+'_Ng'+str(int(sideL))+'_observed_cubex_neb_singleobj_NB_nonoise_z'+str(red)+'_grids'+str(int(grids[0]))+'-'+str(int(grids[len(grids)-1]))+'_R'+str(int(edges[i]))+'-'+str(int(edges[i+1]))+'.npz'
    savez(fout,hist=hist,edges=sb_edges)
    print("Saved histogram to "+fout)
