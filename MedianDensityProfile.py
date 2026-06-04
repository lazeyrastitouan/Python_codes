from numpy import*
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import argparse
import sys, os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

# input parameters
parser  = argparse.ArgumentParser(description='Compute median SBprofile of all grids in a mass bin for Halo10 and Halo1312 (stack 3 LOS) Args: < log10(Mlow) (float) >, < dark matter (str) >, < redshift (float) >, < zoom (str) >, <--Tmax (str) >, < --dM (optional, default = 0.5, float) >')
parser.add_argument('mlow',nargs=1,type=float, help='Lower edge of the mass bin in log10 (float)')
parser.add_argument('DM',nargs=1,type=str, help='Dark matter type (str)')
parser.add_argument('red',nargs=1,type=float, help='Redshift (float)')
parser.add_argument('zoom',nargs=1,type=str, help='Zoom (1p5 or 2) (str)')
parser.add_argument('tmax',nargs=1,type=str, help='Max temperature Xe+Y format (str)')
parser.add_argument('--deltam',nargs=1,type=float, help='Mass bin width in log10 (float)')
args   = parser.parse_args()
mlow   = args.mlow[0]
DM     = args.DM[0]
red    = args.red[0]
zm     = args.zoom[0]
Tmax   = args.tmax[0]
deltaM = 0.5
if args.deltam is not None:
    deltaM = args.deltam[0]
mhi    = mlow+deltaM

# ---

#LOAD
if red == 4.0:
    halos = ['10','1312']
    snaps = ['1','0']
elif red == 3.5:
    halos = ['10','1312']
    snaps = ['4','1']
elif red == 3.0:
    if mlow == 12.0:
        halos = ['10']
        snaps = ['7']
    elif mlow == 13.0:
        halos = ['1312']
        snaps = ['2']
    else:
        halos = ['10','1312']
        snaps = ['7','2']

root   = '/data/titouan/Eagle_cosmo_gizmo/200_512/'
zoom   = 'zoom'+zm+'_hydro_fire2_lvl7'
hydro  = 'sftresh_cool_low_fb_boost'
ns     = '' #'nstar1/'

ODtab, Dtab  = [],[]
for i in range(len(halos)):
    path   = root+'/'+DM+'/'+zoom+'/Halo'+halos[i]+'/'+hydro+'/snapdir_00'+snaps[i]+'/2cMpc/'+ns
    print("Loading grids in "+path)
    gridsf = root+'/'+DM+'/'+zoom+'/Halo'+halos[i]+'/'+hydro+'/AHF_00'+snaps[i]+'_allparts/200_snap_z'+str(red)+'00_halos_sorted_main_hires_Mlow_'+str(mlow)+'_Mhigh_'+str(mhi)+'.txt'
    grids  = genfromtxt(gridsf,names=True,delimiter=';')['gn']
    if shape(grids) == ():
        grids = [grids]
    for g in grids:
        g = int(g)
        f = 'DensityProfile_Tmin_1e+0_Tmax_'+Tmax+'_P2C_grid_'+str(g)+'_256.npz'
        data = load(path+'/'+f)
        ODtab.append(data['overdensity'])
        Dtab.append(data['density'])

print("Loaded "+str(len(ODtab))+" profiles")

#R is same for all, use last SB profile path
R = load(path+'/'+f)['binmean']

# ---
#COMPUTE MEDIAN OF MEDIANS AND PRECENTILES 
#fout = root+"/"+hydro+ns+"Densityprofile_"+zoom+"_"+DM+"_z_"+str(red)+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhi)+".npz"
fout = root+"/"+ns+"Densityprofile_"+zoom+"_"+DM+"_z_"+str(red)+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhi)+"_Tmax_"+Tmax+".npz"
print("Writing profile to "+fout)
#if only 1 halo, keep 3 los
if len(grids) == 1:
    savez(fout,binmean=R,odprofileX=ODtab[0],odprofileY=ODtab[1],odprofileZ=ODtab[2],dprofileX=Dtab[0],dprofileY=Dtab[1],dprofileZ=Dtab[2])
else:
    OD   = nanmedian(ODtab,axis=0)
    OD1  = nanpercentile(ODtab,0.001,axis=0)
    OD2  = nanpercentile(ODtab,0.01,axis=0)
    OD3  = nanpercentile(ODtab,0.1,axis=0)
    OD4  = nanpercentile(ODtab,1,axis=0)
    OD5  = nanpercentile(ODtab,10,axis=0)
    OD6  = nanpercentile(ODtab,25,axis=0)
    OD7  = nanpercentile(ODtab,75,axis=0)
    OD8  = nanpercentile(ODtab,90,axis=0)
    OD9  = nanpercentile(ODtab,99,axis=0)
    OD10 = nanpercentile(ODtab,99.9,axis=0)
    OD11 = nanpercentile(ODtab,99.99,axis=0)
    OD12 = nanpercentile(ODtab,99.999,axis=0)

    D   = nanmedian(Dtab,axis=0)
    D1  = nanpercentile(Dtab,0.001,axis=0)
    D2  = nanpercentile(Dtab,0.01,axis=0)
    D3  = nanpercentile(Dtab,0.1,axis=0)
    D4  = nanpercentile(Dtab,1,axis=0)
    D5  = nanpercentile(Dtab,10,axis=0)
    D6  = nanpercentile(Dtab,25,axis=0)
    D7  = nanpercentile(Dtab,75,axis=0)
    D8  = nanpercentile(Dtab,90,axis=0)
    D9  = nanpercentile(Dtab,99,axis=0)
    D10 = nanpercentile(Dtab,99.9,axis=0)
    D11 = nanpercentile(Dtab,99.99,axis=0)
    D12 = nanpercentile(Dtab,99.999,axis=0)

    #WRITE
    savez(fout,binmean=R,od0001=OD1,od001=OD2,od01=OD3,od1=OD4,od10=OD5,od25=OD6,od50=OD,od75=OD7,od90=OD8,od99=OD9,od999=OD10,od9999=OD11,od99999=OD12,d0001=D1,d001=D2,d01=D3,d1=D4,d10=D5,d25=D6,d50=D,d75=D7,d90=D8,d99=D9,d999=D10,d9999=D11,d99999=D12)
