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
parser  = argparse.ArgumentParser(description='Compute median SBprofile of all grids in a mass bin for Halo10 and Halo1312 (stack 3 LOS) Args: < log10(Mlow) (float >, < dark matter (str) >, < redshift (float) >, < zoom (str) >, < noise (str) >, < --dM (optional, default = 0.5, float) >')
parser.add_argument('mlow',nargs=1,type=float, help='Lower edge of the mass bin in log10 (float)')
parser.add_argument('DM',nargs=1,type=str, help='Dark matter type (str)')
parser.add_argument('red',nargs=1,type=float, help='Redshift (float)')
parser.add_argument('zoom',nargs=1,type=str, help='Zoom (1p5 or 2) (str)')
parser.add_argument('noise',nargs=1,type=str,help='MockOb noise (str)')
parser.add_argument('--deltam',nargs=1,type=float, help='Mass bin width in log10 (float)')
args   = parser.parse_args()
mlow   = args.mlow[0]
DM     = args.DM[0]
red    = args.red[0]
zm     = args.zoom[0]
noise  = args.noise[0]
deltaM = 0.5
if args.deltam is not None:
    deltaM = args.deltam[0]
mhi    = mlow+deltaM

root   = '/data/titouan/Eagle_cosmo_gizmo/200_512/'
zoom   = 'zoom'+zm+'_hydro_fire2_lvl7'
hydro  = 'sftresh_cool_low_fb_boost'
ns     = '' #'metal005/' #'nstar3/' 'nstar100/' 'nstar01/' 'nstar1/'
s_obj  = True

# ---

#LOAD
if red == 4.0:
    if hydro == 'sftresh_cool_low' and mlow == 12.5:
        halos = ['1312']
        snaps = ['0']
    else:
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

SBtab  = []
ng=0
for i in range(len(halos)):
    ng+=1
    path   = root+'/'+DM+'/'+zoom+'/Halo'+halos[i]+'/'+hydro+'/snapdir_00'+snaps[i]+'/2cMpc/'+ns
    print("Loading grids in "+path)
    gridsf = root+'/'+DM+'/'+zoom+'/Halo'+halos[i]+'/'+hydro+'/AHF_00'+snaps[i]+'_allparts/200_snap_z'+str(red)+'00_halos_sorted_main_hires_Mlow_'+str(mlow)+'_Mhigh_'+str(mhi)+'.txt'
    grids  = genfromtxt(gridsf,names=True,delimiter=';')['gn']
    if shape(grids) == ():
        grids = [grids]
    for g in grids:
        g  = int(g)
        if s_obj:
            SBf = [ f for f in os.listdir(path) if os.path.isfile(os.path.join(path,f)) and f[0:8]=='mockCube' and '_noise'+noise+'_gn'+str(g)+'_SB_s1.5_z'+str(red)+'_Ng256_observed_cubex_neb_singleobj_NB_meanSBprofile.npz' in f]
        else:
            SBf = [ f for f in os.listdir(path) if os.path.isfile(os.path.join(path,f)) and f[0:8]=='mockCube' and '_noise'+noise+'_gn'+str(g)+'_SB_s1.5_z'+str(red)+'_Ng256_observed_cubex_neb_SBprofile.npz' in f] 
        #print(SBf)
        if i == 0:
            #R is same for all, use first SB profile path
            R = load(path+'/'+SBf[0])['binmean']
        for f in SBf:
            print(f)
            data = load(path+'/'+f)
            SBtab.append(data['sbprofile'])

print("Loaded "+str(len(SBtab))+" profiles")

# ---
#COMPUTE MEDIAN OF MEDIANS AND PRECENTILES 
if s_obj:
    fout = root+"/"+ns+"SBprofile_"+zoom+"_"+DM+"_observed_cubex_neb_singleobj_NB_z_"+str(red)+"_s1.5_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhi)+".npz"
else:
    fout = root+"/"+ns+"SBprofile_"+zoom+"_"+DM+"_observed_cubex_neb_z_"+str(red)+"_s1.5_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhi)+".npz"
print("Writing SB profile to "+fout)
#if only 1 halo, keep 3 los
if len(grids) == 1 and ng==1:
    #WRITE
    savez(fout,binmean=R,sbprofileX=SBtab[0],sbprofileY=SBtab[1],sbprofileZ=SBtab[2])
else:
    SB  = nanmedian(SBtab,axis=0)
    P16 = nanpercentile(SBtab,16,axis=0)
    P84 = nanpercentile(SBtab,84,axis=0)

    #WRITE
    savez(fout,binmean=R,sbprofile=SB,p16=P16,p84=P84)
