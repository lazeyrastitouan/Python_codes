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
parser  = argparse.ArgumentParser(description='Compute median SBprofile of all grids in a mass bins (stack 3 LOS) Args: < Path .npz SB profiles (str) >, < File with grids in the given mass bin (str) >, < log10(Mlow) (float) >, < --suffix (optional, str) >, < --dM (optional, default = 0.5, float) >')
parser.add_argument('path',nargs=1,type=str, help='Absolute path to folder containing .npz SB profiles (str)')
parser.add_argument('grids',nargs=1,type=str, help='File with grids in the given mass bin (str)')
parser.add_argument('mlow',nargs=1,type=float, help='Lower edge of the mass bin in log10 (float)')
parser.add_argument('red',nargs=1,type=float, help='Redshift (float)')
parser.add_argument('noise',nargs=1,type=str,help='MockOb noise (str)')
parser.add_argument('--suffix',nargs=1,type=str, help='Additional suffix in SBprofile file (ex: redshift, snapnum) (float)')
parser.add_argument('--dm',nargs=1,type=float, help='Mass bin width in log10 (float)')
args   = parser.parse_args()
root   = args.path[0]
gridsf = args.grids[0]
mlow   = args.mlow[0]
red    = args.red[0]
noise  = args.noise[0]
suffix = ''
if args.suffix is not None:
    suffix = args.suffix[0]
dM     = 0.5
if args.dm is not None:
    dM = args.dm[0]
mhi    = mlow+dM

# ---

#LOAD
grids = genfromtxt(gridsf,names=True,delimiter=';')['gn']
if shape(grids) == ():
    grids = [grids]
SBtab = []
for g in grids:
    g  = int(g)
    SBf = [ f for f in os.listdir(root) if os.path.isfile(os.path.join(root,f)) and '_gn'+str(g)+'_SB_s1_z'+str(red)+'_Ng256_observed_cubex_neb_SBprofile.npz' in f and '_noise'+noise in f and suffix in f] 
    for f in SBf:
        print(f)
        data = load(root+'/'+f)
        SBtab.append(data['sbprofile'])

R = load(root+'/'+SBf[0])['binmean'] #same for all

# ---
#WRITE: if only 1 object, plot all three Los (so save all three), else save median and 84/16 perc.
fout = root+"/SBprofile_Eagle_observed_cubex_neb_z_"+str(red)+"_s1_n"+noise+"_Mmin_"+str(mlow)+"_Mhigh_"+str(mhi)+".npz"
print("Writing SB profile to "+fout)
if len(grids) == 1:
    savez(fout,binmean=R,sbprofileX=SBtab[0],sbprofileY=SBtab[1],sbprofileZ=SBtab[2])
else:
    #COMPUTE MEDIAN OF MEDIANS AND PRECENTILES
    SB  = nanmedian(SBtab,axis=0)
    P16 = nanpercentile(SBtab,16,axis=0)
    P84 = nanpercentile(SBtab,84,axis=0)
    
    savez(fout,binmean=R,sbprofile=SB,p16=P16,p84=P84)
