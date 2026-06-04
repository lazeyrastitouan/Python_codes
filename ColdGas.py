#####################################################################################################################
#                                                                                                                   #
#       This codes plots the cold gas fraction and volume filling factor profiles for a given halo catalog          #
#       This codes does not verify if halos are main or subhalos. This distinction should be done beforhand         #
#       Format of halo catalogs: AHF/P2C ASCII (at least 12 columns: Host, x, y, z, Rvir in col 2, 6, 7, 8, 12)     #  
#                                                                                                                   #
#       T. Lazeyras, 2025                                                                                           #
#                                                                                                                   #
#####################################################################################################################

from numpy import*
import sys
import os
import h5py
import argparse
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib
#from pylab import *
from matplotlib.colors import* 
from matplotlib.colors import LogNorm


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

def insphere(self, part, halo, radius):
    diff = part-halo
    dist = linalg.norm(diff)

    return dist <= radius

# What to do
parser = argparse.ArgumentParser(description='Get cold gas fraction and volume filling factor profilesfor a given halo catalog. Args: < path to snapshot (hdf5) >, < halo catalog (str) >, < simulation type (str) >, < Optional, density threshold (float) >')
parser.add_argument('snappath1',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('Hcat1',nargs=1,type=str, help='Halo catalog (ascii), in that case plots TvsD inside/outside of halos')
parser.add_argument('simtype',nargs=1,type=str,help='Simulation type (dalya/eagle) (str)')
parser.add_argument('--nstar',nargs=1,type=float,help='Density threshold in atom per cc (optional) (float)')
parser.add_argument('--snappath2',nargs=1,type=str, help='Optional for dalya, second path to snapshot (hdf5) (str)')
parser.add_argument('--Hcat2',nargs=1,type=str, help='Optional for dalya, second halo catalog (ascii) (str)')
args      = parser.parse_args()
snappaths = [args.snappath1[0]]
sim       = args.simtype[0]
catfiles  = [args.Hcat1[0]]
if args.Hcat2 is not None:
    catfiles.append(args.Hcat2[0])
if args.snappath2 is not None:
    snappaths.append(args.snappath2[0])

datatype = 'HDF5' #'Gadget2'
cattype  = 'P2C' #'AHF' #'P2C'
verbose  = True
savenpy  = True
plot     = False
#fRv      = [0.0,0.2,0.5,0.7,1.0]
Rs       = array([50,160,250,350,700]) #ckpc 
Ro       = ['50','160','250','350','700']
Tmax     = 1e5 #K
root_out = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/'

print("Simulation type = "+sim)
print("Data type = "+datatype)
print("Catalog type = "+cattype)

# ---- LOAD PARTICLES ---- #
print("Loading...")
dens, temp, massG, massDM = empty(0),empty(0),empty(0),empty(0)
posG, posDM = empty((0,3)), empty((0,3))
Hpos = empty((0,3))
Rvir = empty(0)

#define a region as eagle is too big
if sim =='eagle':
    centerH = array([24191.,77352.,61463.]) #ckpc
    dx      = 20000.
    print("Considering particles only in "+str(centerH-dx)+" <= x <= "+str(centerH+dx)+" ckpc")

#get Temperature, masses and density (so volume) from simulation output
for snappath in snappaths:
    sim_i   = snappaths.index(snappath)
    catfile = catfiles[sim_i]

    h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5' and f[0:4] == 'snap']
    h5files  = sorted(h5files)
    h5files  = sorted(h5files,key=len)
    numfiles = len(h5files)
    if verbose:
        print(h5files)

    h  = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['HubbleParam']
    z  = round(h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['Redshift'],1)
    a  = 1./(1.+z)

    if sim == 'eagle':
        Nptot  = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['NumPart_Total'][1]
        massDM = concatenate((massDM,array([h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['MassTable'][1]]*Nptot)))

    # --- HDF5 FILES LOOP --- #
    for i in h5files:
        f   = h5py.File(snappath+"/"+i, "r")
        print(f)

        # --- POSITIONS, MASS AND DENSITY --- #
        if sim == 'dalya':
            posG   = concatenate((posG,array(f['PartType0']['Coordinates'])/h+1e6*sim_i)) #ckpc, the +1e6*sim_i is to translate particles from the second box by 1cGpc to not associate them to halos of the first box
            posDM  = concatenate((posDM,array(f['PartType1']['Coordinates'])/h+1e6*sim_i)) #ckpc
            dens   = concatenate((dens,array(f['PartType0']['Density'])))
            massG  = concatenate((massG,array(f['PartType0']['Masses'])))
            massDM = concatenate((massDM,array(f['PartType1']['Masses'])))
    #–––––––––––––––––––––––––––––––––––––––––––––
        elif sim == 'eagle':
            posG  = concatenate((posG,array(f['PartType0']['Coordinates'])*1e3/h+1e6*sim_i)) #ckpc ; for eagle sim_i = 0 only always so pos no modified
            posDM = concatenate((posDM,array(f['PartType1']['Coordinates'])*1e3/h+1e6*sim_i)) #ckpc
            #get density from entropy, need a few params more for eagle to get density from entropy
            aexp_scale_exponent_S  = f['PartType0']['Entropy'].attrs['aexp-scale-exponent']
            h_scale_exponent_S     = f['PartType0']['Entropy'].attrs['h-scale-exponent']
            CGSConversionFactor_S  = f['PartType0']['Entropy'].attrs['CGSConversionFactor']
            aexp_scale_exponent_U  = f['PartType0']['InternalEnergy'].attrs['aexp-scale-exponent']
            h_scale_exponent_U     = f['PartType0']['InternalEnergy'].attrs['h-scale-exponent']
            CGSConversionFactor_U  = f['PartType0']['InternalEnergy'].attrs['CGSConversionFactor']
            aexp_scale_exponent_NH = f['PartType0']['ElementAbundance/Hydrogen'].attrs['aexp-scale-exponent']
            h_scale_exponent_NH    = f['PartType0']['ElementAbundance/Hydrogen'].attrs['h-scale-exponent']
            CGSConversionFactor_NH = f['PartType0']['ElementAbundance/Hydrogen'].attrs['CGSConversionFactor']

            S     = array(f['PartType0']['Entropy'])
            S     = S * a**aexp_scale_exponent_S * h**h_scale_exponent_S * CGSConversionFactor_S
            U     = array(f['PartType0']['InternalEnergy'])
            U     = U * a**aexp_scale_exponent_U * h**h_scale_exponent_U * CGSConversionFactor_U
            NH    = array(f['PartType0']['ElementAbundance/Hydrogen'])
            NH    = NH * a**aexp_scale_exponent_NH * h**h_scale_exponent_NH * CGSConversionFactor_NH
            rho_S = ((2./3.) * (U/S))**(3./2.)*NH  #gram / physical cm^-3 
            
            dens  = concatenate((dens,rho_S))
            del U, S, NH, rho_S
    
            # --- MASS ---#
            massG  = concatenate((massG,array(f['PartType0']['Mass'])))

        # --- TEMPERATURE --- #
        temp = concatenate((temp,array(f['PartType0']['Temperature'])))

    # --- END OF HDF5 FILES LOOP --- #

    # --- HALO CATALOG --- #
    if cattype == 'AHF':
        cat = loadtxt(catfile,usecols=(5,6,7,11),dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
        #convert to ckpc
        cat['Hx']   = cat['Hx']/h
        cat['Hy']   = cat['Hy']/h
        cat['Hz']   = cat['Hz']/h
        cat['Rvir'] = cat['Rvir']/h
    elif cattype == 'P2C':
        cat = loadtxt(catfile,usecols=(1,2,3,4),delimiter=';',dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
        #convert to ckpc
        cat['Hx']   = cat['Hx']*1000.
        cat['Hy']   = cat['Hy']*1000.
        cat['Hz']   = cat['Hz']*1000.
        cat['Rvir'] = cat['Rvir']*(1.+z)

        Rvir = concatenate((Rvir,cat['Rvir']))
        Hpos = concatenate((Hpos,array([cat['Hx'],cat['Hy'],cat['Hz']]).T+1e6*sim_i))

# --- END OF LOADING --- #

#If density threshold, transform density units to hydrogen atoms per physical cc
if args.nstar is not None:
    nstar = args.nstar[0]
    print("Selecting cells with density < "+str(nstar))
    mH  = 1.6735575e-24 #gram
    mHe = 6.6464731e-24 #gram
    mA  = mH #0.76*mH+0.24*mHe
    if sim == 'dalya':
        kpcincm    = 3.086e21
        Msuningram = 1.989e43
        dens = dens*a**(-3.)*h**2.*Msuningram/kpcincm**3./mA
    elif sim == 'eagle':
        dens = dens/mA

    maskd = where(dens <= nstar)
    dens  = dens[maskd]
    posG  = posG[maskd]
    temp  = temp[maskd]
    massG = massG[maskd]

# For Eagle, keep only particles in a sub-region 
if sim == 'eagle':
    maskEG  = where((posG[:,0]<=centerH[0]+dx) & (posG[:,0]>=centerH[0]-dx) & (posG[:,1]<=centerH[1]+dx) & (posG[:,1]>=centerH[1]-dx) & (posG[:,2]<=centerH[2]+dx) & (posG[:,2]>=centerH[2]-dx))
    maskEDM = where((posDM[:,0]<=centerH[0]+dx) & (posDM[:,0]>=centerH[0]-dx) & (posDM[:,1]<=centerH[1]+dx) & (posDM[:,1]>=centerH[1]-dx) & (posDM[:,2]<=centerH[2]+dx) & (posDM[:,2]>=centerH[2]-dx))
    posG    = posG[maskEG]
    posDM   = posDM[maskEDM]
    dens    = dens[maskEG]
    temp    = temp[maskEG]
    massG   = massG[maskEG]
    massDM  = massDM[maskEDM]

    maskEH = where((Hpos[:,0]<=centerH[0]+dx) & (Hpos[:,0]>=centerH[0]-dx) & (Hpos[:,1]<=centerH[1]+dx) & (Hpos[:,1]>=centerH[1]-dx) & (Hpos[:,2]<=centerH[2]+dx) & (Hpos[:,2]>=centerH[2]-dx))
    Rvir   = Rvir[maskEH]
    Hpos   = Hpos[maskEH]

if verbose:
    np = len(dens)
    print("Number of gas particles: "+str(np))

if shape(Rvir) == ():
    nh = 1
else:
    nh   = len(Rvir)
if verbose:
    print("Got "+str(nh)+" halos")
del cat

# ---- GET DISTANCES ----#
print("Computing distances...")

if nh ==1:
    distG  = [linalg.norm(pos-Hpos) for posp in pos]
    distDM = [linalg.norm(pos-Hpos) for posp in posDM]
else:
    distG  = [linalg.norm(posp-Hpos,axis=1) for posp in posG]
    distDM = [linalg.norm(posp-Hpos,axis=1) for posp in posDM]
distG  = array(distG)
distDM = array(distDM)

del posG, posDM 

# ---- GET PROFILES ----#
print("Computing fV and fcold profiles...")
#nbins     = len(fRv) - 1
nbins     = len(Rs) - 1
fV, fcold = empty(nbins), empty(nbins)

if sim == 'dalya':
    Ob = h5py.File(snappaths[sim_i]+'/'+h5files[0],'r')['Header'].attrs['Omega_Baryon']
    Om = h5py.File(snappaths[sim_i]+'/'+h5files[0],'r')['Header'].attrs['Omega_Matter']
if sim == 'eagle':
    Om = h5py.File(snappaths[0]+'/'+h5files[0],'r')['Header'].attrs['Omega0']
    Ob = h5py.File(snappaths[0]+'/'+h5files[0],'r')['Header'].attrs['OmegaBaryon']

for i in range(nbins):
    #print('bin ' +str(fRv[i])+' - '+str(fRv[i+1]))
    print('bin ' +str(Rs[i])+' - '+str(Rs[i+1]))
    #inoroutG   = logical_and(fRv[i]*Rvir < distG, distG <= fRv[i+1]*Rvir)
    #first condition to make sure to be in cgm
    inoroutG   = (distG <= Rvir) & ((Rs[i] <= distG) & (distG < Rs[i+1]))
    selectedG  = any(inoroutG,axis=1)
    #inoroutDM  = logical_and(fRv[i]*Rvir < distDM, distDM <= fRv[i+1]*Rvir)
    inoroutDM  = (distDM <= Rvir) & ((Rs[i] <= distDM) & (distDM < Rs[i+1]))
    selectedDM = any(inoroutDM,axis=1)
    mG     = massG[selectedG]
    mtot   = sum(massDM[selectedDM])+sum(mG)
    voltot = mG/dens[selectedG] #units of mg are not same as mass units in dens but it doesnt matter as we will tke ratios
    
    #get cold gas
    cold     = where(temp[selectedG] <= Tmax)
    mcold    = sum(mG[cold])
    fV[i]    = sum(voltot[cold])/sum(voltot)
    fcold[i] = mcold/(Ob/Om*mtot)

    del inoroutG, selectedG, inoroutDM, selectedDM, mG, mtot, voltot, cold, mcold

#save if wanted
if savenpy:
    print("Saving profiles to "+root_out+'/fV_fcold_in_annulickpc_z'+str(z)+'.npz')
    if args.nstar is not None:
        savez(root_out+'/fV_fcold_in_Rvirfrac_Hcat_nstar_'+str(nstar), bins=fRv, fV=fV, fcold=fcold)
    else:
        #savez(snappath+'/fV_fcold_in_Rvirfrac_Hcat', bins=fRv, fV=fV, fcold=fcold)
        savez(root_out+'/fV_fcold_in_annulickpc_z'+str(z), bins=Rs, fV=fV, fcold=fcold)
    
if plot:
    #plt.stairs(fV,fRv,color='b',fill=False,label='fV')
    plt.stairs(fV,Rs,color='b',fill=False,label='fV')
    #plt.stairs(fcold,fRv,color='cyan',fill=False,label='fcold')
    plt.stairs(fcold,Rs,color='cyan',fill=False,label='fcold')
    plt.legend(loc=0)
    plt.xlabel('r [Rvir]')
    plt.ylabel('fraction')
    plt.show()
