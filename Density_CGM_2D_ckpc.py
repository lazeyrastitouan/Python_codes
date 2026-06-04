#####################################################################################################################
#                                                                                                                   #
#       This codes gets the 2D projected density of CGM (inside virial radii of halos)                              #
#       particles in a temperature range in annuli ckpc                                                             #
#                                                                                                                   #
#       This codes does not verify if halos are main or subhalos. This distinction should be done beforhand         #
#                                                                                                                   #
#       Format of sims: gizmo/eagle hdf5 (needs keys tempeature, position and density in PartType0 group)           #
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
from decimal import Decimal

def insphere(self, part, halo, radius):
    diff = part-halo
    dist = linalg.norm(diff)

    return dist <= radius

# What to do
parser = argparse.ArgumentParser(description='Get projected 2D density in temperature range for CGM particles of a halo catalog. Args: < path to snapshot >, < simulation type >, < Halo catalog >')
parser.add_argument('snappath',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('simtype',nargs=1,type=str, help='Simulation type (dalya/eagle/colibre). Hello = dalya, rename hello snapshot by snap_<hello snapshot> (str)')
parser.add_argument('Hcat',nargs=1,type=str, help='Halo catalog (ascii)')
args     = parser.parse_args()
snappath = args.snappath[0]
sim      = args.simtype[0]

# settings
cattype  = 'P2C' #'AHF' #'P2C'
verbose  = True
ptype    = 'PartType0'
Tmin     = 1e3
Tmax     = 1e5
Rs       = array([100,200]) #ckpc 
Ro       = ['100','200']
LOS      = ['X','Y','Z']

print("Simulation type = "+sim)
print("Catalog type = "+cattype)
print("Tmin: "+str(Tmin)+", Tmax: "+str(Tmax))

#get density/T/position from simulation output
if sim == 'colibre':
    h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5' and f[0:7] == 'colibre' and len(f) > 17]
else:
    h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5' and f[0:4] == 'snap']

h5files  = sorted(h5files)
h5files  = sorted(h5files,key=len)
numfiles = len(h5files)
if verbose:
    print(h5files)

if sim == 'colibre':
    h    = h5py.File(snappath+'/'+h5files[0],'r')['Cosmology'].attrs['h'][0]
    #for units conversion
    U_M  = float(h5py.File(snappath+'/'+h5files[0],'r')['Units'].attrs['Unit mass in cgs (U_M)'][0])
    U_L  = float(h5py.File(snappath+'/'+h5files[0],'r')['Units'].attrs['Unit length in cgs (U_L)'][0])
    U_t  = float(h5py.File(snappath+'/'+h5files[0],'r')['Units'].attrs['Unit time in cgs (U_t)'][0])
    Lbox = float(h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['BoxSize'][0])*1e3 #ckpc
else:
    h = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['HubbleParam']
z = h5py.File(snappath+'/'+h5files[0],'r')['Header'].attrs['Redshift']
a = 1./(1.+z)

#convert bins to ckpc/h
Rs = Rs*h

if sim == 'eagle':
    #define a region as eagle is too big
    centerH = array([24191.,77352.,61463.])*h #ckpc/h
    dx      = 20000.*h
    print("Considering particles only in "+str(centerH-dx)+" <= x <= "+str(centerH+dx))
elif sim == 'colibre':
    #define a region as colibre is too big
    #centerH = array([70363.29,59403.1,69764.46])*h #ckpc/h
    centerH = array([Lbox/2.,Lbox/2.,Lbox/2.])*h #ckpc/h
    dx      = Lbox/2.*h
    print("Considering particles only in "+str(centerH-dx)+" <= x <= "+str(centerH+dx))

# Load particles
dens, temp, mass = empty(0),empty(0),empty(0)
pos = empty((0,3))

# --- HDF5 FILES LOOP --- #
print("Loading...")
for i in h5files:
    f   = h5py.File(snappath+"/"+i, "r")
    print(f)

    # --- POSITIONS --- #
    pos  = concatenate((pos,array(f[ptype]['Coordinates']))) #kpc/h if dalya, mpc/h if eagle, mpc if colibre
    
    # --- MASSES --- #
    mass = concatenate((mass,array(f[ptype]['Masses']))) #Msun/h if dalya, Msun if Colibre

    if sim == 'dalya':
        # --- DENSITY --- #
        dens = concatenate((dens,array(f[ptype]['Density'])))
        # --- TEMPERATURE --- #
        temp = concatenate((temp,array(f[ptype]['Temperature'])))

    elif sim == 'eagle':
        #get density from entropy, need a few params more for eagle to get density from entropy
        aexp_scale_exponent_S  = f[ptype]['Entropy'].attrs['aexp-scale-exponent']
        h_scale_exponent_S     = f[ptype]['Entropy'].attrs['h-scale-exponent']
        CGSConversionFactor_S  = f[ptype]['Entropy'].attrs['CGSConversionFactor']
        aexp_scale_exponent_U  = f[ptype]['InternalEnergy'].attrs['aexp-scale-exponent']
        h_scale_exponent_U     = f[ptype]['InternalEnergy'].attrs['h-scale-exponent']
        CGSConversionFactor_U  = f[ptype]['InternalEnergy'].attrs['CGSConversionFactor']
        aexp_scale_exponent_NH = f[ptype]['ElementAbundance/Hydrogen'].attrs['aexp-scale-exponent']
        h_scale_exponent_NH    = f[ptype]['ElementAbundance/Hydrogen'].attrs['h-scale-exponent']
        CGSConversionFactor_NH = f[ptype]['ElementAbundance/Hydrogen'].attrs['CGSConversionFactor']

        S     = array(f[ptype]['Entropy'])
        S     = S * a**aexp_scale_exponent_S * h**h_scale_exponent_S * CGSConversionFactor_S
        U     = array(f[ptype]['InternalEnergy'])
        U     = U * a**aexp_scale_exponent_U * h**h_scale_exponent_U * CGSConversionFactor_U
        NH    = array(f[ptype]['ElementAbundance/Hydrogen'])
        NH    = NH * a**aexp_scale_exponent_NH * h**h_scale_exponent_NH * CGSConversionFactor_NH
        rho_S = ((2./3.) * (U/S))**(3./2.)*NH  #gram / physical cm^-3 
            
        dens  = concatenate((dens,rho_S))
        del U, S, NH, rho_S
        # --- TEMPERATURE --- #
        temp = concatenate((temp,array(f[ptype]['Temperature'])))

    elif sim == 'colibre':
        a_exp_D   = float(f[ptype]['Densities'].attrs['a-scale exponent'][0])
        U_L_exp_D = float(f[ptype]['Densities'].attrs['U_L exponent'][0])
        U_M_exp_D = float(f[ptype]['Densities'].attrs['U_M exponent'][0])
        D_factor  = a**a_exp_D * U_L**U_L_exp_D * U_M**U_M_exp_D #a^-3 U_M U_L^-3  [ g pcm^-3 ]

        D     = array(f[ptype]['Densities'])
        D     = D * D_factor

        dens  = concatenate((dens,D))

        # --- TEMPERATURE --- #
        temp = concatenate((temp,array(f[ptype]['Temperatures'])))

# --- END OF HDF5 FILES LOOP --- #

# For Eagle, keep only particles in a sub-region 
if sim == 'eagle':
    pos   = pos*1e3 #kpc/h now
    maskE = where((pos[:,0]<=centerH[0]+dx) & (pos[:,0]>=centerH[0]-dx) & (pos[:,1]<=centerH[1]+dx) & (pos[:,1]>=centerH[1]-dx) & (pos[:,2]<=centerH[2]+dx) & (pos[:,2]>=centerH[2]-dx))
    pos   = pos[maskE]
    dens  = dens[maskE]
    temp  = temp[maskE]
    mass  = mass[maskE]
elif sim == 'colibre':
    pos   = pos*1e3*h #kpc/h now
    maskE = where((pos[:,0]<=centerH[0]+dx) & (pos[:,0]>=centerH[0]-dx) & (pos[:,1]<=centerH[1]+dx) & (pos[:,1]>=centerH[1]-dx) & (pos[:,2]<=centerH[2]+dx) & (pos[:,2]>=centerH[2]-dx))
    pos   = pos[maskE]
    dens  = dens[maskE]
    temp  = temp[maskE]
    mass  = mass[maskE]

# --- KEEP PARTICLES IN CORRECT TEMPERATURE RANGE --- #
maskT = where((temp >= Tmin) & (temp <= Tmax))
dens  = dens[maskT]
pos   = pos[maskT]
mass  = mass[maskT]

if verbose:
    np = len(dens)
    print("Number of particles: "+str(np))

if cattype == 'AHF':
    cat = loadtxt(args.Hcat[0],usecols=(5,6,7,11),dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
elif cattype == 'P2C':
    cat = loadtxt(args.Hcat[0],usecols=(1,2,3,4),delimiter=';',dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
    #convert to ckpc/h
    cat['Hx']   = cat['Hx']*h*1000.
    cat['Hy']   = cat['Hy']*h*1000.
    cat['Hz']   = cat['Hz']*h*1000.
    cat['Rvir'] = cat['Rvir']*h*(1.+z)

Rvir = cat['Rvir']
Hpos = array([cat['Hx'],cat['Hy'],cat['Hz']]).T
if sim == 'eagle' or sim == 'colibre':
    maskEH = where((Hpos[:,0]<=centerH[0]+dx) & (Hpos[:,0]>=centerH[0]-dx) & (Hpos[:,1]<=centerH[1]+dx) & (Hpos[:,1]>=centerH[1]-dx) & (Hpos[:,2]<=centerH[2]+dx) & (Hpos[:,2]>=centerH[2]-dx))
    Rvir   = Rvir[maskEH]
    Hpos   = Hpos[maskEH]
if shape(Rvir) == ():
    nh = 1
else:
    nh   = len(Rvir)
if verbose:
    print("Kept "+str(nh)+" halos")
del cat

# --- KEEP PARTICLES IN VIRIAL RADII --- #
print("Keeping particles inside the CGM")
if nh == 1:
    distall = [linalg.norm(posall-Hpos) for posall in pos]
else:
    distall = [linalg.norm(posall-Hpos, axis=1) for posall in pos]
distall = array(distall)
incgm   = where(distall<=Rvir,1,0)
if nh == 1:
    maskR   = where(incgm==1)
else:
    maskR   = where(any(incgm==1,axis=1))
pos     = pos[maskR]
dens    = dens[maskR]
mass    = mass[maskR]
print("Kept "+str(len(dens))+" particles")

# --- LOOP OVER LOS AND RVs ---#
for los in LOS:
    il = LOS.index(los)

    #compute projected distance
    print("Computing projected distances along "+los)
    newpos  = array([pos[:,(il+1)%3],pos[:,(il+2)%3]]).T
    if nh == 1:
        newHpos = array([Hpos[(il+1)%3],Hpos[(il+2)%3]]).T #keep only indices perpendicular to LOS
        dist    = [linalg.norm(posp-newHpos) for posp in newpos]
    else:
        newHpos = array([Hpos[:,(il+1)%3],Hpos[:,(il+2)%3]]).T 
        dist    = [linalg.norm(posp-newHpos,axis=1) for posp in newpos]
    
    dist    = array(dist)

    print("Selecting particles...")
    for i in range(len(Rs)-1):
        #restrict to wanted particles
        inorout = where((dist>=Rs[i]) & (dist<=Rs[i+1]),1,0) #1 if within the projected annulus, 0 otherwise
        if nh == 1:
            selected = where(inorout==1)
        else:
            selected = where(any(inorout==1,axis=1))
        if verbose:
            print("shape(selected) = "+str(shape(selected)))
        newdens = dens[selected]  
        newmass = mass[selected]

        #Now let's transform density untis to hydrogen atoms per physical cc
        mH         = 1.6735575e-24 #gram
        mHe        = 6.6464731e-24 #gram
        mA         = mH #0.76*mH+0.24*mHe
        Msuningram = 1.989e43 # 10^10 Msun in g
        kpcincm    = 3.086e21
        if sim == 'dalya':
            newdens = newdens*a**(-3.)*h**2.*Msuningram/kpcincm**3./mA 
            newmass = newmass/h*Msuningram/mA
        elif sim == 'eagle'or sim =='colibre':
            newdens = newdens/mA
            newmass = newmass*Msuningram/mA

        #let's get volume as needed for weighted histograms later
        vol = newmass/newdens #physical cc
        vol = vol/(kpcincm**3.) #physical kpc
        
        #save
        fout = 'density_volume_Tmin'+'{:.0e}'.format(Decimal(Tmin))+'_Tmax'+'{:.0e}'.format(Decimal(Tmax))+'_in_'+Ro[i]+'-'+Ro[i+1]+'ckpc_'+str(los)
        print("Saving density to "+snappath+'/'+fout)
        savez(snappath+'/'+fout, density=newdens, volume=vol)

        del newdens, newmass, vol

    del newHpos, newpos, dist
