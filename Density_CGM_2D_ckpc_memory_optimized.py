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
dustcut  = 1e-5 #dust mass fraction threshold. Set to None for no cut
ptype    = 'PartType0'
Tmin     = 1e3
Tmax     = 1e5
Rs       = array([50,160]) #ckpc 
Ro       = ['50','160']
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
# Store one array per HDF5 file and concatenate once at the end.
# This avoids repeated full-array reallocations/copies inside the file loop.
dens_list, temp_list, mass_list, pos_list, dust_list = [], [], [], [], []

# --- HDF5 FILES LOOP --- #
print("Loading...")
for i in h5files:
    with h5py.File(snappath+"/"+i, "r") as f:
        print(f)

        # --- POSITIONS --- #
        pos_list.append(array(f[ptype]['Coordinates'])) #kpc/h if dalya, mpc/h if eagle, mpc if colibre
        
        # --- MASSES --- #
        mass_list.append(array(f[ptype]['Masses'])) #Msun/h if dalya, Msun if Colibre

        if sim == 'dalya':
            # --- DENSITY --- #
            dens_list.append(array(f[ptype]['Density']))
            # --- TEMPERATURE --- #
            temp_list.append(array(f[ptype]['Temperature']))

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
                
            dens_list.append(rho_S)
            del U, S, NH, rho_S
            # --- TEMPERATURE --- #
            temp_list.append(array(f[ptype]['Temperature']))

        elif sim == 'colibre':
            a_exp_D   = float(f[ptype]['Densities'].attrs['a-scale exponent'][0])
            U_L_exp_D = float(f[ptype]['Densities'].attrs['U_L exponent'][0])
            U_M_exp_D = float(f[ptype]['Densities'].attrs['U_M exponent'][0])
            D_factor  = a**a_exp_D * U_L**U_L_exp_D * U_M**U_M_exp_D #a^-3 U_M U_L^-3  [ g pcm^-3 ]

            D     = array(f[ptype]['Densities'])
            D     = D * D_factor

            dens_list.append(D)

            # --- TEMPERATURE --- #
            temp_list.append(array(f[ptype]['Temperatures']))

            # --- DUST --- #
            if dustcut is not None:
                dust_list.append(array(f[ptype]['TotalDustMassFractions']))
                

# Concatenate once, after reading all HDF5 files.
pos  = concatenate(pos_list, axis=0)
mass = concatenate(mass_list, axis=0)
dens = concatenate(dens_list, axis=0)
temp = concatenate(temp_list, axis=0)
if dustcut is not None:
    dust = concatenate(dust_list, axis=0)
del pos_list, mass_list, dens_list, temp_list, dust_list

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
    if dustcut is not None:
        dust = dust[maskE]

# --- KEEP PARTICLES IN CORRECT TEMPERATURE RANGE --- #
maskT = where((temp >= Tmin) & (temp <= Tmax))
dens  = dens[maskT]
pos   = pos[maskT]
mass  = mass[maskT]
if dustcut is not None:
    dust = dust[maskT]
    
if verbose:
    n_particles = len(dens)
    print("Number of particles: "+str(n_particles))

# --- KEEP PARTICLES W/ CORRECT DUST MASS FRACTION --- #
if sim == 'colibre' and dustcut is not None:
    maskdst = where(dust <= dustcut)
    pos   = pos[maskdst]
    dens  = dens[maskdst]
    temp  = temp[maskdst]
    mass  = mass[maskdst]

# --- HALO CATALOGUE --- #
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

# Force a consistent one-dimensional halo catalogue, even if only one halo is present.
Rvir = atleast_1d(Rvir)
Hpos = atleast_2d(Hpos)
nh   = len(Rvir)
if verbose:
    print("Kept "+str(nh)+" halos")
del cat

# --- KEEP PARTICLES IN VIRIAL RADII --- #
print("Keeping particles inside the CGM")

# Avoid building an enormous (N_particles, N_halos) distance array.
# Instead, loop over halos and only compute distances for particles inside
# a cheap cubic preselection around each halo.
keep = zeros(len(pos), dtype=bool)

for ih in range(nh):
    hx, hy, hz = Hpos[ih]
    rv = Rvir[ih]

    box = (
        (pos[:,0] >= hx-rv) & (pos[:,0] <= hx+rv) &
        (pos[:,1] >= hy-rv) & (pos[:,1] <= hy+rv) &
        (pos[:,2] >= hz-rv) & (pos[:,2] <= hz+rv)
    )

    idx = where(box)[0]
    if len(idx) == 0:
        continue

    dxp = pos[idx,0] - hx
    dyp = pos[idx,1] - hy
    dzp = pos[idx,2] - hz

    inside = dxp*dxp + dyp*dyp + dzp*dzp <= rv*rv
    keep[idx[inside]] = True

maskR = where(keep)[0]
pos   = pos[maskR]
dens  = dens[maskR]
mass  = mass[maskR]
print("Kept "+str(len(dens))+" particles")

# --- LOOP OVER LOS AND RVs ---#
for los in LOS:
    il = LOS.index(los)

    ax1 = (il+1)%3
    ax2 = (il+2)%3

    #compute projected distance
    print("Computing projected distances along "+los)

    print("Selecting particles...")
    for i in range(len(Rs)-1):
        # Avoid building an enormous (N_particles, N_halos) projected-distance array.
        # Loop over halos and compute projected distances only after a cheap 2D box preselection.
        selected_mask = zeros(len(pos), dtype=bool)

        rmin2 = Rs[i]*Rs[i]
        rmax2 = Rs[i+1]*Rs[i+1]
        rmax  = Rs[i+1]

        for ih in range(nh):
            h1 = Hpos[ih,ax1]
            h2 = Hpos[ih,ax2]

            box = (
                (pos[:,ax1] >= h1-rmax) & (pos[:,ax1] <= h1+rmax) &
                (pos[:,ax2] >= h2-rmax) & (pos[:,ax2] <= h2+rmax)
            )

            idx = where(box)[0]
            if len(idx) == 0:
                continue

            d1  = pos[idx,ax1] - h1
            d2  = pos[idx,ax2] - h2
            rp2 = d1*d1 + d2*d2

            in_annulus = (rp2 >= rmin2) & (rp2 <= rmax2)
            selected_mask[idx[in_annulus]] = True

        selected = where(selected_mask)[0]

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

        del selected_mask, selected, newdens, newmass, vol
