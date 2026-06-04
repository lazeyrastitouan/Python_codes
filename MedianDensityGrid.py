from numpy import*
import h5py
import glob
import argparse
from decimal import Decimal
from collections import defaultdict

def weighted_median(values, weights):
    i = argsort(values)
    c = cumsum(weights[i])
    return values[i[searchsorted(c, 0.5 * c[-1])]]

parser = argparse.ArgumentParser(description='Get median density of particles in a temperature range along a LOS CGM particles of a halo catalog. Args: < simulation type >, < Halo catalog >, < sub-grid 1D length >, < sub-grid 1D size >, < Tmin >, < Tmax >')
parser.add_argument('-simtype',nargs=1,type=str, help='Simulation type (dalya/colibre/hello). For hello, rename napshot by snap_<hello snapshot> and make sure only one is named like this (str)')
parser.add_argument('-subbox',nargs=1,type=float,default=2000., help='sub-box 1D length (ckpc)')
parser.add_argument('-ngrid',nargs=1,type=int,default=256., help='subgrid 1D size')
parser.add_argument('-tmin',nargs=1,type=float,default=1e3, help='min temperature (K)')
parser.add_argument('-tmax',nargs=1,type=float,default=1e5, help='max temperature (K)')
parser.add_argument('-mmin',nargs=1,type=float,default=12.0, help='lower mass edge (log10)')
args     = parser.parse_args()
sim      = args.simtype[0]
subL     = args.subbox[0]
grid     = args.ngrid[0]
Tmin     = args.tmin[0]
Tmax     = args.tmax[0]
mmin     = args.mmin[0]
cell     = subL/grid #ckpc

Rs       = array([50,160]) #ckpc 
LOS      = ['X','Y','Z']
verbose  = False

#Units transformation
mH         = 1.6735575e-24 #gram
mHe        = 6.6464731e-24 #gram
mA         = mH #0.76*mH+0.24*mHe
Msuningram = 1.989e43 # 10^10 Msun in g
kpcincm    = 3.086e21
ptype      = 'PartType0'


print("Tmin: "+str(Tmin)+", Tmax: "+str(Tmax))

#LOAD
if sim == 'colibre':
    path    = '/data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/'
    h5files = glob.glob(path+'colibre*.hdf5')
    
    #halos
    Hcat        = '/data/titouan/COLIBRE/L0050N0752/z3p5/200_snap_z3.5_halos_sorted_main_hires_Mlow_'+str(mmin)+'_Mhigh_'+str(mmin+0.5)+'.txt'
    cat         = loadtxt(Hcat,usecols=(1,2,3,4),delimiter=';',dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
    #convert to ckpc (except Rvir for which we need z, done later)
    cat['Hx']   =   cat['Hx']*1000.
    cat['Hy']   = cat['Hy']*1000.
    cat['Hz']   = cat['Hz']*1000.
    cat['Rvir'] = cat['Rvir']
    Rvir = cat['Rvir']
    Hpos = array([cat['Hx'],cat['Hy'],cat['Hz']]).T

elif sim == 'dalya':
    path    = '/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/'
    h5files = glob.glob(path+'Halo10/sftresh_cool_low_fb_boost/snapdir_004/snap*.hdf5')+glob.glob(path+'Halo1312/sftresh_cool_low_fb_boost/snapdir_001/snap*.hdf5')
    Hcat    = [path+'Halo10/sftresh_cool_low_fb_boost/AHF_004_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_'+str(mmin)+'_Mhigh_'+str(mmin+0.5)+'.txt',path+'Halo1312/sftresh_cool_low_fb_boost/AHF_001_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_'+str(mmin)+'_Mhigh_'+str(mmin+0.5)+'.txt']

    #halos
    Hpos = empty((0,3))
    Rvir = empty(0)
    for hcat in Hcat:
        cat = loadtxt(hcat,usecols=(1,2,3,4),delimiter=';',dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
        #convert to ckpc (except Rvir for which we need z, done later)
        Hx = cat['Hx']*1000.
        Hy = cat['Hy']*1000.
        Hz = cat['Hz']*1000.
        HR = cat['Rvir']

        Htemp = array([Hx,Hy,Hz]).T

        Rvir  = concatenate((Rvir,HR))
        Hpos  = concatenate((Hpos,Htemp))

else:#hello
    path    = '/data/titouan/HELLO/CDM/g9.20e12/'
    h5files = glob.glob(path+'snap_*.hdf5')
    
    #halos
    Hcat      = '/data/titouan/HELLO/CDM/g9.20e12/g9.20e12.00128.z3.596.AHF_halos_sorted_main_hires_Mlow_'+str(mmin)+'_Mhigh_'+str(mmin+0.5)+'.txt'
    cat       = loadtxt(Hcat,usecols=(1,2,3,4),delimiter=';',dtype={'names':('Hx','Hy','Hz','Rvir'), 'formats':(float,float,float,float)})
    #convert to ckpc (except Rvir for which we need z, done later)
    cat['Hx'] = cat['Hx']*1000.
    cat['Hy'] = cat['Hy']*1000.
    cat['Hz'] = cat['Hz']*1000.
    Rvir      = cat['Rvir']
    Hpos      = array([cat['Hx'],cat['Hy'],cat['Hz']]).T

if shape(Rvir) == ():
    nh = 1
else:
    nh   = len(Rvir)
print('Found '+str(nh)+' halos')
del cat

h5files  = sorted(h5files)
h5files  = sorted(h5files,key=len)
numfiles = len(h5files)
if verbose:
        print(h5files)

if sim == 'colibre':
    h    = h5py.File(h5files[0],'r')['Cosmology'].attrs['h'][0]
    #for units conversion
    U_M  = float(h5py.File(h5files[0],'r')['Units'].attrs['Unit mass in cgs (U_M)'][0])
    U_L  = float(h5py.File(h5files[0],'r')['Units'].attrs['Unit length in cgs (U_L)'][0])
    U_t  = float(h5py.File(h5files[0],'r')['Units'].attrs['Unit time in cgs (U_t)'][0])
else:
    h = h5py.File(h5files[0],'r')['Header'].attrs['HubbleParam']
z = h5py.File(h5files[0],'r')['Header'].attrs['Redshift']
a = 1./(1.+z)

#Now that we have z we can convert Rvir to ckpc
Rvir = Rvir*(1.+z)

# Load particles
dens, temp, mass = empty(0),empty(0),empty(0)
pos = empty((0,3))

# --- HDF5 FILES LOOP --- #
print("Loading...")
for i in h5files:
    f   = h5py.File(i, "r")
    print(f)

    if sim == 'colibre':
        # --- POSITIONS --- #
        pos  = concatenate((pos,array(f[ptype]['Coordinates'])*U_L/kpcincm)) #ckpc

        # --- MASSES --- #
        mass = concatenate((mass,array(f[ptype]['Masses']))) #10^10 Msun

        # --- DENSITIES --- #
        a_exp_D   = float(f[ptype]['Densities'].attrs['a-scale exponent'][0])
        U_L_exp_D = float(f[ptype]['Densities'].attrs['U_L exponent'][0])
        U_M_exp_D = float(f[ptype]['Densities'].attrs['U_M exponent'][0])
        D_factor  = a**a_exp_D * U_L**U_L_exp_D * U_M**U_M_exp_D #a^-3 U_M U_L^-3  [ g pcm^-3 ]

        D     = array(f[ptype]['Densities'])
        D     = D * D_factor #physical cgs

        dens  = concatenate((dens,D))

        # --- TEMPERATURE --- #
        temp = concatenate((temp,array(f[ptype]['Temperatures'])))

    else:
        # --- POSITIONS --- #
        pos  = concatenate((pos,array(f[ptype]['Coordinates'])/h)) #ckpc

        # --- MASSES --- #
        mass = concatenate((mass,array(f[ptype]['Masses'])/h)) #1e10 Msun

        # --- DENSITY --- #
        dens = concatenate((dens,array(f[ptype]['Density'])))
        
        # --- TEMPERATURE --- #
        temp = concatenate((temp,array(f[ptype]['Temperature'])))


# --- KEEP PARTICLES IN CORRECT TEMPERATURE RANGE --- #
maskT = where((temp >= Tmin) & (temp <= Tmax))
dens  = dens[maskT]
pos   = pos[maskT]
mass  = mass[maskT]

if verbose:
    np = len(dens)
    print("Number of particles: "+str(np))

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

#convert mass and density to atom per physical cm^3
if sim == 'dalya' or sim == 'hello':
    dens = dens*a**(-3.)*h**2.*Msuningram/kpcincm**3./mA
else:
    dens = dens/mA
mass = mass*Msuningram/mA

# --- LOOP OVER LOS AND HALOS ---#
#final table of medians and std_dev that will be filled for each 2D pixel (along each LOS) for each halo
med_in_pix, mean_in_pix, std_in_pix = [], [], []

for los in LOS:
    il = LOS.index(los)

    #compute projected distance
    print("Computing projected distances along "+los)
    newpos  = array([pos[:,(il+1)%3],pos[:,(il+2)%3]]).T #keep only indices perpendicular to LOS
    if nh == 1:
        newHpos = array([Hpos[(il+1)%3],Hpos[(il+2)%3]]).T #keep only indices perpendicular to LOS
        dist    = [linalg.norm(posp-newHpos) for posp in newpos]
    else:
        newHpos = array([Hpos[:,(il+1)%3],Hpos[:,(il+2)%3]]).T
        dist    = [linalg.norm(posp-newHpos,axis=1) for posp in newpos]

    dist    = array(dist).T #distance of each particle to each center

    #find which particles are in which halo 
    inorout = where((dist>=Rs[0]) & (dist<=Rs[1]),1,0) #1 if within the projected annulus, 0 otherwise
    if nh == 1:
        inorout = [inorout]
        newHpos = [newHpos]

    #NGP grid assignment
    print("NGP assignment for each  halo")
    for ih in range(nh):
        in_this_halo   = inorout[ih] == 1
        
        #position of particles in ih-th halo in ckpc
        part_this_halo = newpos[in_this_halo]

        #shift coord to have halo at 0,0
        part_this_halo = part_this_halo - newHpos[ih]

        #ij pair corresponding to position on grid
        ij             = (part_this_halo/cell).round() #maybe add a np.clip 

        if verbose:
            print("There are "+str(shape(ij)[0])+" particles in this halo's annulus")
    
        #get density and volume
        densin   = dens[in_this_halo]
        massin   = mass[in_this_halo]
        volume   = massin/densin

        #compute volume-weighted median density for all particles with same ij
        same_cells = defaultdict(list)
        for ic, pair in enumerate(map(tuple, ij)):
            same_cells[pair].append(ic)

        index_ij = list(same_cells.values())

        #compute volume-weighted median density and std dev for all particles with same ij (i.e. in each projected cell) and add to table
        for ic in index_ij:
            med_in_pix.append(weighted_median(densin[ic],volume[ic]))
            mean = average(densin[ic], weights=volume[ic])
            mean_in_pix.append(mean)
            # Fast and numerically precise:
            std_in_pix.append(average((densin[ic]-mean)**2, weights=volume[ic]))
        
        del part_this_halo, ij, densin, massin, volume, same_cells, index_ij

    del newpos, newHpos, dist, inorout

fout = path+'medians_stds_Mlow_'+str(mmin)+'_Mhigh_'+str(mmin+0.5)+'_R_'+str(Rs[0])+'-'+str(Rs[1])+'.npz'
print("Saving medians and stds to "+fout)
savez(fout,median=med_in_pix,mean=mean_in_pix,std=std_in_pix)
