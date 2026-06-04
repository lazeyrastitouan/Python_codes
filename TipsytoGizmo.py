from numpy import*
import pynbody as pn
import sys
import os
import h5py
import argparse
from multiprocessing import Pool
from tipsyutils import*

# SETTINGS
parser = argparse.ArgumentParser(description='Transform Tipsy format to Gizmo hdf5 format. Args: < Input/output folder (str) >, < Simulation name (tipsy) (str) >, < Tipsy parameterfile (str) >, < Omega baryon (float) >')
parser.add_argument('folder',nargs=1,type=str, help='Input/output folder (str)')
parser.add_argument('simulation',nargs=1,type=str, help='Tipsy simulation (str)')
parser.add_argument('paramfile',nargs=1,type=str, help='Tipsy parameter file (str)')
parser.add_argument('omegab',nargs=1,type=str, help='Omega baryon (float)')
args      = parser.parse_args()
folder    = args.folder[0]
simname   = args.simulation[0]
paramfile = args.paramfile[0]
Ob 	      = float(args.omegab[0])
verbose   = True

# READ
print("Loading Tipsy file and obtaining units conversion factors...")
sim        = pn.load(folder+'/'+simname)
params     = readparams(folder+'/'+paramfile) #type str:str for the moment...
props      = sim.properties
h          = props['h']
box        = float(params['dKpcUnit'])*h #boxsize in ckpc/h
cool       = int(bool(int(params['iGasModel'][0])))
a          = props['a']
red        = 1./a - 1.
ng         = len(sim.gas)
ndm        = len(sim.dm)
ns         = len(sim.star)
Lincgs     = 3.0856776e+21/h #kpc/h in cgs units. We hardcode the units parameters since we transform tipsy to these units
Mincgs     = 1.989e+43/h # 10^10 Msol in cgs
Vincgs     = 100000.0 #km/s in cgs
timeincgs  = Lincgs/Vincgs
Eincgs     = Mincgs*Lincgs**2/timeincgs**2
mproton    = 1.6726e-24 #[g]
kB         = 1.3806e-16 #erg/K
gamma      = 5./3.
HFrac      = 0.76
mH         = 1.6735575e-24
Gcode      = 6.672696243742867e-08/Vincgs**2./Lincgs*Mincgs
sfrden     = float(params['dPhysDenMin'])*mH/Mincgs*Lincgs**3.
tipsyMunit = float(params['dMsolUnit'][0:14])
Mmin       = float(params['dMinGasMass'])*tipsyMunit*1e-10*h

# WRITE
print("Writing Gizmo file...")
of     = h5py.File(folder+"/"+simname+'.000.hdf5', 'w')
groups = ['Header','PartType0','PartType1','PartType4'] #Type1 will be DM at all resolution levels. Type4 is stars and BH

# create groups
for g in groups:
    of.create_group(g)

# header
print("Header")
header = {'BoxSize':box,\
        'ComovingIntegrationOn':int(params['bComove']),\
        'Density_Threshold_For_SF_CodeUnits':sfrden,\
        #'Effective_Kernel_NeighborNumber':,\
        #'Effective_Kernel_NeighborNumber_CollisionlessParticles':,\
        #'Fixed_ForceSoftening_Keplerian_Kernel_Extent':,\
        'Flag_Cooling':cool,\
        #'Flag_DoublePrecision':,\
        'Flag_Feedback':int(params['bFeedBack'][0]),\
        'Flag_IC_Info':int(params['bStandard'][0]),\
        'Flag_Metals':int(params['nSpeciesExtraTrack'])+10,\
        'Flag_Sfr':int(params['bStarForm'][0]),\
        #'Flag_StellarAge':,\
        'GIZMO_version':-1,\
        'Gravitational_Constant_In_Code_Inits':Gcode,\
        'HubbleParam':h,\
        #'Kernel_Function_ID':,\
        'MassTable':array([0., 0., 0., 0., 0., 0.], dtype=float64),\
        #'Maximum_Mass_For_Cell_Split',\
        'Metals_Atomic_Number_Or_Key':array([ 0,  2,  6,  7,  8, 10, 11, 12, 13, 14, 16, 20, 21, 22, 23, 26], dtype=int32),\
        'Minimum_Mass_For_Cell_Merge':Mmin,\
        'NumFilesPerSnapshot':1,\
        'NumPart_ThisFile':array([ ng,  ndm, 0,       0,   ns,       0], dtype=int32),\
        'NumPart_Total':array([ ng,  ndm, 0,       0,   ns,       0], dtype=int32),\
        'NumPart_Total_HighWord':array([0, 0, 0, 0, 0, 0], dtype=uint32),\
        'Omega_Baryon':Ob,\
        'Omega_Lambda':float(params['dLambda']),\
        'Omega_Matter':float(params['dOmega0']),\
        'Omega_Radiation':0.0,\
        'Redshift':red,\
        #'SF_Timescale_At_Density_Threshold_CodeUnits':,\
        #'Solar_Abundances_Adopted':,\
        'Time':a,\
        'UnitLength_In_CGS':Lincgs,\
        'UnitMass_In_CGS':Mincgs,\
        'UnitVelocity_In_CGS':Vincgs}

for k in header.keys():
    v = header[k]
    y = type(array([1]))#trick for isinstance to work
    if isinstance(v,y):
        of['Header'].attrs.create(k,v,shape(v),v.dtype)
    else:
        of['Header'].attrs.create(k,v,shape(v),type(v))
if verbose:
    print("Header")
    print(list(of['Header'].attrs.items()))

# gas
print("Gas")
posg   = sim.gas['pos'].in_units('kpc a h**-1')+box/2. #comoving kpc, origin is now llc instead of centre
idg    = arange(0,ng)
velg   = array(sim.gas['vel'].in_units('km a**0.5 s**-1'))
epsg   = array(sim.gas['eps'].in_units('kpc a h**-1'))
massg  = array(sim.gas['mass'].in_units('Msol h**-1'))*1e-10 #Unit mass is 10^10 solar masses
temp   = array(sim.gas['temp']) #K
densg  = array(sim.gas['rho'].in_units('Msol kpc**-3 a**-3 h**2'))*1e-10 #Unit mass is 10^10 solar masses
metalg = array(sim.gas['metals']) # all but H and HeI, in mass fraction
HI     = 1.0 - array(sim.gas['HeI']) - metalg #in mass fraction, but here it is all H not just neutral
elec   = 1.0 - HI #H1 + elec must be one but Tipsy doesn't give separate fields. P2C needs it because it reads it in readgizmo, but it is only used to convert IE to T that is then not used
inte   = temp*(3.*HFrac + 1.+4.*HFrac*elec)*kB*Mincgs/(4.0*mproton*(gamma-1.)*Eincgs)
softg  = (3./(4.*pi)*massg/densg)**(1/3)/2. #get softening length from volume assuming sphere of radius 2h

of['PartType0'].create_dataset('Coordinates',data=posg,dtype='f8')
of['PartType0'].create_dataset('Masses',data=massg,dtype='f8')
of['PartType0'].create_dataset('ParticleIDs',data=idg,dtype='u8')
of['PartType0'].create_dataset('SmoothingLength',data=softg,dtype='f8')
of['PartType0'].create_dataset('Velocities',data=velg,dtype='f8')
of['PartType0'].create_dataset('AGS-Softening',data=epsg,dtype='f8')
of['PartType0'].create_dataset('Density',data=densg,dtype='f8')
of['PartType0'].create_dataset('Metallicity',data=metalg,dtype='f8')
of['PartType0'].create_dataset('ElectronAbundance',data=elec,dtype='f8')
of['PartType0'].create_dataset('InternalEnergy',data=inte,dtype='f8')
of['PartType0'].create_dataset('NeutralHydrogenAbundance',data=HI,dtype='f8') #FOR FEDE this has to be NHI/NH and at the moment it is NH/NTOT (probably)
of['PartType0'].create_dataset('Temperature',data=temp,dtype='f8')

# dm
print("DM")
posdm  = array(sim.dm['pos'].in_units('kpc a h**-1')) + box/2. #comoving kpc, origin is now llc instead of centre
iddm   = arange(ng,ng+ndm)
softdm = array(sim.dm['eps'].in_units('kpc a h**-1'))
veldm  = array(sim.dm['vel'].in_units('km a**0.5 s**-1'))
massdm = array(sim.dm['mass'].in_units('Msol h**-1'))*1e-10

of['PartType1'].create_dataset('Coordinates',data=posdm,dtype='f8')
of['PartType1'].create_dataset('Masses',data=massdm,dtype='f8')
of['PartType1'].create_dataset('ParticleIDs',data=iddm,dtype='u8')
of['PartType1'].create_dataset('SmoothingLength',data=softdm,dtype='f8') #is this correct??
of['PartType1'].create_dataset('Velocities',data=veldm,dtype='f8')

# stars
print("Stars")
poss   = array(sim.star['pos'].in_units('kpc'))*h/a + box/2.
ids    = arange(ng+ndm,ng+ndm+ns)
tform  = array(sim.star['tform'].in_units('a**-1 kpc s km**-1')) #not sure units are correct here 
vels   = array(sim.star['vel'].in_units('km a**0.5 s**-1'))
metals =  array(sim.star['metals']) #Hopefully this is mass fraction, so solar is 0.02
masss  = array(sim.star['mass'].in_units('Msol h**-1'))*1e-10
softs  = array(sim.star['eps'].in_units('kpc a h**-1'))

of['PartType4'].create_dataset('Coordinates',data=poss,dtype='f8')
of['PartType4'].create_dataset('Masses',data=masss,dtype='f8')
of['PartType4'].create_dataset('ParticleIDs',data=ids,dtype='u8')
of['PartType4'].create_dataset('SmoothingLength',data=softs,dtype='f8') #is this correct?
of['PartType4'].create_dataset('Velocities',data=vels,dtype='f8')
of['PartType4'].create_dataset('Metallicity',data=metals,dtype='f8')
of['PartType4'].create_dataset('StellarFormationTime',data=tform,dtype='f8')

of.close()
print("Done!")
