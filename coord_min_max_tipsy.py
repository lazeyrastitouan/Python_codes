from numpy import*
import pynbody as pn
import h5py, sys, os
import argparse
from astropy.io import fits
import pylab
from tipsyutils import readparams

# input parameters
parser     = argparse.ArgumentParser(description='Routine to find minimum and maximum coordinate of particles in a zoom-in to select a region. Args: < simulation (str) >, < Tipsy parameter file(str) >, < Desired P2Cregion side length in cMpc (float) >')
parser.add_argument('simulation',nargs=1,type=str, help='Simulation file')
parser.add_argument('paramfile',nargs=1,type=str, help='Tipsy parameter file')
parser.add_argument('desiredl',nargs=1,type=float, help='Desired P2C region length (cMpc)')
args       = parser.parse_args()
simulation = args.simulation[0]
paramfile  = args.paramfile[0]
desiredl   = args.desiredl[0]

sim     = pn.load(simulation)
a       = sim.properties['a']
params  = readparams(paramfile) #type str:str for the moment...
h       = sim.properties['h']
box     = float(params['dKpcUnit'])/1e3 #boxsize in cMpc
posorig = sim.gas['pos'] #code units
pos     = posorig.in_units('Mpc') #physical
pos     = (array(pos)/a+box/2.)*h #comoving, origin is now llc instead of centre
middle  = array([max((pos[:,i])+min(pos[:,i]))/2. for i in range(3)])

norm     = linalg.norm(pos,axis=1)
normmin  = argmin(norm)
normmax  = argmax(norm)
#llc      = posorig[normmin] #lower left corner
#urc      = posorig[normmax]
llc      = pos[normmin] #lower left corner
urc      = pos[normmax]


#cut_le = ((middle - desiredl/2.) - box/2.)/box #code units
#cut_re = ((middle + desiredl/2.) - box/2.)/box #code units
cut_le = middle - desiredl*h/2.
cut_re = middle + desiredl*h/2.

print("llc: "+str(llc))
print("urc: "+str(urc))
print("Box sides length: "+str(urc-llc))

print("For P2C:\n \
cut_le = "+str(cut_le)+"\n\
cut_re = "+str(cut_re))
