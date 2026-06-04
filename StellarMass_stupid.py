from numpy import*
import sys
import os
import h5py
import argparse
#from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# What to do
parser = argparse.ArgumentParser(description='Compute and plot the stellar mass to halo mass relation. Args: < path to snapshots (hdf5, str) >, < path to groups (hdf5, str) >')
parser.add_argument('snap',nargs=1,type=str, help='Path to snapshot files (hdf5)')
parser.add_argument('group',nargs=1,type=str, help='Path to group files (hdf5)')
args      = parser.parse_args()
snappath  = args.snap[0]
grouppath = args.group[0]

verbose = 0
hubble  = 0.677

h5snaps  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5']
h5snaps  = sorted(h5snaps)
h5snaps  = sorted(h5snaps,key=len)
h5groups = [f for f in os.listdir(grouppath+"/") if os.path.isfile(os.path.join(grouppath+"/", f)) and f[-4:] == 'hdf5']
h5groups = sorted(h5groups)
h5groups = sorted(h5groups,key=len)
if verbose:
	print("Snapshots "+str(h5snaps))
	print("Groups "+str(h5groups))

fg0    = h5py.File(grouppath+'/'+h5groups[0],'r')
fs0    = h5py.File(snappath+'/'+h5snaps[0],'r')
attrsg = fg0['Header'].attrs
attrss = fs0['Header'].attrs

Nh     = attrsg['Ngroups_Total']
red    = round(attrsg['Redshift'],1)
out    = 'logMh_logMstar_z'+str(red)+'.dat'

grouppos, R = [],[]
starpos     = []
Mh          = []
Mstar       = zeros(Nh)

#loop over group files to get positions, radii and masses of groups
print("Loading halos...")
for i in h5groups:
	print("File "+str(h5groups.index(i)+1)+"/"+str(len(h5groups)))
	
	fg   = h5py.File(grouppath+'/'+i,'r')
	Nhtf = fg['Header'].attrs['Ngroups_ThisFile']

	for h in range(Nhtf):
		Mh.append(fg['Group/Group_M_Mean200'][h]*1e10/hubble)
		R.append(fg['Group/Group_R_Mean200'][h])
		grouppos.append(fg['Group/GroupPos'][h])

Mh       = array(Mh)
R        = array(R)
grouppos = array(grouppos)

#loop over snapshot files to load positions and masses of stellar particules
print("Loading and selecting stars...")
for j in h5snaps:
	print("File "+str(h5snaps.index(j)+1)+"/"+str(len(h5snaps)))
	
	fs   = h5py.File(snappath+'/'+j,'r')
	Nstf = fs['Header'].attrs['NumPart_ThisFile'][4]

	for s in range(Nstf):
		starpos.append(fs['PartType4/Coordinates'][s])
	
		#test wether stellar particule belongs to a group and add its mass to the stellar mass of that group
		for h in range(Nh):
			normdiff = linalg.norm(grouppos[h]-starpos[s])
			if normdiff <= R[h]:
				Mstar[h]+=(fs['PartType4/Masses'][s]*1e10/hubble)

print("Done!")

#write
#fout = open(out,'W')
#for i in range(Nh):
#	fout.write(str(Mh[i])+"\t"+str(Mstar[i])+"\n") 
#fout.close()

#plot Mstar vs Mh
plt.plot(log10(Mh),log10(Mstar),'r.')
plt.xlabel(r'$\log(M_{200{\textrm mean}}) \; [M_\odot]$',fontsize=25,labelpad=15)
plt.ylabel(r'$\log(M_\star) \; [M_\odot]$',fontsize=25,labelpad=15)
plt.title(r'$z = $'+str(red),fontsize=20)
#plt.savefig(bbox_inches='tight')
plt.show()
