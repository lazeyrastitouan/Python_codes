from numpy import*
import matplotlib.pyplot as plt
import argparse, sys, os, h5py

# What to do
parser = argparse.ArgumentParser(description='Find halos above a mass threshold and returns mass, idx, pos, radius. Args: < path to catalogs (hdf5) (str) >, < log10(lower Mass)  (float) >, < log10(upper Mass)  (float) >')
parser.add_argument('path',nargs=1,type=str, help='Path to halo catalogs')
parser.add_argument('Mlow',nargs=1,type=float, help='Log lower mass threshold [Msun]')
parser.add_argument('Mhigh',nargs=1,type=float, help='Log higher mass threshold [Msun]')
args      = parser.parse_args()
path      = args.path[0]
Mlow	  = args.Mlow[0]
Mhigh     = args.Mhigh[0]
simtype   = 'Gizmo' #'Arepo'
DMonly    = True

#load
M, R, pos, ID, FS  = [],[],[],[],[]
h5files            = [f for f in os.listdir(path+"/") if os.path.isfile(os.path.join(path+"/", f)) and f[-4:] == 'hdf5']
h5files            = sorted(h5files)
h5files            = sorted(h5files,key=len)
for i in h5files:
	f = h5py.File(path+"/"+i, "r")
	print(f)
	Mt   = array(f['/Group/GroupMass'])
	post = array(f['/Group/GroupPos'])
	FSt  = array(f['Group/GroupFirstSub'])
	lm = len(Mt)
	if lm == 0:
		continue
	if simtype == 'Arepo':
		Rt  = array(f['Group/Group_R_Mean200'])
	else:
		Rt   = array(f['Group/Group_vsDelta_Radius'])
		#IDst = array(f['IDs/PID '])
		#if len(IDst) == 0:
		#	continue
		if DMonly:
			IDt = array(f['Subhalo/SubhaloIDMostbound'])
		else:
			if lm == 1:
				#IDt = array([f['IDs/PID '][1]])
				IDt = array([f['Subhalo/SubhaloIDMostbound'][0]])	
			else:
				#IDt = array(f['IDs/PID '][:,1])
				IDt = array(f['Subhalo/SubhaloIDMostbound'][:,0]) #to have the ID of the most bound particle of the first subhalo, which is the main halo center
	for m in range(lm):
		M.append(Mt[m]*1e10)
		pos.append(post[m])
		R.append(Rt[m])
		FS.append(FSt[m])
	nsub_thisfile = len(IDt)
	for n in range(nsub_thisfile):
		ID.append(IDt[n])

M   = array(M)
R   = array(R)
pos = array(pos)
ID  = array(ID)
FS  = array(FS)

redshift = round(f['Header'].attrs['Redshift'],1)
h        = f['Header'].attrs['HubbleParam']
M        = log10(M/h)

#select halos and write to a file
idsup  = where(M >= Mlow)
idslow = where(M <= Mhigh)
ids    = nonzero(in1d(idsup[0], idslow[0]))[0]
print("Kept "+str(len(ids))+" halos")

f=open(path+"/halos_Mlow"+str(Mlow)+"_Mhigh"+str(Mhigh)+".txt","w")
for idx in ids:
	f.write(str(M[idx])+"\t"+str(idx)+"\t"+str(pos[idx,0])+"\t"+str(pos[idx,1])+"\t"+str(pos[idx,2])+"\t"+str(R[idx])+"\t"+str(ID[FS[idx]])+"\n")
f.close()
