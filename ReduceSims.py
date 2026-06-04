from numpy import*
import sys
import os
import h5py
import argparse
from multiprocessing import Pool

def write_head(fp, header, size):
        fp.write(int32(8))
        fp.write(header.encode())
        fp.write(int32(size+8))
        fp.write(int32(8))


def write_header(fl, attr, newparams):
        fl.write(int32(256))
        fl.write(newparams['NumPart_ThisFile'][:6])
        fl.write(array(attr['MassTable'][:6], dtype=float64))
        fl.write(asarray([attr['Time'], attr['Redshift']], dtype=float64))
        fl.write(array([0, 0], dtype=int32))
        fl.write(newparams['NumPart_Total'][:6])
        fl.write(array([1, attr['NumFilesPerSnapshot']], dtype=int32))
        fl.write(asarray([newparams['BoxSize'], newparams['Omega0'], attr['OmegaLambda'], attr['HubbleParam']], dtype=float64))
        fl.write(zeros(int32(256/4)-6-12-4-2-6-2-8, dtype=int32))
        fl.write(int32(256))

# What to do
parser = argparse.ArgumentParser(description='Select particle of a given simulation region only. Args: < path to snapshot (hdf5) (str) >, < output folder (str) >, < xmin, ymin, zmin (floats) >, < new box size (float) >, < redshift >')
parser.add_argument('snappath',nargs=1,type=str, help='Path to snapshot (hdf5)')
parser.add_argument('folderout',nargs=1,type=str, help='Output folder (str)')
parser.add_argument('xyz',nargs=3,type=float, help='xmin, ymin, zmin [kpc/h]')
parser.add_argument('size',nargs=1,type=float, help='new box size [kpc/h]')
parser.add_argument('redshift',nargs=1,type=float, help='redshift')
args     = parser.parse_args()
snappath = args.snappath[0]
outpath  = args.folderout[0]
xmin 	 = args.xyz[0]
ymin 	 = args.xyz[1]
zmin 	 = args.xyz[2]
newsize  = args.size[0]
red      = args.redshift[0]
xmax     = xmin+newsize
ymax     = ymin+newsize
zmax     = zmin+newsize

# settings
datatype = 'Gizmo' #'Arepo' 'Gadget2'
add_b    = False #Need Arepo for true
verbose  = False
hydro    = True #will reduce and save only baryons
if hydro:
	ptype = 'PartType0'
else:
	ptype = 'PartType1'
		

print("Datatye = "+datatype)
print("Hydro sim = "+str(hydro)) 
print("Adding baryons = "+str(add_b))

h5files  = [f for f in os.listdir(snappath+"/") if os.path.isfile(os.path.join(snappath+"/", f)) and f[-4:] == 'hdf5']
h5files  = sorted(h5files)
h5files  = sorted(h5files,key=len)
if verbose:
	print(h5files)

print("Keeping particule in cube\n"+str(xmin/1000.)+" < x < "+str(xmax/1000.)+" Mpc/h\n"+str(ymin/1000.)+" < y < "+str(ymax/1000.)+" Mpc/h\n"+str(zmin/1000.)+" < z < "+str(zmax/1000.)+" Mpc/h")
 
print("Selecting particles in volume (only considering highest refinment ones)...")
numfiles = len(h5files)
npthisfile = []
ids, pos, vel, mas = [[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)]
if hydro:
	if datatype == 'Arepo':
		AlR, CR, Den, Xe, HRM, IE, NH, Prs, Sof, SFR = [[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)]
	elif datatype == 'Gizmo':
		AGS, Den, Xe, IE, NH, PcID, PIGN, SL, SFR, Temp = [[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)],[[] for x in range(numfiles)]

np = 0
for i in h5files:
	f     = h5py.File(snappath+"/"+i, "r")
	fi    = h5files.index(i)
	print(f)
	
	if i == 0 and 'PartType0' in f.keys():
		hydro = True

	idst  = array(f[ptype]['ParticleIDs'])
	post  = array(f[ptype]['Coordinates'])
	velt  = array(f[ptype]['Velocities'])
	if 'Masses' in f[ptype].keys():
		mast = array(f[ptype]['Masses'])
		
	if hydro:
		if datatype == 'Arepo':
			AlRt = array(f[ptype]['AllowRefinement'])
			CRt  = array(f[ptype]['CoolingRate'])
			Dent = array(f[ptype]['Density'])
			Xet  = array(f[ptype]['ElectronAbundance'])
			HRMt = array(f[ptype]['HighResGasMass'])
			IEt  = array(f[ptype]['InternalEnergy'])
			NHt  = array(f[ptype]['NeutralHydrogenAbundance'])
			Prst = array(f[ptype]['Pressure'])
			Soft = array(f[ptype]['Softenings'])
			SFRt = array(f[ptype]['StarFormationRate'])
		elif datatype == 'Gizmo':
			AGSt  = array(f[ptype]['AGS-Softening'])
			Dent  = array(f[ptype]['Density'])
			Xet   = array(f[ptype]['ElectronAbundance'])
			IEt   = array(f[ptype]['InternalEnergy']) 
			NHt   = array(f[ptype]['NeutralHydrogenAbundance'])
			PcIDt = array(f[ptype]['ParticleChildIDsNumber'])
			PIGNt = array(f[ptype]['ParticleIDGenerationNumber'])
			SLt   = array(f[ptype]['SmoothingLength'])
			SFRt  = array(f[ptype]['StarFormationRate'])
			Tempt = array(f[ptype]['Temperature'])

	if hydro:	
		numparts = f['Header'].attrs['NumPart_ThisFile'][0]
	else:
		numparts = f['Header'].attrs['NumPart_ThisFile'][1]	
	for j in range(numparts):
		if xmin <= post[j,0] <= xmax and ymin <= post[j,1] <= ymax and zmin <= post[j,2] <= zmax:		
			ids[fi].append(idst[j])
			#Translate pos to be consistent with newsize
			pos[fi].append(array([post[j,0]-xmin,post[j,1]-ymin,post[j,2]-zmin]))
			vel[fi].append(velt[j])
			if 'Masses' in f[ptype].keys():
				mas[fi].append(mast[j])
	if hydro:
		if datatype == 'Arepo':
			for j in range(numparts):
				if xmin <= post[j,0] <= xmax and ymin <= post[j,1] <= ymax and zmin <= post[j,2] <= zmax:
					AlR[fi].append(AlRt[j])
					CR[fi].append(CRt[j])
					Den[fi].append(Dent[j])
					Xe[fi].append(Xet[j])
					HRM[fi].append(HRMt[j])
					IE[fi].append(IEt[j])
					NH[fi].append(NHt[j])
					Prs[fi].append(Prst[j])
					Sof[fi].append(Soft[j]) 
					SFR[fi].append(SFRt[j])
		elif datatype == 'Gizmo':
			for j in range(numparts):
				if xmin <= post[j,0] <= xmax and ymin <= post[j,1] <= ymax and zmin <= post[j,2] <= zmax:
					AGS[fi].append(AGSt[j])
					Den[fi].append(Dent[j])
					Xe[fi].append(Xet[j])
					IE[fi].append(IEt[j])
					NH[fi].append(NHt[j])
					PcID[fi].append(PcIDt[j])				
					PIGN[fi].append(PIGNt[j])
					SL[fi].append(SLt[j])
					SFR[fi].append(SFRt[j])
					Temp[fi].append(Tempt[j])
	print(len(ids[fi]))
	npthisfile.append(len(ids[fi]))
	np += len(ids[fi])
	
print("Done. Kept "+str(np)+" particles in total")

#adjusting cosmology (only for non hydro runs)
if hydro == False:
	ft     = h5py.File(snappath+"/"+h5files[0], "r") 
	mp     = f['/Header'].attrs['MassTable'][1]*1e10
	Omega0 = f['/Header'].attrs['Omega0']
	rho_c  = 27.75e10 #h^2 Msun/Mpc^3
	rho_m  = np*mp/(newsize/1000.)**3. 
	Omegam = rho_m/rho_c
	print("Old Omega_m="+str(Omega0)+"\nNew Omega_m="+str(Omegam))

#transforming lists to arrays shape=((numfiles,npthifile))
ids = array(ids)
pos = array(pos)
vel = array(vel)
mas = array(mas)
if hydro:
	Den = array(Den)
	Xe  = array(Xe)
	IE  = array(IE)
	NH  = array(NH)
	SFR = array(SFR)
	if datatype == 'Arepo':
		AlR = array(AlR)
		CR  = array(CR)
		HRM = array(HRM)
		Prs = array(Prs)
		Sof = array(Sof)
	elif datatype == 'Gizmo':
		AGS  = array(AGS)
		PcID = array(PcID)
		PIGN = array(PIGN)
		SL   = array(SL)
		Temp = array(Temp)

for i in range(numfiles):
	ids[i] = array(ids[i])
	pos[i] = array(pos[i])
	vel[i] = array(vel[i])
	mas[i] = array(mas[i])
	if hydro:
		Den[i] = array(Den[i])
		Xe[i]  = array(Xe[i])
		IE[i]  = array(IE[i])
		NH[i]  = array(NH[i])
		SFR[i] = array(SFR[i])
		if datatype == 'Arepo': 
			AlR[i] = array(AlR[i])
			CR[i]  = array(CR[i])
			HRM[i] = array(HRM[i])
			Prs[i] = array(Prs[i])
			Sof[i] = array(Sof[i])
		elif datatype == 'Gizmo':
			AGS[i]  = array(AGS[i])
			PcID[i] = array(PcID[i])
			PIGN[i] = array(PIGN[i])
			SL[i]   = array(SL[i])
			Temp[i] = array(Temp[i])

#one3 = pos/pos #array of ones with correct irregular shape for baryons arrays
one1 = (ids+1)/(ids+1) #array of ones with correct irregular shape for baryons arrays

if add_b:
	print("Creating baryons...")
	normv   = [[] for x in range(numfiles)]
	for i in range(numfiles):
		for j in range(len(vel[i])):
			normv[i].append(linalg.norm(vel[i][j]))
		normv[i]=array(normv[i])
	normv   = array(normv)
	Omegab  = Omega0/6.362444 #hard coded ratio to keep ratio of the base sim
	mb      = mp/6.362444
	masb    = mb*1e-10*one1
	mc      = mp-mb
	MassTab = [mb*1e-10,mc*1e-10,0.0,0.0,0.0,0.0]
	H       = newsize/np**(1./3.) #mean interdistance particle
	volcell = 4./3.*pi*H**3. #vol of cell
	IE      = 0.5*normv**2./(1.+red) #need physical to be consistent with arepo output 
	soft    = H/40.*one1
	dens    = mb/volcell*one1
	pres    = mb*normv**2./volcell
	SFR     = one1*0.0 #no SFR in fake baryons
	CR      = one1*0.0 #no CR in fake baryons
	ElecAb  = 1.157*one1 #comes from readgadget 
	NHAb	= 0.76*one1 #cosmic mean based on Planck

print("Writing files...")
if datatype == 'Gadget2': #DM only
	for i in range(numfiles):
		print(str(i+1)+"/"+str(numfiles))
	
		of        = open(outpath+h5files[i][:-5], 'wb')
		f         = h5py.File(snappath+"/"+h5files[i], "r")
		attrs     = f['/Header'].attrs
		newparams = {'BoxSize' : newsize, 'Omega0' : Omegam, 'NumPart_ThisFile' : array([0,npthisfile[i],0,0,0,0],dtype=int32) , 'NumPart_Total' : array([0,np,0,0,0,0],dtype=int32)}
	
		# HEAD
		write_head(of, "HEAD", 256)
		write_header(of, attrs, newparams)
		# POS
		write_head(of, "POS ", uint32(npthisfile[i]*4*3))
		of.write(uint32(npthisfile[i]*4*3))
		of.write(float32(pos[i]))
		of.write(uint32(npthisfile[i]*4*3))
		# VEL
		write_head(of, "VEL ", uint32(npthisfile[i]*4*3))
		of.write(uint32(npthisfile[i]*4*3))
		of.write(float32(vel[i]))
		of.write(uint32(npthisfile[i]*4*3))
		# ID
		write_head(of, "ID  ", uint32(npthisfile[i]*4)) #*8 for longids
		of.write(uint32(npthisfile[i]*4)) #*8 for longids
		of.write(uint32(ids[i]))
		of.write(uint32(npthisfile[i]*4)) #*8 for longids	
		# MASS
		if 'Masses' in f['PartType1'].keys():
			if verbose:
				print("Also writing mass table")
			write_head(of, "MASS", uint32(npthisfile[i]*4))
			of.write(uint32(npthisfile[i]*4))
			of.write(float32(mas[i]))
			of.write(uint32(npthisfile[i]*4))
	
		of.close()	

#Arepo
elif datatype == 'Arepo':
	if add_b != hydro : #xor
		groups = ['Config','Header','Parameters','PartType0','PartType1']
	else:
		groups = ['Config','Header','Parameters','PartType1']
	
	for i in range(numfiles):
		print(str(i+1)+"/"+str(numfiles))

		of        = h5py.File(outpath+h5files[i][:-5]+'_reduced.hdf5', 'w')
		f         = h5py.File(snappath+"/"+h5files[i], "r")
		attrs     = f['/Header'].attrs
		if add_b:
			newparams = {'BoxSize' : newsize, 'Omega0' : Omegam, 'OmegaBaryon' : Omegab, 'MassTable' : array(MassTab,dtype=float32), 'NumPart_ThisFile' : array([npthisfile[i],npthisfile[i],0,0,0,0],dtype=int32) , 'NumPart_Total' : array([np,np,0,0,0,0],dtype=int32)} 
		elif hydro:
			newparams = {'BoxSize' : newsize, 'MassTable' : array([0.,0.,0.,0.,0.,0.],dtype=float32), 'NumPart_ThisFile' : array([npthisfile[i],0,0,0,0,0],dtype=int32) , 'NumPart_Total' : array([np,0,0,0,0,0],dtype=int32)} 
		else:
			newparams = {'BoxSize' : newsize, 'Omega0' : Omegam, 'NumPart_ThisFile' : array([0,npthisfile[i],0,0,0,0],dtype=int32) , 'NumPart_Total' : array([0,np,0,0,0,0],dtype=int32)}

		# create groups
		for g in groups:
			of.create_group(g)
		
		#CONFIG
		for k in f['Config'].attrs.keys():
			v = f['Config'].attrs[k]
			t = type(v)
			if v == '':
				of['Config'].attrs.create(k,'',(),'S10')
			else:
				of['Config'].attrs.create(k,v,(),t)	
		
		if add_b:
			of['Config'].attrs.create('GENERATE_GAS_IN_ICS','',(),'S10')	
		
		if verbose:
			print("\n Config")
			print(of['Config'].attrs.items())		
	
		#HEADER
		for k in f['Header'].attrs.keys():
			v=f['Header'].attrs[k]
			if k in newparams.keys():
				of['Header'].attrs.create(k,newparams[k],shape(v),v.dtype)
			else:
				if type(v)==string_:
					of['Header'].attrs.create(k,v,(),'S10')
				else:
					of['Header'].attrs.create(k,v,shape(v),v.dtype)	
		if verbose:
			print("\n Header")
			print(of['Header'].attrs.items())	

		#PARAMETERS - not corrected from the base sim 
		for k in f['Parameters'].attrs.keys():
			v = f['Parameters'].attrs[k]
			if type(v) == string_:
				of['Parameters'].attrs.create(k,v,(),'S10')
			else:
				of['Parameters'].attrs.create(k,v,shape(v),v.dtype)
		if verbose:
			print("\n Parameters")
			print(of['Parameters'].attrs.items())

		#PARTICLES
		#datasets for DM and hydro		
		of[ptype].create_dataset('Coordinates',data=pos[i],dtype='f8')
		if 'Masses' in f[ptype].keys():
			of[ptype].create_dataset('Masses',data=mas[i],dtype='f8')
		of[ptype].create_dataset('ParticleIDs',data=ids[i],dtype='u8')
		of[ptype].create_dataset('Velocities',data=vel[i],dtype='f8')
		if hydro:
			of[ptype].create_dataset('AllowRefinement',data=AlR[i],dtype='f8')
			of[ptype].create_dataset('CoolingRate',data=CR[i],dtype='f8')
			of[ptype].create_dataset('Density',data=Den[i],dtype='f8')
			of[ptype].create_dataset('ElectronAbundance',data=Xe[i],dtype='f8')
			of[ptype].create_dataset('HighResGasMass',data=HRM[i],dtype='f8')
			of[ptype].create_dataset('InternalEnergy',data=IE[i],dtype='f8')
			of[ptype].create_dataset('NeutralHydrogenAbundance',data=NH[i],dtype='f8')
			of[ptype].create_dataset('Pressure',data=Prs[i],dtype='f8')
			of[ptype].create_dataset('Softenings',data=Sof[i],dtype='f8')
			of[ptype].create_dataset('StarFormationRate',data=SFR[i],dtype='f8')		
			
		#datasets for fake baryons
		# same pos, vel, and IDs as dm ; IE=vel^2, dens=mb/vol cell , soft = mean dist baryons, pressure = vol cell*mb*vel**2. 
		if add_b:
			of['PartType0'].create_dataset('CoolingRate',data=CR[i],dtype='f8')
			of['PartType0'].create_dataset('Coordinates',data=pos[i],dtype='f8')	
			of['PartType0'].create_dataset('Density',data=dens[i],dtype='f8')
			of['PartType0'].create_dataset('ElectronAbundance',data=ElecAb[i],dtype='f8')
			of['PartType0'].create_dataset('InternalEnergy',data=IE[i],dtype='f8') 
			of['PartType0'].create_dataset('Masses',data=masb[i],dtype='f8')
			of['PartType0'].create_dataset('NeutralHydrogenAbundance',data=NHAb[i],dtype='f8')
			of['PartType0'].create_dataset('ParticleIDs',data=ids[i],dtype='u8')
			of['PartType0'].create_dataset('Pressure',data=pres[i],dtype='f8')
			of['PartType0'].create_dataset('Softenings',data=soft[i],dtype='f8')
			of['PartType0'].create_dataset('StarFormationRate',data=SFR[i],dtype='f8')
			of['PartType0'].create_dataset('Velocities',data=vel[i],dtype='f8')
		
		of.close()
#Gizmo
else:
	if hydro:
		groups = ['Header','PartType0']
	else:
		groups = ['Header','PartType1']

	for i in range(numfiles):
		print(str(i+1)+"/"+str(numfiles))

		of        = h5py.File(outpath+"/"+h5files[i][:-5]+'_reduced.hdf5', 'w')
		f         = h5py.File(snappath+"/"+h5files[i], "r")
		attrs     = f['/Header'].attrs
		if hydro:
			newparams = {'BoxSize' : newsize, 'MassTable' : array([0.,0.,0.,0.,0.,0.],dtype=float32), 'NumPart_ThisFile' : array([npthisfile[i],0,0,0,0,0],dtype=int32) , 'NumPart_Total' : array([np,0,0,0,0,0],dtype=uint32)}
		else:
			newparams = {'BoxSize' : newsize, 'Omega0' : Omegam, 'NumPart_ThisFile' : array([0,npthisfile[i],0,0,0,0],dtype=int32) , 'NumPart_Total' : array([0,np,0,0,0,0],dtype=int32)}

		# create groups
		for g in groups:
			of.create_group(g)

		#HEADER
		for k in f['Header'].attrs.keys():
			v=f['Header'].attrs[k]
			if k in newparams.keys():
				of['Header'].attrs.create(k,newparams[k],shape(v),v.dtype)
			else:
				if type(v)==string_:
					of['Header'].attrs.create(k,v,(),'S10')
				else:
					of['Header'].attrs.create(k,v,shape(v),v.dtype)
		if verbose:
			print("\n Header")
			print(list(of['Header'].attrs.items()))

		#PARTICLES
		#datasets for DM and hydro              
		of[ptype].create_dataset('Coordinates',data=pos[i],dtype='f8')
		if 'Masses' in f[ptype].keys():
			of[ptype].create_dataset('Masses',data=mas[i],dtype='f8')
		of[ptype].create_dataset('ParticleIDs',data=ids[i],dtype='u8')
		of[ptype].create_dataset('Velocities',data=vel[i],dtype='f8')
		if hydro:
			of[ptype].create_dataset('AGS-Softening',data=AGS[i],dtype='f8')
			of[ptype].create_dataset('Density',data=Den[i],dtype='f8')
			of[ptype].create_dataset('ElectronAbundance',data=Xe[i],dtype='f8')
			of[ptype].create_dataset('InternalEnergy',data=IE[i],dtype='f8')
			of[ptype].create_dataset('NeutralHydrogenAbundance',data=NH[i],dtype='f8')
			of[ptype].create_dataset('ParticleChildIDsNumber',data=PcID[i],dtype='f8')
			of[ptype].create_dataset('ParticleIDGenerationNumber',data=PIGN[i],dtype='f8')
			of[ptype].create_dataset('SmoothingLength',data=SL[i],dtype='f8')
			of[ptype].create_dataset('StarFormationRate',data=SFR[i],dtype='f8')
			of[ptype].create_dataset('Temperature',data=Temp[i],dtype='f8')	
		
		of.close()
print("Done!")
