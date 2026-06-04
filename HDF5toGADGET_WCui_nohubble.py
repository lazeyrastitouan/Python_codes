import h5py
import numpy as np
import sys
import os


def write_head(fp, header, size):
	fp.write(np.int32(8))
	fp.write(header.encode())
	fp.write(np.int32(size+8))
	fp.write(np.int32(8))


def write_header(fl, attr,h):
	fl.write(np.int32(256))
	fl.write(np.array(attr['NumPart_ThisFile'][:6], dtype=np.int32))
	fl.write(np.array(attr['MassTable'][:6]/h, dtype=np.float64))
	fl.write(np.asarray([attr['Time'], attr['Redshift']], dtype=np.float64))
	fl.write(np.array([1, 1], dtype=np.int32))
	fl.write(np.array(attr['NumPart_Total'][:6], dtype=np.int32))
	fl.write(np.array([1, attr['NumFilesPerSnapshot']], dtype=np.int32))
	fl.write(np.asarray([attr['BoxSize']/h, attr['Omega0'], attr['OmegaLambda'], attr['HubbleParam']/h], dtype=np.float64))
	fl.write(np.zeros(np.int32(256/4)-6-12-4-2-6-2-8, dtype=np.int32))
	fl.write(np.int32(256))


if len(sys.argv) != 3:
	print('I require two arguments -- the input folder and the output folder')
	print('I exit, because you providing is not fit!', str(sys.argv))
	exit(0)
else:	
	folder = str(sys.argv[1])+"/"
	outpath = str(sys.argv[2])+"/"
print(folder, outpath)
h5files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

# loop to get all required information
for i in h5files:
	if i[-4:] == 'hdf5' or i[-4:] == 'HDF5':
		f = h5py.File(folder+i, "r")
	else:
		continue
	of = open(outpath+i[0:8]+'_noh'+i[8:-5], 'wb')
	attrs = f['/Header'].attrs
	h = attrs['HubbleParam']
	TNp = attrs['NumPart_ThisFile'][:6]
	TNs = np.sum(TNp)
	
	# HEAD
	write_head(of, "HEAD", 256)
	write_header(of, attrs, h)
	
	# POS
	write_head(of, "POS ", np.uint32(TNs*4*3))
	of.write(np.uint32(TNs*4*3))
	if 'PartType0' in f.keys():
		of.write(np.float32(f['/PartType0/Coordinates'].value)/h)
	if 'PartType1' in f.keys():
		of.write(np.float32(f['/PartType1/Coordinates'].value)/h)
	if 'PartType2' in f.keys():
		of.write(np.float32(f['/PartType2/Coordinates'].value)/h)
	if 'PartType3' in f.keys():
		of.write(np.float32(f['/PartType3/Coordinates'].value)/h)
	if 'PartType4' in f.keys():
		of.write(np.float32(f['/PartType4/Coordinates'].value)/h)
	if 'PartType5' in f.keys():
		of.write(np.float32(f['/PartType5/Coordinates'].value)/h)
	of.write(np.uint32(TNs*4*3))
	
	# VEL
	write_head(of, "VEL ", np.uint32(TNs*4*3))
	of.write(np.uint32(TNs*4*3))
	if 'PartType0' in f.keys():
		of.write(np.float32(f['/PartType0/Velocities'].value))
	if 'PartType1' in f.keys():
		of.write(np.float32(f['/PartType1/Velocities'].value))
	if 'PartType2' in f.keys():
		of.write(np.float32(f['/PartType2/Velocities'].value))
	if 'PartType3' in f.keys():
		of.write(np.float32(f['/PartType3/Velocities'].value))
	if 'PartType4' in f.keys():
		of.write(np.float32(f['/PartType4/Velocities'].value))
	if 'PartType5' in f.keys():
		of.write(np.float32(f['/PartType5/Velocities'].value))
	of.write(np.uint32(TNs*4*3))
	
	# ID
	write_head(of, "ID  ", np.uint32(TNs*4))
	of.write(np.uint32(TNs*4))
	if 'PartType0' in f.keys():
		of.write(np.uint32(f['/PartType0/ParticleIDs'].value))
	if 'PartType1' in f.keys():
		of.write(np.uint32(f['/PartType1/ParticleIDs'].value))
	if 'PartType2' in f.keys():
		of.write(np.uint32(f['/PartType2/ParticleIDs'].value))
	if 'PartType3' in f.keys():
		of.write(np.uint32(f['/PartType3/ParticleIDs'].value))
	if 'PartType4' in f.keys():
		of.write(np.uint32(f['/PartType4/ParticleIDs'].value))
	if 'PartType5' in f.keys():
		of.write(np.uint32(f['/PartType5/ParticleIDs'].value))
	of.write(np.uint32(TNs*4))
	
	# MASS
	if 'Masses' in f['/PartType1'].keys(): #and 'Masses' in f['/PartType2'].keys():
		write_head(of, "MASS", np.uint32(TNs*4))
		of.write(np.uint32(TNs*4))
		if 'PartType0' in f.keys():
			of.write(np.float32(f['/PartType0/Masses'].value)/h)
		if 'PartType1' in f.keys():
			of.write(np.float32(f['/PartType1/Masses'].value)/h)
		if 'PartType2' in f.keys():
			of.write(np.float32(f['/PartType2/Masses'].value)/h)
		if 'PartType3' in f.keys():
			of.write(np.float32(f['/PartType3/Masses'].value)/h)
		if 'PartType4' in f.keys():
			of.write(np.float32(f['/PartType4/Masses'].value)/h)
		if 'PartType5' in f.keys():
			of.write(np.float32(f['/PartType5/Masses'].value)/h)
		of.write(np.uint32(TNs*4))
	
	# U <only gas>
	if 'PartType0' in f.keys():
		write_head(of, "U   ", np.uint32(TNp[0]*4))
		of.write(np.uint32(TNp[0]*4))
		of.write(np.float32(f['/PartType0/InternalEnergy'].value))
		of.write(np.uint32(TNp[0]*4))
	
	# AGE <only Star> stellar formation time !!No BHs
	if 'PartType4' in f.keys():
		write_head(of, "AGE ", np.uint32(TNp[4]*4))
		of.write(np.uint32(TNp[4]*4))
		of.write(np.float32(f['/PartType4/StarFormationTime'].value))
		of.write(np.uint32(TNp[4]*4))
	
	if 'PartType0' in f.keys(): 
		# Z <gas + Star>
		if 'PartType4' in f.keys():
			Tbp = TNp[0]+TNp[4]
		else:
			Tbp = TNp[0]
		write_head(of, "Z   ", np.uint32(Tbp*4))
		of.write(np.uint32(Tbp*4))
		of.write(np.float32(f['/PartType0/Metallicity'].value))
		if 'PartType4' in f.keys():
			of.write(np.float32(f['/PartType4/Metallicity'].value))
		of.write(np.uint32(Tbp*4))
		
		# Rho <only gas>
		write_head(of, "RHO ", np.uint32(TNp[0]*4))
		of.write(np.uint32(TNp[0]*4))
		of.write(np.float32(f['/PartType0/Density'].value))
		of.write(np.uint32(TNp[0]*4))
		
		# NE <only gas>
		write_head(of, "NE  ", np.uint32(TNp[0]*4))
		of.write(np.uint32(TNp[0]*4))
		of.write(np.float32(f['/PartType0/ElectronAbundance'].value))
		of.write(np.uint32(TNp[0]*4))
		
		# NH <only gas>
		write_head(of, "NH  ", np.uint32(TNp[0]*4))
		of.write(np.uint32(TNp[0]*4))
		of.write(np.float32(f['/PartType0/NeutralHydrogenAbundance'].value))
		of.write(np.uint32(TNp[0]*4))
		
		# HSML <only gas>  !!Note no hsml for Arepo
		#Calculate
		dens   = np.array(f[u'PartType0/Density'])
		mass   = np.array(f[u'PartType0/Masses'])
		volume = np.divide(mass, dens)
		hsml   = np.power(3.*volume/(4.*np.pi), 1./3.)
		
		write_head(of, "HSML", np.uint32(TNp[0]*4))
		of.write(np.uint32(TNp[0]*4))
		#of.write(np.float32(f['/PartType0/AllowRefinement'].value))
		of.write(hsml)
		of.write(np.uint32(TNp[0]*4))
		
		# SFR <only gas>
		write_head(of, "SFR ", np.uint32(TNp[0]*4))
		of.write(np.uint32(TNp[0]*4))
		of.write(np.float32(f['/PartType0/StarFormationRate'].value))
		of.write(np.uint32(TNp[0]*4))
		
