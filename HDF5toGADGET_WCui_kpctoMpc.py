# Last Changelog
# AKh 18.04.18 all units are converted to kpc/h
# read data, x1000,convert to float, this works only for Hestia-Hydro-double snapshots
import h5py
import numpy as np
import sys
import os


def write_head(fp, header, size):
	fp.write(np.int32(8))
	fp.write(header.encode())
	fp.write(np.int32(size+8))
	fp.write(np.int32(8))


def write_header(fl, attr):
	fl.write(np.int32(256))
	fl.write(np.array(attr['NumPart_ThisFile'][:6], dtype=np.int32))
	fl.write(np.array(attr['MassTable'][:6], dtype=np.float64))
	fl.write(np.asarray([attr['Time'], attr['Redshift']], dtype=np.float64))
	fl.write(np.array([1, 1], dtype=np.int32))
	fl.write(np.array(attr['NumPart_Total'][:6], dtype=np.int32))
	fl.write(np.array([1, attr['NumFilesPerSnapshot']], dtype=np.int32))
	fl.write(np.asarray([attr['BoxSize']/1000.0, attr['Omega_Matter'], attr['Omega_Lambda'], attr['HubbleParam']], dtype=np.float64))
	fl.write(np.zeros(np.int32(256/4)-6-12-4-2-6-2-8, dtype=np.int32))
	fl.write(np.int32(256))


if len(sys.argv) != 2:
	print('I require only one argument -- the folder name!')
	print('I exit, because you providing is not fit!', str(sys.argv))
	exit(0)
else:
	folder = str(sys.argv[1])+"/"
print(folder, sys.argv[1])
h5files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

# loop to get all required information
for i in h5files:
	if i[-4:] == 'hdf5' or i[-4:] == 'HDF5':
		f = h5py.File(folder+i, "r")
		#of = open(folder+i[:-4]+"G3", 'wb')
		of = open(folder+i[:-5], 'wb')
		attrs = f['/Header'].attrs
		TNp = attrs['NumPart_ThisFile'][:6]
		TNs = np.sum(TNp)

		# HEAD
		write_head(of, "HEAD", 256)
		write_header(of, attrs)

		# POS
		write_head(of, "POS ", np.uint32(TNs*4*3))
		of.write(np.uint32(TNs*4*3))
		if 'PartType0' in f.keys():
			of.write(np.float32(np.array(f['/PartType0/Coordinates'])/1000.0))
		if 'PartType1' in f.keys():
			of.write(np.float32(np.array(f['/PartType1/Coordinates'])/1000.0))
		if 'PartType2' in f.keys():
			of.write(np.float32(np.array(f['/PartType2/Coordinates'])/1000.0))
		if 'PartType3' in f.keys():
			of.write(np.float32(np.array(f['/PartType3/Coordinates'])/1000.0))
		if 'PartType4' in f.keys():
			of.write(np.float32(np.array(f['/PartType4/Coordinates'])/1000.0))
		if 'PartType5' in f.keys():
			of.write(np.float32(np.array(f['/PartType5/Coordinates'])/1000.0))
		of.write(np.uint32(TNs*4*3))

		# VEL
		write_head(of, "VEL ", np.uint32(TNs*4*3))
		of.write(np.uint32(TNs*4*3))
		if 'PartType0' in f.keys():
			of.write(np.float32(np.array(f['/PartType0/Velocities'])))
		if 'PartType1' in f.keys():
			of.write(np.float32(np.array(f['/PartType1/Velocities'])))
		if 'PartType2' in f.keys():		
			of.write(np.float32(np.array(f['/PartType2/Velocities'])))
		if 'PartType3' in f.keys():
			of.write(np.float32(np.array(f['/PartType3/Velocities'])))
		if 'PartType4' in f.keys():
			of.write(np.float32(np.array(f['/PartType4/Velocities'])))
		if 'PartType5' in f.keys():
			of.write(np.float32(np.array(f['/PartType5/Velocities'])))
		of.write(np.uint32(TNs*4*3))

		# ID
		write_head(of, "ID  ", np.uint32(TNs*4))
		of.write(np.uint32(TNs*4))
		if 'PartType0' in f.keys():
			of.write(np.uint32(np.array(f['/PartType0/ParticleIDs'])))
		if 'PartType1' in f.keys():		
			of.write(np.uint32(np.array(f['/PartType1/ParticleIDs'])))
		if 'PartType2' in f.keys():
			of.write(np.uint32(np.array(f['/PartType2/ParticleIDs'])))
		if 'PartType3' in f.keys():
			of.write(np.uint32(np.array(f['/PartType3/ParticleIDs'])))
		if 'PartType4' in f.keys():
			of.write(np.uint32(np.array(f['/PartType4/ParticleIDs'])))
		if 'PartType5' in f.keys():
			of.write(np.uint32(np.array(f['/PartType5/ParticleIDs'])))
		of.write(np.uint32(TNs*4))

		# MASS
		write_head(of, "MASS", np.uint32(TNs*4))
		of.write(np.uint32(TNs*4))
		if 'PartType0' in f.keys():
			of.write(np.float32(np.array(f['/PartType0/Masses'])))
		if 'PartType1' in f.keys():
			of.write(np.float32(np.array(f['/PartType1/Masses'])))
		if 'PartType2' in f.keys():
			of.write(np.float32(np.array(f['/PartType2/Masses'])))
		if 'PartType3' in f.keys():
			of.write(np.float32(np.array(f['/PartType3/Masses'])))
		if 'PartType4' in f.keys():
			of.write(np.float32(np.array(f['/PartType4/Masses'])))
		if 'PartType5' in f.keys():
			of.write(np.float32(np.array(f['/PartType5/Masses'])))
		of.write(np.uint32(TNs*4))

		# U <only gas>
		if 'PartType0' in f.keys():
			write_head(of, "U   ", np.uint32(TNp[0]*4))
			of.write(np.uint32(TNp[0]*4))
			of.write(np.float32(np.array(f['/PartType0/InternalEnergy'])))
			of.write(np.uint32(TNp[0]*4))

		# AGE <only Star> stellar formation time !!No BHs
		if 'PartType4' in f.keys():
			write_head(of, "AGE ", np.uint32(TNp[4]*4))
			of.write(np.uint32(TNp[4]*4))
			of.write(np.float32(np.array(f['/PartType4/StellarFormationTime'])))
			of.write(np.uint32(TNp[4]*4))

		# Z <gas + Star>
		if 'PartType0' in f.keys():
			if 'PartType4' in f.keys():
				Tbp = TNp[0]+TNp[4]
			else:
				Tbp = TNp[0]
			write_head(of, "Z   ", np.uint32(Tbp*4))
			of.write(np.uint32(Tbp*4))
			of.write(np.float32(np.array(f['/PartType0/Metallicity'])))
			if 'PartType4' in f.keys():
				of.write(np.float32(np.array(f['/PartType4/Metallicity'])))
			of.write(np.uint32(Tbp*4))

		if 'PartType0' in f.keys():
			# Rho <only gas>
			write_head(of, "RHO ", np.uint32(TNp[0]*4))
			of.write(np.uint32(TNp[0]*4))
			of.write(np.float32(np.array(f['/PartType0/Density'])))
			of.write(np.uint32(TNp[0]*4))

			# NE <only gas>
			write_head(of, "NE  ", np.uint32(TNp[0]*4))
			of.write(np.uint32(TNp[0]*4))
			of.write(np.float32(np.array(f['/PartType0/ElectronAbundance'])))
			of.write(np.uint32(TNp[0]*4))

			# NH <only gas>
			write_head(of, "NH  ", np.uint32(TNp[0]*4))
			of.write(np.uint32(TNp[0]*4))
			of.write(np.float32(np.array(f['/PartType0/NeutralHydrogenAbundance'])))
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
			of.write(hsml/1000.0)
			of.write(np.uint32(TNp[0]*4))

			# SFR <only gas>
			write_head(of, "SFR ", np.uint32(TNp[0]*4))
			of.write(np.uint32(TNp[0]*4))
			of.write(np.float32(np.array(f['/PartType0/StarFormationRate'])))
			of.write(np.uint32(TNp[0]*4))

