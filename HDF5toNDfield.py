import h5py
from numpy import *
import sys
import os
import struct

def write_ftype(fp,ftype,size):
	l     = len(ftype)
	dummy = b'\0'*(size-l)
	fp.write(struct.pack("i",int(size)))
	fp.write(ftype)
	fp.write(dummy)
	fp.write(struct.pack("i",int(size)))

def write_header(fp, head, size):
	dummy = b'\0'*160
	fill  = 18
	fp.write(struct.pack("i",int(size)))
	fp.write(head['comment'])
	fp.write(b'\0'*(80-len(comment)))
	fp.write(struct.pack("i",int(head['ndims'])))
	fp.write(struct.pack("i",int(head['ndims'])))
	fp.write(struct.pack("i",int(head['np'])))
	if ndims < NDFIELD_MAX_DIMS:
		for i in range(fill): 
			fp.write(struct.pack("i",int(0)))#to fill to 20
	fp.write(struct.pack("i",int(head['fdims'])))
	fp.write(struct.pack("i",int(head['datatype'])))
	fp.write(struct.pack("d",float(head['x0'])))
	fp.write(struct.pack("d",float(head['x0'])))
	if ndims == 3:
		fp.write(struct.pack("d",float(head['x0'])))
	if ndims < NDFIELD_MAX_DIMS:
		for i in range(fill):
			fp.write(struct.pack("d",float(0)))#to fill to 20
	fp.write(struct.pack("d",float(head['BoxSize'])))
	fp.write(struct.pack("d",float(head['BoxSize'])))
	if ndims == 3:
	        fp.write(struct.pack("d",float(head['BoxSize'])))
	if ndims < NDFIELD_MAX_DIMS:
		for i in range(fill):
	        	fp.write(struct.pack("d",float(0)))#to fill to 20
	fp.write(dummy) #to fill to size
	fp.write(struct.pack("i",int(size)))


if len(sys.argv) != 3:
	print('I require thwo arguments -- the input folder, and the dimensions (2D or 3D)')
	print('I exit, because you providing is not fit!', str(sys.argv))
	exit(0)
else:	
	folder  = str(sys.argv[1])+"/"
	ndims   = int(sys.argv[2])

verbose          = False
NDFIELD_MAX_DIMS = 20
l                = 2.0
BoxSize          = 2.0*l #135.54 
xc               = 66.2
yc               = 68.4
zc               = 68.1 
xmin             = (xc-l)
xmax             = (xc+l)
ymin             = (yc-l)
ymax             = (yc+l)
zmin             = (zc-l)
zmax             = (zc+l)

if verbose:
        print(folder)


# LOAD
h5files = [f for f in os.listdir(folder+"/") if os.path.isfile(os.path.join(folder+"/", f)) and f[0:4] == 'snap' and f[-4:] == 'hdf5']
h5files  = sorted(h5files)
h5files  = sorted(h5files,key=len)
if verbose:
        print(h5files)

print("Loading positions")
pos = []
np = 0
for i in h5files:
	print(str(i))	
	f     = h5py.File(folder+"/"+i, "r")
	
	numparts = f['Header'].attrs['NumPart_ThisFile'][1]
	postemp = array(f['/PartType1/Coordinates'],dtype=float32)*0.001
	
	#fill positions
	for j in range(numparts):
		if xmin <= postemp[j,0] <= xmax and ymin <= postemp[j,1] <= ymax and zmin <= postemp[j,2] <= zmax:
			pos.append(array([postemp[j,0]-xmin,postemp[j,1]-ymin,postemp[j,2]-zmin]))
			np+=1

pos=array(pos)


# WRITE
of = open(folder+'/'+h5files[0][:-7]+'.ND', 'wb')

# FILETYPE
print("Writing file type")
ftype=b'NDFIELD'
write_ftype(of,ftype,16)

# HEADER
print("Writing header")
f        = h5py.File(folder+h5files[0],"r")
attrs    = f['/Header'].attrs
h        = attrs['HubbleParam']
datatype = 256
x0       = 0.0
comment  = b'Arepo DM particles'
fdims    = 1

head = {'comment':comment, 'ndims':ndims, 'np':np, 'fdims':fdims, 'datatype':datatype, 'x0':x0, 'BoxSize':BoxSize}

size = 4*(NDFIELD_MAX_DIMS+3) + 8*(2*NDFIELD_MAX_DIMS) +(160+80)*1 #652
write_header(of,head, size)

# POS
print("Writing positions")
if ndims == 2:
	size = np*2*4	
	of.write(struct.pack("i",int(size)))
	for i in range(np):
		of.write(struct.pack("f",float(pos[i,1])))
		of.write(struct.pack("f",float(pos[i,0])))
		#of.write(struct.pack("f",float(pos[i,1])))
	of.write(struct.pack("i",int(size)))
else:
	size = np*3*4
	of.write(struct.pack("i",int(size)))
	for i in range(np):
		of.write(struct.pack("f",float(pos[i,0])))
		of.write(struct.pack("f",float(pos[i,1])))
		of.write(struct.pack("f",float(pos[i,2])))
	of.write(struct.pack("i",int(size)))

of.close()

print("Done")
