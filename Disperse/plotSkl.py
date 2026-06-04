from numpy import *
import sys
import os
import struct
import matplotlib.pyplot as plt
import matplotlib

if len(sys.argv) != 2:
	print('I require one argument -- the input file')
	print('I exit, because you providing is not fit!', str(sys.argv))
	exit(0)
else:	
	infile  = str(sys.argv[1])
	if infile[-8:] == '.a.NDskl':
		filetype = 'ascii'
	else:
		filetype = 'bin'

verbose              = True
if verbose:
	print(filetype)

# LOAD
# ---------- Binary File ---------- #
# Not working
if filetype == 'bin':
	NDFIELD_MAX_DIMS     = 20
	NDSKEL_DATA_STR_SIZE = 20

	with open(infile,"rb") as file:
		fC=file.read()

	ndims     = struct.unpack("i",fC[108:112])[0]
	nsegdata  = struct.unpack("i",fC[520:524])[0]
	nsegs     = struct.unpack("i",fC[512:516])[0]
	nnodes    = struct.unpack("i",fC[516:520])[0]
	nnodedata = struct.unpack("i",fC[524:528])[0]

	name      = struct.unpack(str(20*nsegdata)+"c",fC[536:536+20*nsegdata])
	name1=''
	for n in name:
		if n == b'\x00':
			name1+=' '
		else:
			name1+=n.decode("utf8")

	j       = 2*nsegs*ndims
	segpos  = array(struct.unpack(str(j)+"i",fC[552+20*nsegdata+20*nnodedata:552+20*nsegdata+20*nnodedata+2*nsegs*ndims*4])) 
	j       = nnodes*ndims 
	nodepos = array(struct.unpack(str(j)+"i",fC[560+20*(nsegdata+nnodedata)+8*nsegs*ndims:560+20*(nsegdata+nnodedata)+8*nsegs*ndims+nnodes*ndims*4]))
	nodes   = reshape(nodepos,(2,nnodes))
	print(nodepos)
	print(nodes)

	if verbose: #must be 0
		j = struct.unpack("i",fC[580+20*(nsegdata+nnodedata)+8*nsegs*(ndims+nsegdata)+4*nnodes*(ndims+2*nnodedata):584+20*(nsegdata+nnodedata)+8*nsegs*(ndims+nsegdata)+4*nnodes*(ndims+2*nnodedata)])[0]
		print(j)


# ---------- ASCII File ---------- #

else:
	delimiters = ['[CRITICAL POINTS]\n','[FILAMENTS]\n','[CRITICAL POINTS DATA]\n','[FILAMENTS DATA]\n']
	with open(infile, 'r') as fp:
		fC = list(enumerate(fp))
	fp.close()
	
	#jump header
	li = 0
	while fC[li][1] != delimiters[0]:
		# jump line
		li+=1
	li+=1 #to account for delimiters
	
	# CRITICAL POINTS
	#number of critical points
	num_cp = int(fC[li][1])
	
	#read positions of density maxima
	nodes  = []
	saddle = []
	critpt = [] #to store all critical points, will be used for filaments later
	while fC[li][1] != delimiters[1]:
		lC = fC[li][1].split()
		if len(lC) == 6:
			critpt.append([float(lC[1]),float(lC[2])])
			if int(lC[0]) == 2: #meaning we have a maximum
				nodes.append([float(lC[1]),float(lC[2])])
			elif int(lC[0]) == 1: #saddle
				saddle.append([float(lC[1]),float(lC[2])])
		del lC
		li+=1
	
	# FILAMENTS
	#number of critical points
	li+=1 #to account for delimiters
	#number of filaments
	num_f = int(fC[li][1])

	filaments = []	
	while fC[li][1] != delimiters[2]:
		lC = fC[li][1].split()
		if len(lC) == 3:
			filaments.append([int(lC[0]),int(lC[1])])
		li+=1	
	
nodes     = array(nodes)
saddle    = array(saddle)
filaments = array(filaments)
critpt    = array(critpt)

# PLOT
plt.plot(nodes[:,0],nodes[:,1],'r.')
plt.plot(saddle[:,0],saddle[:,1],'g.')
for i in range(num_f):
	plt.plot([critpt[filaments[i,0]][0],critpt[filaments[i,1]][0]],[critpt[filaments[i,0]][1],critpt[filaments[i,1]][1]],'k-')
plt.show()
