#Receives an array of masses and a bin of masses and return a count chart
import scipy as sp
from numpy import*
import math as m
import h5py
import sys 
import os

def countfct(mass,Mbin):
	"""
	Count the number of halos per mass bin and return it to a table
	"""
	#comptage #count
	count = zeros_like(Mbin)

	
	j=0
	for i in range(shape(mass)[0]):
		while log10(mass[i]) < Mbin[j] :
			j=j+1
			if j>=shape(Mbin)[0] :
				return count
		
		count[j] = count[j]+1

	return count
	
def loadhalos(ctype,path,filein):
	"""
	AHF halos in ASCII (subhalos must have been removed previously if necessary, Mass column must be 3)
	FOF/SUBFIND halos in hdf5
	"""
	if ctype == 'AHF':
		Mcat = loadtxt(path+"/"+filein,usecols=(3,))
	
	elif ctype == 'FOF':
		Mcat = []
		h5files = [f for f in os.listdir(path+"/") if os.path.isfile(os.path.join(path+"/", f)) and f[-4:] == 'hdf5']
		for i in h5files:
			f = h5py.File(path+"/"+i, "r")
			print(f)

			Mt = array(f['/Group/GroupMass'])
			lm = len(Mt)
			for m in xrange(lm): 
				Mcat.append(Mt[m]*1e10)
		Mcat.sort()
		Mcat=array(Mcat)		
		Mcat=Mcat[::-1]
	else:
		print("Catalogue type not recognised")
		sys.exit()

	return Mcat
