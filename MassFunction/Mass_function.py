#WRITING A FILE WITH THE AVERAGE NUMBER OF HALOS IN BINS OF MASS AND PLOT
#COSMOLOGY SPECIDFIE IN cosmo.py

from numpy import*
import matplotlib.pyplot as plt
import matplotlib
import argparse
import sys
import os
import h5py

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# What to do
parser = argparse.ArgumentParser(description='This script computes the mass function of halos. Cosmology specified in cosmo.py. Do you want to check the plot? Default value is False. Args: Catalog type (AHF/FOF), NumofFile (int), Path to halo catalog (str), halos catalog name (str)')
parser.add_argument('--plot', action='store_true', help='Do check plot?')
parser.add_argument('cattype',nargs=1,type=str, help='Catalog type (AHF/FOF)')
parser.add_argument('numfiles',nargs=1,type=int, help='Number of catalogs (realisations)')
parser.add_argument('path',nargs=1,type=str, help='Path to halo catalogs')
parser.add_argument('filename',nargs=1,type=str, help='Halo catalogs name')

args 	  = parser.parse_args()
doplot	  = args.plot
ctype	  = args.cattype[0]
path	  = args.path[0]
filein	  = args.filename[0]
NumOfFile = args.numfiles[0]

print("Plotting: "+str(doplot))

#function that receives the array of masses and returns a count chart
from fonction_n_M import*

#Cosmology [Mpc/h]
import cosmo
V  	 = cosmo.V
Om 	 = cosmo.Om
h 	 = cosmo.h
rho_crit = cosmo.rho_crit

Mmax_tab = empty(NumOfFile)
Mmin_tab = empty(NumOfFile)
M=[]

print("Loading...")
for i in range(NumOfFile):
	Mtemp = loadhalos(ctype,path,filein)
	M.append(Mtemp)
	Mmax_tab[i]=max(Mtemp)
	Mmin_tab[i]=min(Mtemp)

Mmax=log10(max(Mmax_tab))
Mmin=log10(min(Mmin_tab))
M=array(M)

print(Mmax)
print(Mmin)

#Table of mass bins
dlogM = 0.33
Mbin = arange(Mmin,Mmax,dlogM)
Mbin = Mbin[::-1]
SMb = shape(Mbin)[0]

#to plot the middle of the bin 
Mbin_milieu = Mbin+dlogM/2. 

#count: # halos / mass bin for each file
count = empty((NumOfFile,SMb))

#the dndM for each reals
dndM = empty((NumOfFile,SMb))

print("Counting...")

#average of data points file
for i in range(NumOfFile) :
	count[i]=countfct(M[i],Mbin)
	#only needed for errorbars
	if NumOfFile > 1:
		for j in range(SMb) :
			dndM[i,j] = 10.0**(Mbin_milieu[j])/(Om*rho_crit)*count[i,j]/(V**3.0*dlogM)*1.0/log(10.0)	

#calculate the average of each column (#halos moy / bin mass)
countmoy = empty(SMb, dtype = 'd')
for i in range(SMb) :
	countmoy[i] = sum(count[:,i])/NumOfFile

#for comparison to Tinker
dndM_moy = empty(SMb)
dndM_moy = 10.0**(Mbin_milieu)/(Om*rho_crit)*countmoy/(V**3.0*dlogM)*1.0/log(10.0)

########################################################################
#error bars

#variance on each population of count (1 popu = 1 column)
s_sq_dndM = empty(SMb)

#error on each element of countmoy (error bars a ploter)
s_dndMmoy = empty(SMb)
s_countmoy = empty(SMb)

if NumOfFile > 1:
	for i in range(SMb) :
		s_sq_dndM[i] = sum((dndM_moy[i]-dndM[:,i])**2.0)/(NumOfFile-1.0)
		s_dndMmoy[i] = sqrt(s_sq_dndM[i]/NumOfFile)	#1sigma error 
        
		s_countmoy[i] = sqrt(sum((countmoy[i]-count[:,i])**2.0)/((NumOfFile-1.0)*NumOfFile))

########################################################################
#Tinker theo (if doplot==True)

#if doplot==True:
	#loader
	#tinker = loadtxt(path+"/Tinker_MF_z=0.dat").T

########################################################################
#plot

if doplot==True:
    fout="logM_MsqorhomdndM_"+filein+".pdf"
    plt.plot(Mbin_milieu, log10(dndM_moy), 'r.',label=r'Data points')
    plt.errorbar(Mbin_milieu,log10(dndM_moy), xerr=0.0, yerr=log10((dndM_moy+s_dndMmoy)/dndM_moy),fmt='r.')
    #plt.plot(log10(tinker[0]),log10(tinker[1]), 'k-',label=r'Tinker et al. (2008)')
    plt.xlabel(r'$M [\log M_\odot h^{-1}]$')
    plt.ylabel('$\log [M^2 (\rho_m*dn/dM)^{-1}]$')
    #plt.axis([12,16,-5.0,-0.8])
    plt.legend(loc=0)
    savefig(fout, bbox_inches='tight')
    plt.show()

########################################################################
#write the mass bin, the average number of halos, the error

# N_M
#f=open(path+"/logM_NofM_"+filein+".dat","w")
f=open(path+"/logM_MsqorhomdndM_sMsqorhomdndM_"+filein+".dat","w")

for i in range(SMb):
	#new_line = str(log10(10.0**Mbin_milieu[i]/h)) +"\t"+str(countmoy[i])+"\n"#+str(s_countmoy[i])+"\n"
	new_line = str(log10(10.0**Mbin_milieu[i])) +"\t"+str(dndM_moy[i])+"\t"+str(s_dndMmoy[i])+"\n"

	f.write(new_line)
f.close()
