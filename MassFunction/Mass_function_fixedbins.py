#ECRITURE D'UN FICHIER AVEC LE NOMBRE MOYEN DE HALOS PAR BIN DE MASSE ET PLOT
#WRITING A FILE WITH THE AVERAGE NUMBER OF HALOS BY BASS OF MASS AND PLOT

import numpy as np
from numpy import*
import matplotlib.pyplot as plt
import math as m
import matplotlib
import argparse

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# What to do
parser = argparse.ArgumentParser(description='This script computes the mass function of halos from data dn(M)/dM. Do you want to check the plot? Default value is False. Args: NumofFile (int), Path to halo catalog (str), halos catalog name (str)')
parser.add_argument('--plot', action='store_true', help='Do check plot?')
parser.add_argument('numfiles',nargs=1,type=int, help='Number of catalogs')
parser.add_argument('path',nargs=1,type=str, help='Path to halo catalogs')
parser.add_argument('filename',nargs=1,type=str, help='Halo catalogs name')

args = parser.parse_args()
doplot=args.plot
path=args.path[0]
filein=args.filename[0]
NumOfFile=args.numfiles[0]

print("Plotting: "+str(doplot))

#fonction qui recoit le tableau des masses et retourne un tableau de count
#function that receives the array of masses and returns a count chart
import fonction_n_M
reload(fonction_n_M)
from fonction_n_M import*

#Cosmology [Mpc/h]
#V=1000.
V=250.   #boxsize Mpc
Om=0.3111
rho_crit=27.75e10
#V=1e-6
#rho_crit=277.5

Mmax_tab = empty(NumOfFile)
Mmin_tab = empty(NumOfFile)
M=[]

print("Loading...")
for i in xrange(NumOfFile):
	#Mtemp = loadtxt(path+"/"+str(i)+"/"+filein,usecols=(0,))
	Mtemp = loadtxt(path+"/"+filein,usecols=(3,)) 
        #Mtemp = loadtxt(path+"/"+filein)
        M.append(Mtemp)
	Mmax_tab[i]=max(Mtemp)
	Mmin_tab[i]=min(Mtemp)

#Mmax=log10(max(Mmax_tab))
#Mmin=log10(min(Mmin_tab))
M=asarray(M)

Mmax=15.4500000
Mmin=10.9500000

print Mmax
print Mmin

#Tableau des bins de masse
#Table of mass bins
dlogM=0.5
Mbin = arange(Mmin,Mmax,dlogM)
Mbin = Mbin[::-1]
SMb = shape(Mbin)[0]


#pour ploter le milieu du bin
#to plot the middle of the bin 
Mbin_milieu = Mbin+dlogM/2. 

#count: #halos/bin masse pour chaque fichier
#count: # halos / mass bin for each file
count = empty((NumOfFile,SMb))

#les dndM pour chaques realisations
#the dndM for each achievement
dndM = empty((NumOfFile,SMb))

print("Counting...")

#moyenne des data points file
#average of data points file
for i in xrange(NumOfFile) :
	count[i]=countfct(M[i],Mbin)
	for j in xrange(SMb) :
		dndM[i,j] = 10.0**(Mbin_milieu[j])/(Om*rho_crit)*count[i,j]/(V**3.0*dlogM)*1.0/log(10.0)	

#calculer la moyenne de chaque colonne (#halos moy/bin de masse)
#calculate the average of each column (#halos moy / bin mass)
countmoy = empty(SMb, dtype = 'd')
for i in xrange(SMb) :
	countmoy[i] = m.fsum(count[:,i])/NumOfFile

dndM_moy = empty(SMb)
dndM_moy = 10.0**(Mbin_milieu)/(Om*rho_crit)*countmoy/(V**3.0*dlogM)*1.0/log(10.0)

########################################################################
#error bars

#variance sur chaque population de count (1 popu = 1 colonne)
#variance on each population of count (1 popu = 1 column)
##s_sq_dndM = empty(SMb)

#erreur sur chaque element de countmoy (error bars a ploter)
#error on each element of countmoy (error bars a ploter)
##s_dndMmoy = empty(SMb)

##for i in range(SMb) :
##	s_sq_dndM[i] = m.fsum((dndM_moy[i]-dndM[:,i])**2.0)/(NumOfFile-1.0)
##	s_dndMmoy[i] = m.sqrt(s_sq_dndM[i]/NumOfFile)	#1sigma error 
        
########################################################################
#Tinker theo (if toplot==True)

##if doplot==True:
	#loader
##	tinker = loadtxt(path+"/Tinker_MF_z=0.dat").T

########################################################################
#ploter

##if doplot==True:
##	plt.plot(Mbin_milieu, log10(dndM_moy), 'r.',label=r'Data points')
##	plt.errorbar(Mbin_milieu,log10(dndM_moy), xerr=0.0, yerr=log10((dndM_moy+s_dndMmoy)/dndM_moy),fmt='r.')
##	plt.plot(log10(tinker[0]),log10(tinker[1]), 'k-',label=r'Tinker et al. (2008)')
##	plt.xlabel(r'$M [\log M_\odot h^{-1}]$')
##	plt.ylabel('$\log [M^2 (\rho_m*dn/dM)^{-1}]$')
##	plt.axis([12,16,-5.0,-0.8])
##	plt.legend(loc=0)
##	plt.show()

########################################################################
#ecrire le bin de masse, la moyenne du nombre de halos, l'erreur
#write the mass bin, the average number of halos, the error
#f=open(path+"/logM_dndlogM_sdndlogM.dat", "w")
#f=open(path+"/fid1_logM_dndlogM_.dat", "w")

#f=open("/u/titouan/MassFunction_titouan/sbf_logM_ax_0.4_z=0_withbaryons.dat", "w")  

#onlyfor N_M
f=open("/u/titouan/MassFunction_titouan/2f_fid4_logM_dndM_0.5_z=3_fixed.dat","w")

for i in range(SMb):
	#new_line = str(Mbin_milieu[i]) +"\t"+str((Om*rho_crit)*dndM_moy[i]/10.0**Mbin_milieu[i])+"\t"+str((Om*rho_crit)*s_dndMmoy[i]/10.0**Mbin_milieu[i])+"\n"
        #new_line = str(Mbin_milieu[i]) +"\t"+str(dndM_moy[i])+"\n"
#for dndM without anything else(ax)
        new_line = str(Mbin_milieu[i]) +"\t"+str((Om*rho_crit)*dndM_moy[i]/(10.0**Mbin_milieu[i]*10.0**Mbin_milieu[i]))+"\n"
#only for N_M
        #new_line = str(Mbin_milieu[i]) +"\t"+str(countmoy[i])+"\n"

	f.write(new_line)
f.close()
