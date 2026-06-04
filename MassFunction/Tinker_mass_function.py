#FIT DE LA REALISATION 1 DE AHF AVEC LA FONCTION DE MASSE DE TINKER ET AL.
#FIT OF REALIZATION 1 OF AHF WITH THE MASS FUNCTION OF TINKER ET AL.

import numpy as np
from numpy import*
import matplotlib.pyplot as plt
from matplotlib.pyplot import*
import scipy as sp
#import math as m
#from math import*

########################################################################
#definition des parametres (Tinker et al (2008)) + cosmo fiducielle
#definition of the parameters (Tinker et al (2008)) + fiducial cosmo
A = 0.186
a = 1.47
b = 2.57
c = 1.19
R = 8.0 #Mpc/h
#R = 8e10 #kpc/h 
#h = 0.7
h=0.6766
rho_crit = 27.75e10 #h^2 Msun/Mpc^3
#rho_crit = 277.5 #h^2 Msun/kpc^3
Omega_m = 0.3111 #0.27
rho_m = rho_crit*Omega_m #h^2 Msun/Mpc^3

#charger le PS et creer la masse
#load the PS and Create the mass
PStxt = loadtxt("/cobra/ptmp/titouan/baryon-cdm/CAMB-matterpows/test_matterpower_z=0_CDM_fid.dat")
#("Pk_lcdm_transfer_out.dat_s8_0.7913_ns_0.9500")
k = PStxt[:,0] #h/Mpc
PS = PStxt[:,1] #(Mpc/h)^3

Mlog = zeros(1200,dtype = 'd')
M = zeros(1200,dtype = 'd')
for i in range(shape(M)[0]) :
	Mlog[i] = 10.0 + i/200.0 
	
M = 10**Mlog #Msun/h	
R_M = (3.0*M/(4.0*pi*rho_m))**(1.0/3.0) #Mpc/h
#R_M va de 0.68 a 14.67
#R_M goes from 0.68 to 14.67

Sk = shape(k)[0]
SM = shape(M)[0]

########################################################################
#test : calcul de sigma_8
#calculation of sigma_8

#definition du filtre Tophat
#definition of the Tophat Filter 
W = zeros(Sk,dtype = 'd')
W = 3.0/(k*R)**3.0*(sin(k*R)-k*R*cos(k*R))

#integration => sigma_8_sq
d_sigma = 0.0
sigma_8_sq = 0.0
u=0

while u < (Sk-1) :
	d_sigma = (k[u+1]-k[u])*1.0/(2.0*(pi)**2.0)*(W[u]**2.0*PS[u]*k[u]**2.0+W[u+1]**2.0*PS[u+1]*k[u+1]**2.0)/2.0
	sigma_8_sq = sigma_8_sq + d_sigma
	u=u+1

print "sigma_8_sq = {}" .format(sigma_8_sq) #0.626159 au lieu 0.626155

sigma = sqrt(sigma_8_sq)
f_sigma = A*((sigma/b)**a+1.0)*exp(-c/sigma_8_sq)

########################################################################
#calcul de la mass function pour tout R
#mass function calculation for all R

#filtre tophat pour tout R
#tophat filter for all R
Cc = 4.0*pi*rho_m/3.0
W_M = zeros((SM, Sk),dtype = 'd')
for i in range(SM) :
	for l in range(Sk) :
		W_M[i,l] = 3.0/(k[l]*R_M[i])**3.0*(sin(k[l]*R_M[i])-k[l]*R_M[i]*cos(k[l]*R_M[i]))

#~ plot(k, W_M[200,:],'k-',k,W_M[350,:],'m-',k, W_M[450,:],'b-',k,W_M[500],'r-',k,W_M[650,:],'g-',k, W_M[900,:], 'y-')
#~ xlabel('k [h/Mpc]')
#~ ylabel('W(kR)')
#~ axis([0,15,-0.2,1.2])
#~ show()

#calcul de sigma pour tous les R_M
#sigma calculation for All R_M
s_m_sq = zeros(SM, dtype = 'd')
d_s_m = zeros(SM, dtype = 'd')
s_m = zeros(SM, dtype = 'd')
f_s_m = zeros(SM, dtype = 'd')

for i in range(SM) : 
	s = 0
	while s < Sk-1 :
		d_s_m[i] = (k[s+1]-k[s])*1.0/(2.0*(pi)**2.0)*(W_M[i,s]**2.0*PS[s]*k[s]**2.0+W_M[i,s+1]**2.0*PS[s+1]*k[s+1]**2.0)/2.0
		s_m_sq[i] = s_m_sq[i] + d_s_m[i]
		s=s+1

s_m = sqrt(s_m_sq)
f_s_m = A*((s_m/b)**(-a)+1.0)*exp(-c/s_m_sq) 
	
#~ plot(log10(M), log10(s_m_sq), 'r-')
#~ show()
#~ 
#~ plot(log10(1.0/s_m),log10(f_s_m), 'r-')
#~ xlabel('log(1/sigma)')
#~ ylabel('log(f(sigma))') 
#~ axis([-0.6,0.5,-3,0])
#~ show()

#dsigma/dM
#refaire une integrale pour d_sigma/d_M
#redo an integral for d_sigma / d_M
d_dSdM = zeros(SM, dtype = 'd')
dSdM = zeros(SM, dtype = 'd')
 
for i in range(shape(M)[0]) :
	j=0
	while j < (shape(k)[0]-1) :
		d_dSdM[i] = (k[j+1]-k[j])/(2.0*(pi)**2.0*M[i]**2.0)*PS[j]*2.0*W_M[i,j]/k[j]*Cc*(3.0*(M[i]*k[j]**3.0/Cc)**(1.0/3.0)*cos((M[i]*k[j]**3.0/Cc)**(1.0/3.0))+((M[i]*k[j]**3.0/Cc)**(2.0/3.0)-3.0)*sin((M[i]*k[j]**3.0/Cc)**(1.0/3.0)))
		dSdM[i] = dSdM[i] + d_dSdM[i]
		j=j+1

#mass function
dndM = zeros(SM, dtype = 'd')
dndM = - f_s_m/(2.0*s_m_sq)*rho_m/M*dSdM

########################################################################
#rajouter les data points (ici que d'un fichier)

#aller ouvrir le fichier avec les masses des halos 
#Z=0
#halos=loadtxt("200_snap_LCDM_z0.000_sorted")

#definir un tableau de bin pour compter et ploter (bin de 0.5 en 0.5)
#Mbin = array([15.3,15.2,15.1,15.0,14.9,14.8,14.7,14.6,14.5,14.4,14.3,14.2,14.1,14.0,13.9,13.8,13.7,13.6,13.5,13.4,13.3,13.2,13.1,13.0,12.9,12.8,12.7,12.6,12.5,12.4,12.3,12.2,12.1])
#SMb = shape(Mbin)[0]
#count = zeros(SMb)

#faire des test sur les elements de la collone masse et ajouter a chaque bin de masse
#M_dat=halos[:,3]#compte des colonnes start a 0, ligne des en-tete entre pas dans le compte des lignes
#Sm = shape(M_dat)[0]
#j=0

#for i in xrange(Sm):	
#	while log10(M_dat[i]) < Mbin[j] :
#		j=j+1
#	count[j] = count[j]+1

#print "count {}".format(count)

#pour ploter le milieu du bin
#Mbin_milieu = zeros(SMb, dtype='d')
#Mbin_milieu = Mbin+0.05

#n = zeros(SMb, dtype='d') #correspond au M^2/rho_m dn/dM de Tinker
#n = 10.0**(Mbin_milieu)/rho_m*count/(500.0**3.0*0.1)*1.0/log(10.0)
 
########################################################################
#comparaison avec Christian

#~ #loader
#~ chris = loadtxt("M_sigmaM_dlogsigma_dlogM.dat")
#~ M_chris = chris[:,0]
#~ sigma_chris = chris[:,1]
#~ dlnsdlnM_chris = chris[:,2]
#~ 
#~ #mass function
#~ SMc = shape(M_chris)[0] 
#~ 
#~ f_sigma_chris = zeros(SMc)
#~ f_sigma_chris = A*((sigma_chris/b)**(-a)+1.0)*exp(-c/sigma_chris**2.0)
#~ 
#~ toplot_chris = zeros(SMc)
#~ toplot_chris = f_sigma_chris*dlnsdlnM_chris

#perfect agreement

########################################################################
#ploter
toplot = zeros(SM, dtype = 'd')
toplot = M**2.0/rho_m*dndM
###
##plot(log10(M), log10(toplot), 'b-')#, Mbin_milieu, log10(n), 'r.')#, log10(M_chris),log10(toplot_chris), 'g-')
##xlabel('log(M [Msun/h])')
##ylabel('log(M^2/rho_m*dndM)')
##axis([10,16,-3.6,-0.8])
##show()

########################################################################
#ecriture dans un fichier
f=open("/u/titouan/MassFunction_titouan/Tinker-dndM_CDM_z=0_fid.dat", "w") 

for i in range(SM):
	#if i < SMb :
	new_line=str(M[1199-i])+"\t"+str(dndM[1199-i])+"\n"# +"	"+str(10.0**(Mbin_milieu[i])) +"	"+str(n[i])+"\n"
	#else :
		#new_line=str(M[1199-i])  +"	"+str(toplot[1199-i]) +"\n"
	 
	f.write(new_line)
 
f.close()
