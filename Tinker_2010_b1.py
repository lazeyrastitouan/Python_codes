# -*- coding: utf-8 -*-
import numpy as np
from numpy import*
import matplotlib.pyplot as plt
from matplotlib.pyplot import*
import scipy as s
import math

#definition of the parameters (Tinker et al (2010))

D = 200.
y = log10(D)
A = 1.0 + 0.24*y*math.exp(-(4.0/y)**4)
a = 0.44*y - 0.88
B = 0.183
b = 1.5
C = 0.019 + 0.107*y + 0.19*math.exp(-(4/y)**4)
c = 2.4

R = 8.0 #Mpc/h  

h=0.6766
Omega_m = 0.3111
rho_crit = 27.75e10 #h^2 Msun/Mpc^3 
rho_m = rho_crit*Omega_m #h^2 Msun/Mpc^3

#load the PS and Create the mass
PStxt = loadtxt("/cobra/ptmp/titouan/baryon-cdm/CAMB-matterpows/test_matterpower_z=3_CDM+baryon_fid.dat")

k = PStxt[:,0] #h/Mpc 
PS = PStxt[:,1] #(Mpc/h)^3

Mlog = zeros(1200,dtype = 'd')
M = zeros(1200,dtype = 'd')

for i in range(shape(M)[0]) :
    Mlog[i] = 10.0 + i/200.0

M = 10**Mlog #Msun/h  
R_M = (3.0*M/(4.0*pi*rho_m))**(1.0/3.0) #Mpc/h

#R_M goes from 0.68 to 14.67 


Sk = shape(k)[0]
SM = shape(M)[0]


########################################################################   
#calculation of sigma_8 

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
#f_sigma = A*((sigma/b)**a+1.0)*exp(-c/sigma_8_sq)
delta_c = 1.686
nu8 = delta_c / sigma
b_nu8 = 1.0 - A*(nu8**a/(nu8**a+delta_c**a))+B*nu8**b+C*nu8**c

######################################################################## 

#mass function calculation for all R  

Cc = 4.0*pi*rho_m/3.0
W_M = zeros((SM, Sk),dtype = 'd')
for i in range(SM) :
    for l in range(Sk) :
        W_M[i,l] = 3.0/(k[l]*R_M[i])**3.0*(sin(k[l]*R_M[i])-k[l]*R_M[i]*cos(k[l]*R_M[i]))

#sigma calculation for All R_M 
s_m_sq = zeros(SM, dtype = 'd')
d_s_m = zeros(SM, dtype = 'd')
s_m = zeros(SM, dtype = 'd')
nu = zeros(SM, dtype = 'd')
b_nu = zeros(SM, dtype = 'd')

for i in range(SM) :
    s = 0
    while s < Sk-1 :
        d_s_m[i] = (k[s+1]-k[s])*1.0/(2.0*(pi)**2.0)*(W_M[i,s]**2.0*PS[s]*k[s]**2.0+W_M[i,s+1]**2.0*PS[s+1]*k[s+1]**2.0)/2.0
        s_m_sq[i] = s_m_sq[i] + d_s_m[i]
        s=s+1

s_m = sqrt(s_m_sq)
#f_s_m = A*((s_m/b)**(-a)+1.0)*exp(-c/s_m_sq)
nu = delta_c / s_m
b_nu = 1.0 - A*((nu**a)/(nu**a+delta_c**a))+B*nu**b+C*nu**c



toplot = zeros(SM, dtype = 'd')
toplot = b_nu

f=open("/u/titouan/MassFunction_titouan/Tinker_b1_z=3.dat", "w")

for i in range(SM):
    new_line=str(M[1199-i])+"\t"+str(toplot[1199-i])+"\t"+str(nu[1199-i])+"\n"
    f.write(new_line)
f.close()
