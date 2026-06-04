import numpy as np

gamma     = 5./3.
mu        = 0.593
mH        = 1.67e-24  #g
kb        = 1.3086e-16 #cgs
Msol_in_g = 1.989e33 #g
G         = 6.674e-8 #cgs
kpc_in_cm = 3.086e21 #cm

#DALYA
Trhod = np.load('/data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/density_temperature_volume_in_50-160ckpc_Hcat.npz')
Td    = Trhod['temperature']
rhod  = Trhod['density']
vold  = Trhod['volume']

hotd  = np.log10(Td) >= 5.6
coldd = (np.log10(Td) < 5) & (np.log10(Td) > 3) #4.2 3.6

hotrhod  = np.average(rhod[hotd],weights=vold[hotd])
coldrhod = np.average(rhod[coldd],weights=vold[coldd])
hotTd    = np.average(Td[hotd],weights=vold[hotd])
coldTd   = np.average(Td[coldd],weights=vold[coldd])

meanMd = 1.8378467611037332e12 #Msol
meanRd = 84.19818667715968 #pkpc
Tvird  = G * mu*mH * meanMd*Msol_in_g / (gamma * kb *meanRd*kpc_in_cm)

print("DaLya")
print("<T_hot> [1e6 K] = "+str(hotTd/1e6))
print("<T_vir> [1e6 K] = "+str(Tvird/1e6))
print("<T_cold> [1e4 K] = "+str(coldTd/1e4))
print("<n_hot> [cm^-3] = "+str(hotrhod))
print("<n_cold> [cm^-3] = "+str(coldrhod))
print("<n_hot>/<n_cold> = "+str(hotrhod/coldrhod))
print("<T_cold>/<T_hot> = "+str(coldTd/hotTd))
print("ratio P_hot/P_cold = "+str((hotrhod/coldrhod)/(coldTd/hotTd)))


#COLIBRE
Trhoc = np.load('/data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/density_temperature_volume_in_50-160ckpc_Hcat.npz')
Tc    = Trhoc['temperature']
rhoc  = Trhoc['density']
volc  = Trhoc['volume']

hotc  = np.log10(Tc) >= 5.3
coldc = (np.log10(Tc) < 5) & (np.log10(Tc) > 3) #4.4 3.5

hotrhoc  = np.average(rhoc[hotc],weights=volc[hotc])
coldrhoc = np.average(rhoc[coldc],weights=volc[coldc])
hotTc    = np.average(Tc[hotc],weights=volc[hotc])
coldTc   = np.average(Tc[coldc],weights=volc[coldc])

meanMc = 1.388125e12 #Msol
meanRc = 76.9042985 #pkpc
Tvirc  = G * mu*mH * meanMc*Msol_in_g / (gamma * kb *meanRc*kpc_in_cm)

print("COLIBRE")
print("<T_hot> [1e6 K] = "+str(hotTc/1e6))
print("<T_vir> [1e6 K] = "+str(Tvirc/1e6))
print("<T_cold> [1e4 K] = "+str(coldTc/1e4))
print("<n_hot> [cm^-3] = "+str(hotrhoc))
print("<n_cold> [cm^-3] = "+str(coldrhoc))
print("<n_hot>/<n_cold> = "+str(hotrhoc/coldrhoc))
print("<T_cold>/<T_hot> = "+str(coldTc/hotTc))
print("ratio P_hot/P_cold = "+str((hotrhoc/coldrhoc)/(coldTc/hotTc)))

#HELLO
Trhoh = np.load('/data/titouan/HELLO/CDM/g9.20e12/density_temperature_volume_in_50-160ckpc_Hcat.npz')
Th    = Trhoh['temperature']
rhoh  = Trhoh['density']
volh  = Trhoh['volume']

hoth  = np.log10(Th) >= 5.5
coldh = (np.log10(Th) < 5) & (np.log10(Th) > 3) #4.5 3.5

hotrhoh  = np.average(rhoh[hoth],weights=volh[hoth])
coldrhoh = np.average(rhoh[coldh],weights=volh[coldh])
hotTh    = np.average(Th[hoth],weights=volh[hoth])
coldTh   = np.average(Th[coldh],weights=volh[coldh])

meanMh = 1.101989444573e12 #Msol
meanRh = 69.38280732607619 #pkpc
Tvirh  = G * mu*mH * meanMh*Msol_in_g / (gamma * kb *meanRh*kpc_in_cm)

print("HELLO")
print("<T_hot> [1e6 K] = "+str(hotTh/1e6))
print("<T_vir> [1e6 K] = "+str(Tvirh/1e6))
print("<T_cold> [1e4 K] = "+str(coldTh/1e4))
print("<n_hot> [cm^-3] = "+str(hotrhoh))
print("<n_cold> [cm^-3] = "+str(coldrhoh))
print("<n_hot>/<n_cold> = "+str(hotrhoh/coldrhoh))
print("<T_cold>/<T_hot> = "+str(coldTh/hotTh))
print("ratio P_hot/P_cold = "+str((hotrhoh/coldrhoh)/(coldTh/hotTh)))
