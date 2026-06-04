from numpy import*
import matplotlib.pyplot as plt

T     = 1e4
GHI   = 1e-12 #photoionisation rate
lognh = linspace(-6,3,900)
nh    = 10.**lognh

#Rahmati 2013
sivHI = 2.49e-18 #cm^2
fg    = 0.17 #baryon fraction taken to be universal
nhSS  = 6.73e-3*(sivHI/2.49e-18)**(-2./3.)*(T/1e4)**0.17*(GHI/1e-12)**(2./3.)*(fg/0.17)**(1./3.) #cm^-3
qR13  = nh/nhSS
fR13  = 0.98*(1.+qR13**1.64)**(-2.28) + 0.02*(1.+qR13)**(-0.84)

#Hopkins 2021 (Rahmati)
y     = 400.*nh*T**(-0.173)*(GHI/1e-12)**(-0.66)
fH21  = 0.98*(1.+y**1.64)**(-2.28) + 0.02*(1.+y*(1.+(nh/10)**4.))**(-0.84)
nhSSg = 0.0123 #2 times Rahmati
nhSSz = nhSSg*(GHI/1e-12)**(0.66)*(T/1e4)**0.173
qSSg  = nh/nhSSz
fH21g = 0.98/(1.+qSSg**1.64)**2.28 + 0.02/(1. + qSSg*(1.+1.e-4*nh*nh*nh*nh))**0.84

#Hopkins 2018 
sn0H  = 6e-18 #cm^2
pc2cm = 3.0856775812799585e18
lafit = 4.4*pc2cm**(T/1e4)**(-0.173)*(GHI/1e-12)**(-2./3.)
qH18  = sn0H*nh*lafit#\tau_a^ion in H18
fH18  = exp(-qH18)
fH18g = 1./(1.+qSSg*(1.+qSSg/2.*(1.+qSSg/3.*(1.+qSSg/4.*(1.+qSSg/5.*(1.+qSSg/6.*qSSg))))))

#plot
plt.loglog()
plt.plot(nh,fR13,'g--',label='Rahmati+ 13')
plt.plot(nh,fH18,'b--',label='Hopkins+ 18')        
plt.plot(nh,fH21,'r--',label='Hopkins+ 21') 
plt.plot(nh,fH18g,'b-',label='Hopkins +18 in Gizmo')
plt.plot(nh,fH21g,'r-',label='Hopkins+ 21 in Gizmo')
plt.legend(loc=0)
plt.xlabel('n_h [cm^-3]')
plt.ylabel('f')
plt.show()
