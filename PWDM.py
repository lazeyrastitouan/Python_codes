#Compute TF of a WDM model given a WDM mass and a CDM TF

from numpy import*

mWDM  = 1.5 #in keV
TCDM  = loadtxt("/data/titouan/Eagle_cosmo_gizmo/Eagle_cosmo_transfer_out_z0.dat").T
PSCDM = loadtxt("/data/titouan/Eagle_cosmo_gizmo/Eagle_cosmo_matterpower_z0.dat").T
TWDM  = TCDM
PSWDM = PSCDM
kT    = TCDM[0]
kP    = PSCDM[0]

#parameters
h  = 0.677
Om = 0.307
nu = 1.12 #1.2
a  = 0.049 * (1./mWDM)**1.11 * (Om/0.25)**0.11 *(h/0.7)**1.22 #Mpc/h Viel 2005
TofkT = (1.+(a*kT)**(2.*nu))**(-5./nu) #Viel 2005 (Ma 2000: (1+(a*k)**nu)**(-5/nu)
TofkP = (1.+(a*kP)**(2.*nu))**(-5./nu)

PSWDM[1] = PSWDM[1]*TofkP**2.

for i in [1,2,3,4,6,7,8]:
	TWDM[i]=TofkT*TWDM[i]

fT = open("/data/titouan/Eagle_cosmo_gizmo/Eagle_cosmo_transfer_WDM_1p5kev_transfer_out_z0.dat",'w')
fP = open("/data/titouan/Eagle_cosmo_gizmo/Eagle_cosmo_WDM_1p5kev_matterpower_z0.dat",'w')
for i in range(len(kP)):
	fP.write(str(kP[i])+"\t"+str(PSWDM[1,i])+"\n")
for i in range(len(kT)):
	fT.write(str(TWDM[0,i])+"\t"+str(TWDM[1,i])+"\t"+str(TWDM[2,i])+"\t"+str(TWDM[3,i])+"\t"+str(TWDM[4,i])+"\t"+str(TWDM[5,i])+"\t"+str(TWDM[6,i])+"\t"+str(TWDM[7,i])+"\t"+str(TWDM[8,i])+"\t"+str(TWDM[9,i])+"\t"+str(TWDM[10,i])+"\t"+str(TWDM[11,i])+"\t"+str(TWDM[12,i])+"\n")
fT.close()
fP.close() 
