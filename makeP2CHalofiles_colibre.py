from numpy import*
import argparse
import unyt
import swiftsimio as sw
from operator import itemgetter

# input parameters
parser  = argparse.ArgumentParser(description='Prepare P2C compatible Halo.txt file from SOAP files. Keep only main halos and bin in log mass. Args: < halo catalog path (str) >, < snapshot number (int) >')
parser.add_argument('cat',nargs=1,type=str, help='Catalog absolute path (str)')
parser.add_argument('snap_nr',nargs=1,type=int, help='Snapshot number (int)')
args    = parser.parse_args()
halocat = args.cat[0]
snap_nr = args.snap_nr[0]

dM   = 0.5
ms   = arange(12.5,14.0,0.5)
ms   = ms[::-1]

#Load and filter
print("Loading and filtering")

soap_filename = f'{halocat}/halo_properties_{snap_nr:04}.hdf5'
soap          = sw.load(soap_filename)
red           = round(soap.metadata.redshift,1)
mask          = soap.spherical_overdensity_200_mean.total_mass.to('Msun') >= 1e12
M             = array(soap.spherical_overdensity_200_mean.total_mass.to('Msun')[mask]) #these are all main halos as they do not compute spherical overdensity properties for satellites
pos           = array(soap.spherical_overdensity_200_mean.centre_of_mass.to('Mpc')[mask]) #cMpc
vel           = array(soap.spherical_overdensity_200_mean.centre_of_mass_velocity.to('km/s')[mask]) #km/s
Rvir          = array(soap.spherical_overdensity_200_mean.soradius.to('kpc').to_physical()[mask]) #pkpc
P2Ctab        = vstack((pos[:,0],pos[:,1],pos[:,2],Rvir,M,vel[:,0],vel[:,1],vel[:,2],M))
P2Ctab        = array(sorted(P2Ctab.T,key=itemgetter(4),reverse=True)).T #sorted by decreasing order of mass

j = 1
for mhi in ms:
    print(mhi-dM)
    print(mhi)
    filterM = (log10(P2Ctab[4])>=(mhi-dM)) & (log10(P2Ctab[4])<mhi) 

    #write
    fout    = open("/cosma8/data/do019/dc-laze1/COLIBRE/L0100N1504/200_snap_z"+str(red)+"_halos_sorted_main_hires_Mlow_"+str(mhi-dM)+"_Mhigh_"+str(mhi)+".txt",'w')
    #Header
    fout.write("#gn; posX [cMpc]; posY [cMpc]; posZ [cMpc]; rVir [pkpc]; FOF-mass [Msol]; velX [km/s]; velY [km/s]; velZ [km/s]; spherMass [Msol]; S/N\n")
    #grids
    for i in range(len(P2Ctab[4][filterM])):
        new_line=str(j)+";"+str(P2Ctab[0][filterM][i])+";"+str(P2Ctab[1][filterM][i])+";"+str(P2Ctab[2][filterM][i])+";"+str(P2Ctab[3][filterM][i])+";"+str(P2Ctab[4][filterM][i])+";"+str(P2Ctab[5][filterM][i])+";"+str(P2Ctab[6][filterM][i])+";"+str(P2Ctab[7][filterM][i])+";"+str(P2Ctab[8][filterM][i])+";"+str(2)+"\n"
        fout.write(new_line)
        j+=1
    fout.close()

    del filterM
