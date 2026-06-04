from numpy import*
import argparse

# input parameters
parser  = argparse.ArgumentParser(description='Prepare Halo.txt file for genP2cgrid.sh script. Keep only (main) halos with >98% high resolution particles and bin in log mass. Make sure AHF catalog header is line 0, all units must be AHF default ones (i.e. akpc/h or Msun/h). Args: < path+catalog name (str) >, < z (float) >, < h (float) >')
parser.add_argument('cat',nargs=1,type=str, help='Catalog absolute path (str)')
parser.add_argument('red',nargs=1,type=float, help='Redshift (float)')
parser.add_argument('hubble',nargs=1,type=float, help='Hubble param h (H0/100) (float)')
args    = parser.parse_args()
halocat = args.cat[0]
red     = args.red[0]
h       = args.hubble[0]

dM   = 0.5
ms   = arange(12.0,13.5,0.5)
ms   = ms[::-1]

#Load and filter
print("Catalog: "+halocat)
print("Loading and filtering")
f                = genfromtxt(halocat,names=True)
filter           = f['fMhires38']>=0.9
filter1          = f[filter]['hostHalo2']<=0

j = 1
for mhi in ms:
    print(mhi-dM)
    print(mhi)
    filterM          = (log10(f[filter][filter1]['Mhalo4']/h)>=(mhi-dM)) & (log10(f[filter][filter1]['Mhalo4']/h)<mhi)
    cathiresnosub    = f[filter][filter1][filterM]

    print(f[filter][filter1][filterM]['Mhalo4']/h)

    #write
    fout             = open(halocat+"_main_hires_Mlow_"+str(mhi-dM)+"_Mhigh_"+str(mhi)+".txt",'w')
    #Header
    fout.write("#gn; posX [cMpc]; posY [cMpc]; posZ [cMpc]; rVir [pkpc]; FOF-mass [Msol]; velX [km/s]; velY [km/s]; velZ [km/s]; spherMass [Msol]; S/N\n")
    posX = f[filter][filter1][filterM]['Xc6']/h/1000.
    posY = f[filter][filter1][filterM]['Yc7']/h/1000.
    posZ = f[filter][filter1][filterM]['Zc8']/h/1000.
    rVir = f[filter][filter1][filterM]['Rhalo12']/h/(1.+red)
    M    = f[filter][filter1][filterM]['Mhalo4']/h
    VX = f[filter][filter1][filterM]['VXc9']
    VY = f[filter][filter1][filterM]['VYc10']
    VZ = f[filter][filter1][filterM]['VZc11']
    for i in range(len(cathiresnosub)):
        new_line=str(j)+";"+str(posX[i])+";"+str(posY[i])+";"+str(posZ[i])+";"+str(rVir[i])+";"+str(M[i])+";"+str(VX[i])+";"+str(VY[i])+";"+str(VZ[i])+";"+str(M[i])+";"+str(2)+"\n"
        fout.write(new_line)
        j+=1
    fout.close()

    del filterM, cathiresnosub, posX, posY, posZ, rVir, M, VX, VY, VZ
