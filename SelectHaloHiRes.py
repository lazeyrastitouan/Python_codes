from numpy import*
import argparse
from decimal import Decimal

# input parameters
parser  = argparse.ArgumentParser(description='Keep only (main) halos in high resolution region. Make sure AHF catalog header is line 0. Args: < path+catalog name (str) >')
parser.add_argument('cat',nargs=1,type=str, help='Catalog absolute path (str)=')
args    = parser.parse_args()
halocat = args.cat[0]

Mcut    = 0.#5e10
print("Keeping only halos above "+str(Mcut)+" Msun")

#Load and filter
print("Loading and filtering")
f                = genfromtxt(halocat,names=True)
filter           = f['fMhires38']>0.
filter1          = f[filter]['hostHalo2']==0
filter2          = f[filter][filter1]['Mhalo4']>=Mcut
cathires         = f[filter]
cathiresnosub    = f[filter][filter1]
cathiresnosubMhi = f[filter][filter1][filter2] 
#need to treat IDS correctly as long int
IDS                 = loadtxt(halocat,usecols=(0,1),dtype={'names':('ID','Host'), 'formats':(int64,int64)})
subIDShires         = IDS[filter]
subIDShiresnosub    = IDS[filter][filter1]
subIDShiresnosubMhi = IDS[filter][filter1][filter2]

#write
print("Writing")
f1            = open(halocat+"_hires",'w')
f2            = open(halocat+"_main_hires",'w') #_Mcut"+str("{:.0e}".format(Decimal(Mcut))),'w')
print("File 1")
#for i in range(len(cathires)):
#    new_line=str(list(subIDShires[i])[0])+"\t"+str(list(subIDShires[i])[1])
#    for j in range(len(cathires[i])-2):
#        new_line+="\t"+str(list(cathires[i])[j+2])
#    new_line+="\n"
#    f1.write(new_line)
#f1.close()
print("File 2")
for i in range(len(cathiresnosubMhi)):
    new_line=str(list(subIDShiresnosubMhi[i])[0])+"\t"+str(list(subIDShiresnosubMhi[i])[1])
    for j in range(len(cathiresnosubMhi[i])-2):
        new_line+="\t"+str(list(cathiresnosubMhi[i])[j+2])
    new_line+="\n"
    f2.write(new_line)
f2.close()

print("Done")
