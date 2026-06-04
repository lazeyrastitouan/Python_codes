#DIVISE HALO CATALOG EN SOUS CATALOGUES DANS UN BIN DE MASSE
 
from numpy import*
import sys

if len(sys.argv[1:]) < 1 :
        print("Usage : python haloselection.py z \n")
        exit(1)

z=float(sys.argv[1])

#input parameters
NumOfFile = 160

#edges des bins 
Medges=loadtxt("bins_20p.txt",usecols=(0,))[1:4]
M=loadtxt("bins_20p.txt",usecols=(1,))[0:3]
lm=len(M)
	
fileNum =1
while fileNum <= NumOfFile :
	#print fileNum
	red='{0:.1f}'.format(z)	
	#Charger le catalogue correspondant
	catalog = loadtxt(str(fileNum)+"/200_2400Mpc_z"+str(red)+"_sorted",usecols=(1,3,5,6,7))
	Mcat=catalog[:,1]
	halos = catalog[:,0]
	
	jj=0
	j=0
	while jj < lm:
		fsubcat=open(str(fileNum)+"/200_2400Mpc_z"+str(red)+"_M"+str(M[jj])+"_sorted","w")
		while j <= shape(halos)[0]-1 and Mcat[j]>= Medges[jj]:
			if halos[j] == -1:
				#write line  
				new_line=str(Mcat[j])+"\t"+str(catalog[j,2])+"\t"+str(catalog[j,3])+"\t"+str(catalog[j,4])+"\n"
				fsubcat.write(new_line)
			j=j+1
		fsubcat.close()
		jj=jj+1

	print("Done "+str(fileNum))
	
	fileNum = fileNum+1
