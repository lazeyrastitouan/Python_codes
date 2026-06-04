from numpy import*
import argparse

#input parameters
parser = argparse.ArgumentParser(description='Get halo position in a mass range from AHF catalog, save positions in pixel for a MockOb of Ng=512, L=4cMpc,  around a given input halo at xyz. Args: < halo catalog (str) >, < Mmin (float) >, < Mmax (float) >, < Nbins (int) >, < xyz (position of halo at center of MockoB, 3 floats) >')
parser.add_argument('cat',nargs=1,type=str, help='Halo catalog absolute path, str')
parser.add_argument('mmin',nargs=1,type=float, help='Minimum mass for halo selection in solar masses, float')
parser.add_argument('mmax',nargs=1,type=float, help='Maximum mass for halo selection in solar masses, float')
parser.add_argument('nbins',nargs=1,type=int,help="Number of mass bins, int")
parser.add_argument('xyz',nargs=3,type=float, help='x y z positions of halo at center of the MockOb in ckpc/h, 3 floats')
args  = parser.parse_args()
cat   = args.cat[0]
mmin  = args.mmin[0]
mmax  = args.mmax[0]
nbins = args.nbins[0]
xyz   = array(args.xyz[0:3])

h     = 0.6777
xyz   = xyz/h
Ng    = 512
Lbox  = 4000.0

# center of grid
centerpix = int(Ng/2)

#min and max xyz to consider for halos of the catalog to fall into the MockOb
xyzmax = xyz + Lbox/2.
xyzmin = xyz - Lbox/2.

#load catalog, keeping only halos with 100% high resolution DM particles
f      = genfromtxt(cat,usecols=(3,5,6,7,37),names="M, x, y, z, fMhires")
f['x'] = f['x']/h
f['y'] = f['y']/h
f['z'] = f['z']/h
f['M'] = f['M']/h
filter = f['fMhires']==1.0

#mass selection
filterM = (f[filter]['M']<=mmax) & (f[filter]['M']>=mmin)

#selectpositions are within the MockOb
filterx = (f[filter][filterM]['x']>=xyzmin[0]) & (f[filter][filterM]['x']<=xyzmax[0])
filtery = (f[filter][filterM][filterx]['y']>=xyzmin[1]) & (f[filter][filterM][filterx]['y']<=xyzmax[1])
filterz = (f[filter][filterM][filterx][filtery]['z']>=xyzmin[2]) & (f[filter][filterM][filterx][filtery]['z']<=xyzmax[2])
catfil  = f[filter][filterM][filterx][filtery][filterz]
Nh      = len(cat)

print('Keeping '+str(Nh)+' halos')

#get position in pixel in MockOb
catfil['x'] = (catfil['x']-xyz[0]).astype(int)+centerpix
catfil['y'] = (catfil['y']-xyz[1]).astype(int)+centerpix
catfil['z'] = (catfil['z']-xyz[2]).astype(int)+centerpix

#bin in mass and write
Mbins = logspace(log10(mmax),log10(mmin),nbins)
j=0
for i in range(len(Mbins)-1):
	edgelow = round(log10(Mbins[i+1]),2)
	edgeup  = round(log10(Mbins[i]),2)
	f = open(cat+'_Mmin_'+str(edgelow)+'_Mmax_'+str(edgeup)+'_pixelpos','w')
	while catfil['M'][j]>Mbins[i+1] and j<Nh:
		f.write(str(catfil['M'][j])+"\t"+str(int(catfil['x'][j]))+"\t"+str(int(catfil['y'][j]))+"\t"+str(int(catfil['z'][j]))+"\n")
		j+=1
	f.close()

