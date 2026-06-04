import numpy as np
import MAS_library as MASL
import Pk_library as PKL
import sys,os,glob
import argparse
import smoothing_library as SL

parser = argparse.ArgumentParser(description="description of routine. Example python Halo_overdensity.py 1 4 250-512-2fluid-selected-fiducial")   
parser.add_argument("first",   help="first realization number", type=int)
parser.add_argument("last",    help="last  realization number", type=int)
parser.add_argument("cosmo",   help="cosmology folder")
args  = parser.parse_args()
first, last, cosmo = args.first, args.last, args.cosmo


##################################### INPUT ############################################
#root1 = '/cobra/ptmp/titouan/baryon-cdm/2fluids'
#root2 = '/cobra/ptmp/titouan/baryon-cdm/2fluids'

root1 = '/cobra/ptmp/titouan/baryon-cdm/tests'                                                                                                

root2 = '/cobra/ptmp/titouan/baryon-cdm/tests'     

#root1 = '/cobra/ptmp/titouan/baryon-cdm'
#root2 = '/cobra/ptmp/titouan/baryon-cdm'

ptypesb  = [0]
ptypesc = [1]
grid    = 512
MAS     = 'CIC'
do_RSD  = False
axis    = 0
threads = 1

snapnum = 14
BoxSize = 250.0 #Mpc/h
Filter  = 'Top-Hat'
########################################################################################
    
if not(os.path.exists(cosmo)):  os.system('mkdir %s'%cosmo)

# do a loop over the different realizations
numbers = np.arange(first, last+1)
for i in numbers:

	folder_in  = '%s/%s/%d/'%(root1,cosmo,i)
	folder_out = '%s/%s/%d'%(root1,cosmo,i)
	if not(os.path.exists(folder_out)):  os.makedirs(folder_out)

	# compute the CDM+B density field
	snapshot = '%s/%s/%d/snapshot_%03d'%(root2,cosmo,i,snapnum)
	df_b = MASL.density_field_gadget(snapshot, ptypesb, grid, MAS, do_RSD, axis)
	df_b /= np.mean(df_b, dtype=np.float64);  df_b -= 1.0
	df_c = MASL.density_field_gadget(snapshot, ptypesc, grid, MAS, do_RSD, axis)
	df_c /= np.mean(df_c, dtype=np.float64);  df_c -= 1.0
	df_cb = df_b - df_c

	for smoothing in [10.0, 20.0, 30.0, 40.0]:

		# get files names
		name = 'combined_halos_mhalo_z0_nosubhalo'
		fout = '%s/delta_%.1f_%s.txt'%(folder_out, smoothing, name)
		if os.path.exists(fout):  continue

		data = np.loadtxt('%s/%s'%(folder_in,name))
		M, pos, Mc, Mb = data[:,2], data[:,4:7], data[:,7], data[:,8]
		pos = pos.astype(np.float32)/1e3 #Mpc/h

		#smooth
		W_k = SL.FT_filter(BoxSize, smoothing, grid, Filter, threads)
		df_cb_smoothed = SL.field_smoothing(df_cb, W_k, threads)       
		
		#interpolate on halo pos
		halos_overdensity = np.zeros(pos.shape[0], dtype=np.float32)
		MASL.CIC_interp(df_cb_smoothed, BoxSize, pos, halos_overdensity)        

		# save results to file
		np.savetxt(fout, np.transpose([M,pos[:,0],pos[:,1],pos[:,2],Mc,Mb,halos_overdensity]))

        
