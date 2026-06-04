from mpi4py import MPI
import numpy as np
import MAS_library as MASL
import Pk_library as PKL
import sys,os,glob
import argparse
import smoothing_library as SL

###### MPI DEFINITIONS ######                                    
comm   = MPI.COMM_WORLD
nprocs = comm.Get_size()
myrank = comm.Get_rank()

parser = argparse.ArgumentParser(description="description of routine")   
parser.add_argument("first",   help="first realization number", type=int)
parser.add_argument("last",    help="last  realization number", type=int)
parser.add_argument("cosmo",   help="cosmology folder")
args  = parser.parse_args()
first, last, cosmo = args.first, args.last, args.cosmo


##################################### INPUT ############################################
root1 = '/projects/QUIJOTE/titouan/AHF'
root2 = '/projects/QUIJOTE/Snapshots'

#cosmologies = ['Om_p', 'Ob2_p', 'h_p', 'ns_p', 's8_p',
#               'Om_m', 'Ob2_m', 'h_m', 'ns_m', 's8_m',
#               'Mnu_p', 'Mnu_pp', 'Mnu_ppp',
#               'fiducial']

cosmologies = ['Mnu_ppp', 's8_p', 'fiducial']
realizations = 500

ptypes  = [1]
grid    = 512
MAS     = 'CIC'
do_RSD  = False
axis    = 0
threads = 1

snapnum = 4
BoxSize = 1000.0 #Mpc/h
Filter  = 'Top-Hat'
########################################################################################
    
if not(os.path.exists(cosmo)):  os.system('mkdir %s'%cosmo)

# do a loop over the different realizations
numbers = np.where(np.arange(first, last)%nprocs==myrank)[0]
numbers = np.arange(first, last)[numbers]
for i in numbers:

    folder_in  = '%s/%s/%d/'%(root1,cosmo,i)
    folder_out = '%s/%d'%(cosmo,i)
    if not(os.path.exists(folder_out)):  os.makedirs(folder_out)

    # compute the CDM+B density field
    snapshot = '%s/%s/%d/snapdir_%03d/snap_%03d'%(root2,cosmo,i,snapnum,snapnum)
    df_cb = MASL.density_field_gadget(snapshot, ptypes, grid, MAS, do_RSD, axis)
    df_cb /= np.mean(df_cb, dtype=np.float64);  df_cb -= 1.0

    for smoothing in [5.0, 10.0, 20.0, 30.0]:

        for suffix in ['14.25', '14.51', '14.78', '15.04', '15.3']:

            # get files names
            name = '200_snap_z0.000_M%s_sorted'%suffix
            fout = '%s/delta_%.1f_%s.txt'%(folder_out, smoothing, name)
            #fout = '%s/delta_%.1f_m_%s.txt'%(folder_out, smoothing, name)
            if os.path.exists(fout):  continue

            data = np.loadtxt('%s/%s'%(folder_in,name))
            M, pos, L, s, c = data[:,0], data[:,1:4], data[:,4], data[:,5], data[:,6]
            pos = pos.astype(np.float32)/1e3 #Mpc/h

            W_k = SL.FT_filter(BoxSize, smoothing, grid, Filter, threads)
            df_cb_smoothed = SL.field_smoothing(df_cb, W_k, threads)       

            halos_overdensity = np.zeros(pos.shape[0], dtype=np.float32)
            MASL.CIC_interp(df_cb_smoothed, BoxSize, pos, halos_overdensity)        

            # save results to file
            np.savetxt(fout, np.transpose([M,pos[:,0],pos[:,1],pos[:,2],halos_overdensity]))

        
