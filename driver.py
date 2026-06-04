# This script is used to submit jobs to compute Pk from different realizations/cosmologies
import numpy as np
import sys,os

################################## INPUT #############################################
realizations = 500
step         = 100  #number of realizations each cpu will do
offset       = 0    #the count will start from offset
######################################################################################


# do a loop over the different cosmologies
for cosmo in ['Om_p/', 'Ob2_p/', 'h_p/', 'ns_p/', 's8_p/',          
              'Om_m/', 'Ob2_m/', 'h_m/', 'ns_m/', 's8_m/',           
              'Mnu_p/', 'Mnu_pp/', 'Mnu_ppp/', 'fiducial/']:

    nodes = realizations/step
    if nodes==0:  nodes=1

    # do a loop over the different realizations
    for i in xrange(nodes):    

        a = """#!/bin/bash
#SBATCH -J Pk
#SBATCH --exclusive        
#SBATCH -t 05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#######SBATCH --ntasks-per-node=48
#SBATCH --partition=general
#SBATCH --export=ALL

source ~/Pylians
        
srun -n 3 --mpi=pmi2 python Halos_overdensity.py %d %d %s\n
        """%(i*step+offset, (i+1)*step+offset, cosmo)

        # create the script.sh file, execute it and remove it
        f = open('script.sh','w');  f.write(a);  f.close()
        os.system('sbatch script.sh');  os.system('rm script.sh')
