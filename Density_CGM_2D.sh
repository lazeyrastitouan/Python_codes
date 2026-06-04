#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./Density_CGM.out.%j
#SBATCH -e ./Density_CGM.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J DCGM2D
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#
# for OpenMP:
#SBATCH --cpus-per-task=1
#
# Request 950 GB of main memory per node in units of MB:
#SBATCH --mem=950000
#
#SBATCH --mail-type=end
#SBATCH --mail-user=titouan.lazeyras@unimib.it
#
# Wall clock limit:
#SBATCH --time=24:00:00

module purge

python3 /home/titouan/codes/Python_codes/Density_CGM_2D_ckpc_memory_optimized.py /data/titouan/COLIBRE/L0050N1504/z3p5/colibre_0060/ colibre /data/titouan/COLIBRE/L0050N1504/z3p5/200_snap_z3.5_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt >> /data/titouan/COLIBRE/L0050N1504/z3p5/colibre_0060/Density_CGM_2D_M12p1.log

cd /data/titouan/COLIBRE/L0050N1504/z3p5/colibre_0060/
for l in X Y Z 
do
	mv density_volume_Tmin1e+3_Tmax1e+5_in_100-200ckpc_"$l".npz density_volume_Tmin1e+3_Tmax1e+5_in_100-200ckpc_"$l"_M12p1.npz
	#mv density_volume_Tmin1e+3_Tmax1e+5_in_200-300ckpc_"$l".npz density_volume_Tmin1e+3_Tmax1e+5_in_200-300ckpc_"$l"_M12p1.npz
done

# Run the program 32 times (on 28 cores).
#mpirun -np 28 ./parallel_tasks 0 31 "python ./hello.py %d"
