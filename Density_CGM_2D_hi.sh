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
# Request 900 GB of main memory per node in units of MB:
#SBATCH --mem=900000
#
#SBATCH --mail-type=end
#SBATCH --mail-user=titouan.lazeyras@unimib.it
#
# Wall clock limit:
#SBATCH --time=24:00:00

module purge

python3 /home/titouan/codes/Python_codes/Density_CGM_2D_ckpc.py /data/titouan/COLIBRE/L0100N1504/z3/colibre_0064/ colibre /data/titouan/COLIBRE/L0100N1504/z3/200_snap_z3.0_halos_sorted_main_hires_Mlow_12.5_Mhigh_13.5.txt >> /data/titouan/COLIBRE/L0100N1504/z3/colibre_0064/Density_CGM_2D_M12p7.log

cd /data/titouan/COLIBRE/L0100N1504/z3/colibre_0064/
for l in X Y Z 
do
	mv densityS_Tmin1e+3_Tmax1e+5_in_50-160ckpc_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_50-160ckpc_"$l"_M12p7.npz
	mv densityS_Tmin1e+3_Tmax1e+5_in_160-250ckpc_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_160-250ckpc_"$l"_M12p7.npz
	mv densityS_Tmin1e+3_Tmax1e+5_in_250-350ckpc_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_250-350ckpc_"$l"_M12p7.npz
	mv densityS_Tmin1e+3_Tmax1e+5_in_350-700ckpc_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_350-700ckpc_"$l"_M12p7.npz
done

# Run the program 32 times (on 28 cores).
#mpirun -np 28 ./parallel_tasks 0 31 "python ./hello.py %d"
