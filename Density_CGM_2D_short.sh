#!/bin/bash -l
#
#SBATCH --ntasks 1
#SBATCH -J Density_CGM
#SBATCH -o slurm.%J.out
#SBATCH -e slurm.%J.err
#SBATCH -p cosma-analyse
#SBATCH -A do019
#####SBATCH --exclusive
#SBATCH -t 06:00:00
#SBATCH --mail-type=END                          # notifications for job done & fail
#SBATCH --mail-user=titouan.lazeyras@unimib.it

module purge
module load cosma/2024

#activate the colibre venv
source /cosma/apps/do019/dc-laze1/venvs/colibre-workshop/bin/activate

python /cosma/home/do019/dc-laze1/Python_codes/Density_CGM_2D.py /cosma8/data/dp004/colibre/Runs/L0100N1504/Thermal/snapshots/colibre_0064/ colibre /cosma8/data/do019/dc-laze1/COLIBRE/L0100N1504/200_snap_z3.0_halos_sorted_main_hires_Mlow_12.5_Mhigh_13.5.txt >> /cosma8/data/do019/dc-laze1/COLIBRE/L0100N1504/z3/Density_CGM_2D.log 

cd /cosma8/data/do019/dc-laze1/COLIBRE/L0100N1504/z3/
for l in X Y Z
do
        mv densityS_Tmin1e+3_Tmax1e+5_in_0.0-0.2Rvir_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_0.0-0.2Rvir_"$l"_M12p7.npz
        mv densityS_Tmin1e+3_Tmax1e+5_in_0.2-0.5Rvir_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_0.2-0.5Rvir_"$l"_M12p7.npz
        mv densityS_Tmin1e+3_Tmax1e+5_in_0.5-0.7Rvir_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_0.5-0.7Rvir_"$l"_M12p7.npz
        mv densityS_Tmin1e+3_Tmax1e+5_in_0.7-1.0Rvir_"$l".npz densityS_Tmin1e+3_Tmax1e+5_in_0.7-1.0Rvir_"$l"_M12p7.npz
done

#python /cosma/home/do019/dc-laze1/Python_codes/Density_CGM_2D.py /cosma8/data/dp004/colibre/Runs/L0100N1504/Thermal/snapshots/colibre_0064/ colibre /cosma8/data/do019/dc-laze1/COLIBRE/L0100N1504/z3/200_snap_z3.0_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt >> /cosma8/data/do019/dc-laze1/COLIBRE/L0100N1504/z3/Density_CGM_2D_M12p2.log

# Run the program 32 times (on 28 cores).
#mpirun -np 28 ./parallel_tasks 0 31 "python ./hello.py %d"
