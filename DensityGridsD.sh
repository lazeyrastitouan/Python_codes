#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./DensityGrids.out.%j
#SBATCH -e ./DensityGrids.err.%j
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

for m in 12.0 12.5
do
	python3 MedianDensityGrid.py -simtype dalya -subbox 2000. -ngrid 256 -tmin 1e3 -tmax 1e5 -mmin $m
#	python3 MedianDensityGrid.py -simtype colibre -subbox 2000. -ngrid 256 -tmin 1e3 -tmax 1e5 -mmin $m
done
# Run the program 32 times (on 28 cores).
#mpirun -np 28 ./parallel_tasks 0 31 "python ./hello.py %d"
