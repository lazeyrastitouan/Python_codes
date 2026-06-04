#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./dust.out.%j
#SBATCH -e ./dust.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J dust
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

python3 /home/titouan/codes/Python_codes/dust.py 

# Run the program 32 times (on 28 cores).
#mpirun -np 28 ./parallel_tasks 0 31 "python ./hello.py %d"
