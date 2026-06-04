#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./galV.out.%j
#SBATCH -e ./galV.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J TgalV
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#
# for OpenMP:
#SBATCH --cpus-per-task=1
#
#Request 900 GB of main memory per node in units of MB:
#SBATCH --mem=900000
#
#SBATCH --mail-type=end
#SBATCH --mail-user=titouan.lazeyras@unimib.it
#
# Wall clock limit:
#SBATCH --time=48:00:00

cd /home/titouan/codes/Python_codes/
python3 galaxies_in_voids.py > galaxies_in_voids.log 
