#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./coldgas.out.%j
#SBATCH -e ./coldgas.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J coldgas
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
#SBATCH --time=48:00:00

cd /home/titouan/codes/Python_codes/

python3 ColdGas.py /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_007 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/AHF_007_allparts/200_snap_z3.000_halos_sorted_main_hires_Mlow_12.5_Mhigh_13.1.txt dalya --snappath2 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_002/ --Hcat2 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/AHF_002_allparts/200_snap_z3.000_halos_sorted_main_hires_Mlow_12.5_Mhigh_13.1.txt
mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/fV_fcold_in_annulickpc_z3.0.npz /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/fV_fcold_in_annulickpc_z3.0_M12p5_13p1.npz

#python3 ColdGas.py /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_007 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/AHF_007_allparts/200_snap_z3.000_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt dalya 
#mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/fV_fcold_in_annulickpc_z3.0.npz /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/fV_fcold_in_annulickpc_z3.0_M12_12p5.npz
