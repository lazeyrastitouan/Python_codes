#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./TvsD.out.%j
#SBATCH -e ./TvsD.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J Trho
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

#python3 TvsD.py /data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/ colibre --Hcat /data/titouan/COLIBRE/L0050N0752/z3p5/200_snap_z3.5_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in --annulus 100 200
#mv /data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/TvsD_in.pdf /data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/TvsD_in_100_200_M12-12p5.pdf
#python3 TvsD.py --Hcat /data/titouan/HELLO/CDM/g9.20e12/g9.20e12.00128.z3.596.AHF_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in --annulus 50 160 /data/titouan/HELLO/CDM/g9.20e12/ hello
#mv /data/titouan/HELLO/CDM/g9.20e12/TvsD_in.pdf /data/titouan/HELLO/CDM/g9.20e12/TvsD_in_50-160_M12-12p5.pdf

#python3 TvsD.py --Hcat /data/titouan/COLIBRE/L0050N0752/z3p5/200_snap_z3.5_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in --annulus 50 160 /data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/ colibre
#mv /data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/TvsD_in.pdf /data/titouan/COLIBRE/L0050N0752/z3p5/colibre_0060/TvsD_in_50_160_M12-12p5.pdf

#H10 zm2
python3 TvsD.py --Hcat /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/AHF_004_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in  --annulus 50 160 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/ dalya
mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/TvsD_in.pdf /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/TvsD_in_50_160_M12-12p5.pdf

python3 TvsD.py --Hcat /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/AHF_004_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in  --annulus 100 200 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/ dalya
mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/TvsD_in.pdf /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo10/sftresh_cool_low_fb_boost/snapdir_004/TvsD_in_100_200_M12-12p5.pdf

#H1312 zm1p5
python3 TvsD.py --Hcat /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/AHF_001_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in  --annulus 50 160 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/ dalya
mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in.pdf /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in_50_160_M12-12p5.pdf

python3 TvsD.py --Hcat /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/AHF_001_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in  --annulus 100 200 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/ dalya
mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in.pdf /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom1p5_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in_100_200_M12-12p5.pdf

#H1312 zm2
python3 TvsD.py --Hcat /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/AHF_001_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in  --annulus 50 160 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/ dalya
mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in.pdf /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in_50_160_M12-12p5.pdf

python3 TvsD.py --Hcat /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/AHF_001_allparts/200_snap_z3.500_halos_sorted_main_hires_Mlow_12.0_Mhigh_12.5.txt --InOut in  --annulus 100 200 /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/ dalya
mv /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in.pdf /data/titouan/Eagle_cosmo_gizmo/200_512/CDM/zoom2_hydro_fire2_lvl7/Halo1312/sftresh_cool_low_fb_boost/snapdir_001/TvsD_in_100_200_M12-12p5.pdf
