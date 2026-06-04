#Formate catalogs for Mass_function.py (only main halos and M in 1st column) 

import numpy as np
from numpy import*
import matplotlib.pyplot as plt
import math as m
import matplotlib
import argparse

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# What to do
parser = argparse.ArgumentParser(description='This script formates catalogs for Mass_function.py (only main halos and M in 1st column). Arg: Path to halo catalog (str), halos catalog name (str): in AHF output format.')
parser.add_argument('path',nargs=1,type=str, help='Path to halo catalogs')

args = parser.parse_args()
path=args.path[0]

#Charger le catalogue correspondant
#halo/subhalo | M | x | y | z
catalog = loadtxt(path+"/200_snap_z4.000_halos_sorted_main_hires",usecols=(1,3,5,6,7))
Mcat=catalog[:,1]
halos = catalog[:,0]
    
f=open(path+"/200_snap_z4.000_halos_sorted_main_hires_nosubh","w")
for i in range(len(halos)):
    if halos[i]==0:
        f.write(str(Mcat[i])+"\t"+str(catalog[i,2])+"\t"+str(catalog[i,3])+"\t"+str(catalog[i,4])+"\n")
f.close()
