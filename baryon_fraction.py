import numpy as np
from numpy import*
import matplotlib.pyplot as plt
import math as m
import matplotlib
import argparse


particles_file = np.loadtxt('/cobra/ptmp/titouan/baryon-cdm/test-250Mpc-512-samesetup-floor-eq-dmsof/example.AHF_particles')

#print particles_file
parid = particles_file[:,0]
partype = particles_file[:,1]

#print len(partype)
#if len(partype) != 1:
#    print partype

for i in partype:
    if i == 0:
        print partype
        i = i + 1
