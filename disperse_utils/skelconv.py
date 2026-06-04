import numpy as np
from tqdm import tqdm
import math
from sklearn.utils import resample
from astropy.io import fits

from service.EAGLE import snapshot_highN as sn
from utils import geometry


def read_NDskl_ascii(fname):
    sections = ['[CRITICAL POINTS]', '[FILAMENTS]', '[CRITICAL POINTS DATA]', '[FILAMENTS DATA]']
    
    with open(fname, 'r') as file:
        lines = file.read().splitlines()
    
    idxs = []
    for i in range(len(sections)):
        idxs.append(lines.index(sections[i]))

    critical_points_section = lines[idxs[0] + 1 :idxs[1]]
    filaments_section = lines[idxs[1] + 2 :idxs[2]]

    #critical points: Dictionary to store data
    #print(critical_points_section)
    critical_points = {'n_points': int(critical_points_section[0])}

    index_start = 1 
    for i in range(critical_points['n_points']):
        #print(critical_points_section[index_start].split(' '))
        type, posx, posy, posz, value, pairID, boundary = critical_points_section[index_start].split(' ')
        critical_points[i] = {'type' : int(type), 'position': [float(posx),float(posy),float(posz)], 'value': float(value), 'boundary':int(boundary), 'pair_ID':int(pairID)}
        #print(pairID)
        #print(critical_points[i])
        nfil = int(critical_points_section[index_start +1])
        nfil += 2
        index_start += nfil


    #filaments : Dictionary to store filament data
    filaments = {'n_filaments': int(lines[idxs[1] + 1])}
    print('# Filaments: ',filaments['n_filaments'])

    #Read filaments data from section
    idx_filament = 0
    #i : n_filament
    #idx_filament : position of filament info in section
    for i in range(filaments['n_filaments']):
        n_points = int(filaments_section[idx_filament].split()[2])
        idx_start = idx_filament +1
        idx_stop = idx_start + n_points 
        #filament extracted points
        pos = []
        for x in filaments_section[idx_start : idx_stop]:
            tmp = x.split()
            pos_tmp = []
            for n in tmp: 
                pos_tmp.append(float(n))
            pos.append(pos_tmp)

        index_cp_1 = int(filaments_section[idx_filament].split()[0])
        index_cp_2 = int(filaments_section[idx_filament].split()[1])

        filaments['f_' + str(i)] = {'critical_point_1' : critical_points[index_cp_1], 'critical_point_2': critical_points[index_cp_2],
                            'data': pos}
        

        idx_filament = idx_stop
    
    return filaments, critical_points

# #source = {'simulation', 'field'}
# def compute_filaments_DensityProfile(fname, filaments, radii, scale = 5, source='simulation'):
    
    
#     for i in tqdm(range(filaments['n_filaments'])):
#         xx,yy,zz = np.array(filaments['f_'+ str(i)]['data']).T
    
#         #Rescale to match boxsize [cMpc]
#         xx = xx/scale
#         yy = yy/scale
#         zz = zz/scale

#         #Initialisation
#         length = 0.  
#         filaments['f_'+ str(i)]['radii'] = []
#         filaments['f_'+ str(i)]['n_density'] = []
#         filaments['f_'+ str(i)]['n_density_error'] = []
#         tmp = {}
#         for k in range(1, len(radii) +1):
#             tmp[str(radii[k])] = []

#         #Evaluate each filament segment
#         for j in range(len(xx)-1):
#             extreme_1 = np.array([xx[j], yy[j], zz[j]])
#             extreme_2 = np.array([xx[j+1], yy[j+1], zz[j+1]])

#             #Norm of segment
#             dist = np.linalg.norm(extreme_2 - extreme_1) 
#             length += dist
            
#             if source == 'simulation': 
#                 sim = sn.EagleSimulation(fname + '.0.hdf5')
#                 #Select a region to evaluate
#                 region = np.array([
#                         (min([xx[j], xx[j+1]]) - radii[k] ), (max([xx[j], xx[j+1]]) + radii[k] ),
#                         (min([yy[j], yy[j+1]]) - radii[k] ), (max([yy[j], yy[j+1]]) + radii[k] ),
#                         (min([zz[j], zz[j+1]]) - radii[k] ), (max([zz[j], zz[j+1]]) + radii[k] )
#                         ])

#                 gas_coord = sim.read_region(0,'Coordinates', region) #cMpc/h
#                 gas_coord /= sim.h #cMpc

#                 #Downsample
#                 gas_coord = resample(gas_coord, replace = False, n_samples= int(0.1 * len(gas_coord)), random_state=42)

#             if source == 'field':
#                 hdul = fits.open(fname)
#                 data_fits = hdul[0].data
#                 z,y,x = data_fits.shape
#                 z,y,x = np.mgrid[0:z,0:y,0:x]
                
#                 #Coordinates
#                 gas_coord = np.vstack((z.flat,y.flat,x.flat)).T
#                 gas_coord/=scale

#                 #Density values
#                 dd = np.ravel(data_fits)

#             #Actual evaluation of radii
#             number_points = 0 
#             eval = np.full((len(gas_coord)), False)
#             for k in range(1, len(radii)+1):
#                 mask = []

#                 #Evaluate if gas particle inside a cylinder centered in the segment
#                 for q in tqdm(gas_coord): 
#                     mask.append(geometry.isInsideCylinder(extreme_1, extreme_2, radius = radii[k], point_to_evaluate = q))

#                 mask = np.array(mask)
                
#                 #If eval had already True value (point inside a inner cylinder), mask has to not account for it
#                 #Mask : cylindrical shell 
#                 mask[np.where(eval == True)] = False
                
#                 #Eval is updated with the true values of mask: points for this cylinders and inner ones
#                 #Eval : cylinder shell + inner cylinders 
#                 eval[np.where(mask == True)] = True

                
#                 #Compute the number density in the cylindrical shell
#                 if source=='simulation':
#                     number_points = np.sum(mask)
            
#                     volume = 2 * math.pi * dist * (radii[k]**2 - radii[k-1]**2)

#                     tmp[str(radii[k])].append(number_points/volume)


#                 elif source== 'field':
#                     tmp[str(radii[k])].append(np.mean(dd[mask])) #TODO: check mean and error
                                     


#         #Store results
#         for k in range(1, len(radii)+1):
#             filaments['f_'+ str(i)]['radii'].append(radii[k])
#             filaments['f_'+ str(i)]['n_density'].append(np.median(tmp[str(radii[k])]))
#             filaments['f_'+ str(i)]['n_density_error'].append(np.quantile(tmp[str(radii[k])], 0.95) - np.quantile(tmp[str(radii[k])], 0.05))
#         filaments['f_'+ str(i)]['length'] = length
    
#     return filaments



    

