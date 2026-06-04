from numpy import*
import matplotlib.pyplot as plt

verbose = False

posgal = load('/data/titouan/IllustrisTNG/300-1/fof_galaxy_pos.npy') #mpc/h
ngal   = len(posgal)
print('Number of galaxies: '+str(ngal))

voids  = genfromtxt('/data/titouan/IllustrisTNG/300-1/zobov-Voids_cat.txt',names=True,delimiter=' ')
posV   = array((voids['X'],voids['Y'],voids['Z3Mpch'])).T #mpc/h
Reff   = voids['R_effMpch'] #mpc/h
fRv    = [0,1]

#finding galaxies in voids

count = 0
for i in range(len(Reff)):
    j = 0
    while((j<len(posgal)) and (linalg.norm(posgal[j]-posV[i]) > Reff[i])):
        j+=1
    if j < len(posgal):
        count+=1
        posgal = delete(posgal,j,axis=0)

print('Number of galaxies inside voids: '+str(count))
print('Number of galaxies outside voids: '+str(len(posgal)))

fgalV    = count/ngal*100.

#takes too much memory
'''
dist     = [linalg.norm(pos-posV,axis=1) for pos in posgal]
dist     = array(dist)
print('dist computed')
inorout  = where((dist>=fRv[0]*Reff) & (dist<=fRv[1]*Reff),1,0) #1 if in, 0 if out
del dist
selected = where(any(inorout==1,axis=1))
print('selection done')
del inorout
fgalV    = shape(posgal[selected])[0]/shape(posgal)[0]*100.
'''

print('Fraction of galaxies in voids = '+str(fgalV)+'%')


