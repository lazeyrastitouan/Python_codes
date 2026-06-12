import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dfiles',type=str,nargs='*')
parser.add_argument('outpath',type=str,nargs=1)
args   = parser.parse_args()
files  = args.dfiles
out    = args.outpath[0]
print("Loading "+str(files))
outname = out+'/density_temperature_volume_in_Hcat.npz'
print('Outputting in '+outname)

#Need keys
f1   = np.load(files[0])
keys = list(f1.keys())
data = {k: [] for k in keys} 

for f in files:
    chunk = np.load(f)
    if list(chunk.keys()) != keys:
        exit("File "+f+ " does not have the correct keys")
    for k in list(chunk.keys()):
        data[k].append(chunk[k])

for k in list(chunk.keys()):
    data[k]=np.concatenate(data[k])

np.savez(outname,**data)

