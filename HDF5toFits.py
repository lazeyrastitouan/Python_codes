import h5py
import numpy as np
from astropy.io import fits

f = h5py.File('/data/titouan/ObjectID.hdf5')

filamentID = f['filamentID']

filamentID = np.swapaxes(filamentID, 0, 2)

filaments = fits.PrimaryHDU(filamentID)
filaments.writeto('/data/titouan/ObjectID.fits', overwrite=True)


f.close()
