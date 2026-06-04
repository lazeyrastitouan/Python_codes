import numpy as np
import h5py
#from pyread_eagle import EagleSnapshot

def calc2DRadialDistance(sideX, sideY, centreX, centreY):
    """
    calculates the 2D radial distance from a given point for the whole grid
    INPUT:  sideX: the x-side length in cells (int)
            sideY: the y-side length in cells (int)
            centreX:    the x-coordinate in cells of the centre (int)
            centreY:    the y-coordinate in cells of the centre (int)
    RETURN: the 2D radial distance to the centre of each cell (numpy float array)
    """
    
    # need to handle the cases where sideX != sideY
    maxDiff = max(sideX, sideY)
    sideXmod = sideX + np.clip((sideY-sideX), 0, maxDiff)
    sideYmod = sideY + np.clip((sideX-sideY), 0, maxDiff)

    xDist = np.arange(sideXmod) - centreX
    xDist2D = np.tile(xDist, (sideXmod,1))
    xDist2D = xDist2D.T

    yDist = np.arange(sideYmod) - centreY
    yDist2D = np.tile(yDist, (sideYmod,1))

    radDist2D = np.sqrt(xDist2D**2 + yDist2D**2)
    
    return radDist2D[0:sideX, 0:sideY]

def calc3DRadialDistance(sideX, sideY, sideZ, centreX, centreY, centreZ):
    """
    calculates the 3D radial distance from a given point for the whole grid
    INPUT:  sideX: the x-side length in cells (int)
            sideY: the y-side length in cells (int)
            sideZ: the z-side length in cells (int)
            centreX:    the x-coordinate in cells of the centre (int)
            centreY:    the y-coordinate in cells of the centre (int)
            centreZ:    the z-coordinate in cells of the centre (int)
    RETURN: the 3D radial distance to the centre of each cell (numpy float array)
    """
    xDist = np.arange(sideX) - centreX
    xDist3D = np.tile(xDist, (sideX, sideX,1))
    xDist3D = np.swapaxes(xDist3D, 0,2)
    
    yDist = np.arange(sideY) - centreY
    yDist3D = np.tile(yDist, (sideY, sideY,1))
    yDist3D = np.swapaxes(yDist3D, 1,2)

    zDist = np.arange(sideZ) - centreZ
    zDist3D = np.tile(zDist, (sideZ, sideZ, 1))

    radDist3D = np.sqrt(xDist3D**2 + yDist3D**2 + zDist3D**2)
    
    return radDist3D

def radialProfile(valueGrid, binSize, radDistGrid):
    """
    Calculates the radial profile (median) of valueGrid. The dimension of the valueGrid can be anything.
    INPUT:  valueGrid: a grid holding values that will be radially averaged (float numpy array)
            binSize: the size of the bins in cells (int)
            radDistGrid: grid specifying all the distances to the centre in units of cells for each cell of valueGrid (float numpy array)
            Note: np.shape(valueGrid) = np.shape(radDistGrid)
    RETURN: bins: 1D array holding the bin edges (1D float numpy array)
            radialProf: 1D array holding the radial medians (1D float numpy array)
            perc16: 1D array holding the 16th percentiles of the radial bins
            perc84: 1D array holding the 84th percentiles of the radial bins
            Note: len(bins) = len(radialProf) + 1 and len(radialProf) = len(perc16) = len(perc84)
    """
    bins = np.arange(0, int(np.amax(radDistGrid)), binSize)
    radialProf = np.zeros_like(bins[:-1], dtype=float)
    perc16 = np.zeros_like(bins[:-1], dtype=float)
    perc84 = np.zeros_like(bins[:-1], dtype=float)

    for i in range(0, np.shape(bins)[0] - 1):
        Rfilter = np.logical_and(bins[i] <= radDistGrid, radDistGrid < bins[i + 1])
        rSelection = valueGrid[Rfilter]
        if len(rSelection) == 0:
            continue
        else:
            radialProf[i] = np.nanmedian(rSelection)
            perc16[i] = np.nanpercentile(rSelection, 16)
            perc84[i] = np.nanpercentile(rSelection, 84)

    return bins, radialProf, perc16, perc84

def radialProfile_log(valueGrid, binnum, radDistGrid):
    """
    Calculates the radial profile (median) of valueGrid. The dimension of the valueGrid can be anything.
    INPUT:  valueGrid: a grid holding values that will be radially averaged (float numpy array)
            binnum: the number of bins evenly spaced in log (int)
            radDistGrid: grid specifying all the distances to the centre in units of cells for each cell of valueGrid (float numpy array)
            Note: np.shape(valueGrid) = np.shape(radDistGrid)
    RETURN: bins: 1D array holding the bin edges (1D float numpy array)
            radialProf: 1D array holding the radial medians (1D float numpy array)
            perc16: 1D array holding the 16th percentiles of the radial bins
            perc84: 1D array holding the 84th percentiles of the radial bins
            Note: len(bins) = len(radialProf) + 1 and len(radialProf) = len(perc16) = len(perc84)
    """
    bins = np.logspace(np.log10(1.0), np.log10(np.amax(radDistGrid)), binnum)
    radialProf = np.zeros_like(bins[:-1], dtype=float)
    perc16 = np.zeros_like(bins[:-1], dtype=float)
    perc84 = np.zeros_like(bins[:-1], dtype=float)

    for i in range(0, np.shape(bins)[0] - 1):
        Rfilter = np.logical_and(bins[i] <= radDistGrid, radDistGrid < bins[i + 1])
        rSelection = valueGrid[Rfilter]
        if len(rSelection) == 0:
            continue
        else:
            radialProf[i] = np.nanmedian(rSelection)
            perc16[i] = np.nanpercentile(rSelection, 16)
            perc84[i] = np.nanpercentile(rSelection, 84)
    
    return bins, radialProf, perc16, perc84 

def radialProfileMean(valueGrid, binSize, radDistGrid, threshold=0):
    """
    Calculates the radial profile (mean) of valueGrid. The dimension of the valueGrid can be anything.
    INPUT:  valueGrid: a grid holding values that will be radially averaged (float numpy array)
            binSize: the size of the bins in cells (int)
            radDistGrid:    grid specifying all the distances to the centre in units of cells for each cell of valueGrid (float numpy array)
            threshold: SB threshold for pixels to be considered in the profile (float)
            Note: np.shape(valueGrid) = np.shape(radDistGrid)
    RETURN: bins: 1D array holding the bin edges (1D float numpy array)
            radialProf: 1D array holding the radial means (1D float numpy array)
            errProf: 1D array holding the error on the means (1D float numpy array) 
            Note: len(bins) = len(radialProf) + 1 and len(radialProf) = len(errProf)
    """
    bins = np.arange(0, int(np.amax(radDistGrid)), binSize)
    radialProf = np.zeros_like(bins[:-1], dtype=float)
    errProf = np.zeros_like(bins[:-1], dtype=float)

    valueGrid = np.where(valueGrid<threshold,np.nan,valueGrid)

    for i in range(0, np.shape(bins)[0] - 1):
        Rfilter = np.logical_and(bins[i] <= radDistGrid, radDistGrid < bins[i + 1])
        rSelection = valueGrid[Rfilter]
        if len(rSelection) == 0:
            continue
        else:
            radialProf[i] = np.nanmean(rSelection)
            errProf[i] = np.nanstd(rSelection)/np.sqrt(len(rSelection)-1.)
   
    return bins, radialProf, errProf

def radialProfileMean_log(valueGrid, binnum, radDistGrid):
    """
    Calculates the radial profile (mean) of valueGrid. The dimension of the valueGrid can be anything.
    INPUT:  valueGrid: a grid holding values that will be radially averaged (float numpy array)
            binnum: the number of bins evenly spaced in log (int)
            radDistGrid:    grid specifying all the distances to the centre in units of cells for each cell of valueGrid (float numpy array)
            threshold: SB threshold for pixels to be considered in the profile (float)
            Note: np.shape(valueGrid) = np.shape(radDistGrid)
    RETURN: bins: 1D array holding the bin edges (1D float numpy array)
            radialProf: 1D array holding the radial means (1D float numpy array)
            errProf: 1D array holding the error on the means (1D float numpy array) 
            Note: len(bins) = len(radialProf) + 1 and len(radialProf) = len(errProf)
    """
    bins = np.logspace(np.log10(1), np.log10(np.amax(radDistGrid)), binnum)
    radialProf = np.zeros_like(bins[:-1], dtype=float)
    errProf = np.zeros_like(bins[:-1], dtype=float)

    for i in range(0, np.shape(bins)[0] - 1):
        Rfilter = np.logical_and(bins[i] <= radDistGrid, radDistGrid < bins[i + 1])
        rSelection = valueGrid[Rfilter]
        if len(rSelection) == 0:
            continue
        else:
            radialProf[i] = np.nanmean(rSelection)
            errProf[i] = np.nanstd(rSelection)/np.sqrt(len(rSelection)-1.)

    return bins, radialProf, errProf

def radialProfile_av(valueGrid, binSize, radDistGrid, weightsGrid=None):
    """
    Calculates the radial profile (average) of valueGrid. The dimension of the valueGrid can be anything.
    INPUT:  valueGrid: a grid holding values that will be radially averaged (float numpy array)
            binSize: the size of the bins in cells (int)
            radDistGrid:    grid specifying all the distances to the centre in units of cells for each cell of valueGrid (float numpy array)
            weightsGrid:    the weights to use while averaging, shape needs to be broadcastable to to valueGrid, default is None
            Note: np.shape(valueGrid) = np.shape(radDistGrid)
    RETURN: bins: 1D array holding the bin edges (1D float numpy array)
                radialProf: 1D array holding the radial averages (1D float numpy array)
            stdDev: standard deviation of the average
            Note: len(bins) = len(radialProf) + 1 and len(radialProf) = len(stdDev)
    """
    bins = np.arange(0, int(np.amax(radDistGrid)), binSize)
    radialProf = np.zeros_like(bins[:-1], dtype=float)
    stdDev = np.zeros_like(bins[:-1], dtype=float)

    for i in range(0, np.shape(bins)[0] - 1):
        Rfilter = np.logical_and(bins[i] <= radDistGrid, radDistGrid < bins[i + 1])
        nanFilter = ~np.isnan(valueGrid)
        rSelection = valueGrid[np.logical_and(Rfilter, nanFilter)]
        if type(weightsGrid) == type(None): 
            wSelection = None
        else: 
            wSelection = weightsGrid[Rfilter]
        if len(rSelection) == 0:
            continue
        else:
            radialProf[i] = np.average(rSelection, weights=wSelection)
            stdDev[i] = np.std(rSelection)
    
    return bins, radialProf, stdDev

def profileFromScatter(variableA, variableB, binsVarA, avType='median'):
    """
    generates profile of variable B as a function of variable A
    INPUT:  variableA: numpy array holding the values of variable A
            variableB: numpy array holding the values of variable B
            binsVarA: numpy array holding the bin-egdes for variable A
            avType: string describing how to calculate the profile, options: median, mean
    RETURN: profile: the profile of variable B as a function of variable A
            spreadMin: 25th percentile or mean - std
            spreadMax: 75th percentile or mean + std
            binsVarA[:-1]: the bins of variable A adjusted to match the profile
    """

    numBins = np.shape(binsVarA)[0] - 1
    profile = np.zeros_like(binsVarA[:-1], dtype=float)
    spreadMin = np.zeros_like(profile, dtype=float)
    spreadMax = np.zeros_like(profile, dtype=float)

    for n in range(0, numBins):
        mask = np.logical_and(binsVarA[n] <= variableA, variableA < binsVarA[n + 1]) #define mask for the individual bins
        selection = variableB[mask]
        if len(selection) == 0: continue

        #populate the profile
        if avType == 'median':
            profile[n] = np.nanmedian(selection)
            spreadMin[n] = np.nanpercentile(selection, 25)
            spreadMax[n] = np.nanpercentile(selection, 75)
        elif avType == 'mean':
            profile[n] = np.nanmean(selection)
            spreadMin[n] = profile[n] - np.sqrt(np.nanvar(selection))
            spreadMax[n] = profile[n] + np.sqrt(np.nanvar(selection))

        binsVarA[n] += (binsVarA[n + 1] - binsVarA[n]) * 0.5

    return profile, spreadMin, spreadMax, binsVarA[:-1]


def stackProfiles_rVir(keywordList, stackingDictData, stackingDictBins, rVirDict, keywordToGnDict, binSize, rVirLimit, resolution):
    """
    generates stacked, rescaled by the virial radius profiles
    INPUT:  keywordList: list holding the identifiers of the various profiles to be stacked (strings)
            stackingDictData: dictionary holding the profiles to be stacked, keywords are in the keywordList (strings)
            stackingDictBins: dictionary holding the bins of the profiles, same keywords as stackingDictData (strings)
            rVirDict: dictionary holding the virial radii used for the rescaling, the keywords refer to the halo groupnumbers (strings)
            keywordToGnDict: dictionary linking groupnumber to the keywords (strings)
            binSize: binSize of the stacked profile in units of virialRadii (float)
            rVirLimit: max value to where the profile goes out to (float)
            resolution: size of cells (the units of the bins) in the same units as rVir (normally kpc)
    RETURN: rescaledBins: the bins corresponding to the stacked profile (units of rVir), 1D float numpy array
            median: the stacked profile 1D float numpy array (uses the median)
            perc25: 25th percentile of the stacked profile 1D float numpy array
            perc75: 75th percentile of the stacked profile 1D float numpy array
            medianRvir: the median virial radius of all the stacked profiles
    """
    # initialze the necessary variables
    binsRescaled = {}
    virialRadii = np.array([])
    median = np.array([])
    perc25 = np.array([])
    perc75 = np.array([])

    # cycle through the keywordList/List
    for keyword in keywordList:
        # load the virial radii
        gn = keywordToGnDict[keyword]
        virialRadii = np.append(virialRadii, rVirDict[gn])
        binsRescaled[keyword] = (stackingDictBins[keyword] * resolution) / rVirDict[gn]

    # initialize rescaled bins array in units of virial radii
    rescaledBins = np.arange(0, rVirLimit / binSize) * binSize + binSize

    # cycle through the rescaledBins and stack the profiles of each halo in the appropriate bins
    for dist in rescaledBins:
        # re-initialize selection for each bin
        selection = np.array([])
        for keyword in keywordList:
            # re-define masks for each bin and keyword
            mask = np.logical_and(dist - binSize <= binsRescaled[keyword], binsRescaled[keyword] < dist)
            selection = np.append(selection, stackingDictData[keyword][mask])
        if len(selection) == 0:
            median = np.append(median, 1e-20)
            perc25 = np.append(perc25, 1e-20)
            perc75 = np.append(perc75, 1e-20)
        else:
            median = np.append(median, np.nanmedian(selection))
            perc25 = np.append(perc25, np.nanpercentile(selection, 25))
            perc75 = np.append(perc75, np.nanpercentile(selection, 75))

    medianRvir = np.median(virialRadii)

    return rescaledBins, median, perc25, perc75, medianRvir

def stackProfilesAv_rVir(keywordList, stackingDictData, stackingDictBins, rVirDict, keywordToGnDict, binSize, rVirLimit, resolution):
    """
    generates stacked, rescaled by the virial radius profiles
    INPUT:  keywordList: list holding the identifiers of the various profiles to be stacked (strings)
            stackingDictData: dictionary holding the profiles to be stacked, keywords are in the keywordList (strings)
            stackingDictBins: dictionary holding the bins of the profiles, same keywords as stackingDictData (strings)
            rVirDict: dictionary holding the virial radii used for the rescaling, the keywords refer to the halo groupnumbers (strings)
            keywordToGnDict: dictionary linking groupnumber to the keywords (strings)
            binSize: binSize of the stacked profile in units of virialRadii (float)
            rVirLimit: max value to where the profile goes out to (float)
            resolution: size of cells (the units of the bins) in the same units as rVir (normally kpc)
    RETURN: rescaledBins: the bins corresponding to the stacked profile (units of rVir), 1D float numpy array
            median: the stacked profile 1D float numpy array (uses the median)
            perc25: 25th percentile of the stacked profile 1D float numpy array
            perc75: 75th percentile of the stacked profile 1D float numpy array
            medianRvir: the median virial radius of all the stacked profiles
    """
    # initialze the necessary variables
    binsRescaled = {}
    virialRadii = np.array([])
    mean = np.array([])
    stdDev = np.array([])

    # cycle through the keywordList/List
    for keyword in keywordList:
        # load the virial radii
        gn = keywordToGnDict[keyword]
        virialRadii = np.append(virialRadii, rVirDict[gn])
        binsRescaled[keyword] = (stackingDictBins[keyword] * resolution) / rVirDict[gn]

    # initialize rescaled bins array in units of virial radii
    rescaledBins = np.arange(0, rVirLimit / binSize) * binSize + binSize

    # cycle through the rescaledBins and stack the profiles of each halo in the appropriate bins
    for dist in rescaledBins:
        # re-initialize selection for each bin
        selection = np.array([])
        for keyword in keywordList:
            # re-define masks for each bin and keyword
            mask = np.logical_and(dist - binSize <= binsRescaled[keyword], binsRescaled[keyword] < dist)
            selection = np.append(selection, stackingDictData[keyword][mask])
        if len(selection) == 0:
            mean = np.append(mean, 1e-20)
            stdDev = np.append(stdDev, 1e-20)
        else:
            mean = np.append(mean, np.nanmean(selection))
            stdDev = np.append(stdDev, np.nanstd(selection))

    medianRvir = np.median(virialRadii)

    return rescaledBins, mean, stdDev, medianRvir

def stackProfiles(keywordList, stackingDictData, stackingDictBins, keywordToGnDict, binSize, cellLimit):
    """
    generates stacked profile (median) of the profiles stored in dictionary passed to the function
    INPUT:  keywordList: list holding the identifiers of the various profiles to be stacked (strings)
            stackingDictData: dictionary holding the profiles to be stacked, keywords are in the keywordList (strings)
            stackingDictBins: dictionary holding the bin edges of the profiles, same keywords as stackingDictData (strings)
            keywordToGnDict: dictionary linking groupnumber to the keywords (strings)
            binSize: binSize of the stacked profile in units of cells (float)
            cellLimit: max value to where the profile goes out to (float)
            len(stackingDictData[keyword]) = len(stackingDictBins[keyword]) - 1 
    RETURN: stackBins: the bin edges corresponding to the stacked profile (units of cells), 1D float numpy array
            median: the stacked, median profile; 1D float numpy array
            perc25: 25th percentile of the stacked profile; 1D float numpy array
            perc75: 75th percentile of the stacked profile 1D float numpy array
            len(median) = len(perc25) = len(perc75) = len(stackBins) - 1 
    """
    # initialze the necessary variables
    bins = {}

    # cycle through the keywordList/List
    for keyword in keywordList:
        #gn = keywordToGnDict[keyword]
        bins[keyword] = stackingDictBins[keyword] 

    # define the bins for the stack
    stackBins = np.arange(0,cellLimit + 1, binSize, dtype=float)
    # initialize the median, perc25 and perc75 arrays
    median = np.zeros_like(stackBins[:-1], dtype=float)
    perc25 = np.zeros_like(stackBins[:-1], dtype=float)
    perc75 = np.zeros_like(stackBins[:-1], dtype=float)

    # cycle through the rescaledBins and stack the profiles of each halo in the appropriate bins
    for i in range(0, len(median)):
        # re-initialize selection and for each bin
        selection = np.array([])
        for keyword in keywordList:
            rFilter = np.logical_and(np.greater_equal(bins[keyword][:-1], stackBins[i]), np.less(bins[keyword][:-1], stackBins[i+1]))
            selection = np.append(selection, stackingDictData[keyword][rFilter])
        if len(selection) == 0:
            continue
        else:
            median[i] = np.nanmedian(selection)
            perc25[i] = np.nanpercentile(selection, 25)
            perc75[i] = np.nanpercentile(selection, 75)
    
    return stackBins, median, perc25, perc75

def stackProfilesAv(keywordList, stackingDictData, stackingDictBins, keywordToGnDict, binSize, cellLimit):
    """
    generates stacked profile (median) of the profiles stored in dictionary passed to the function
    INPUT:  keywordList: list holding the identifiers of the various profiles to be stacked (strings)
            stackingDictData: dictionary holding the profiles to be stacked, keywords are in the keywordList (strings)
            stackingDictBins: dictionary holding the bin edges of the profiles, same keywords as stackingDictData (strings)
            keywordToGnDict: dictionary linking groupnumber to the keywords (strings)
            binSize: binSize of the stacked profile in units of cells (float)
            cellLimit: max value to where the profile goes out to (float)
            len(stackingDictData[keyword]) = len(stackingDictBins[keyword]) - 1 
    RETURN: stackBins: the bin edges corresponding to the stacked profile (units of cells), 1D float numpy array
            mean: the stacked, mean profile; 1D float numpy array
            stdDev: standard deviation of the stacked profile; 1D float numpy array
            len(median) = len(perc25) = len(perc75) = len(stackBins) - 1 
    """
    # initialze the necessary variables
    bins = {}

    # cycle through the keywordList/List
    for keyword in keywordList:
        #gn = keywordToGnDict[keyword]
        bins[keyword] = stackingDictBins[keyword] 

    # define the bins for the stack
    stackBins = np.arange(0,cellLimit + 1, binSize, dtype=float)
    # initialize the median, perc25 and perc75 arrays
    mean = np.zeros_like(stackBins[:-1], dtype=float)
    stdDev = np.zeros_like(stackBins[:-1], dtype=float)

    # cycle through the rescaledBins and stack the profiles of each halo in the appropriate bins
    for i in range(0, len(mean)):
        # re-initialize selection and for each bin
        selection = np.array([])
        for keyword in keywordList:
            rFilter = np.logical_and(np.greater_equal(bins[keyword][:-1], stackBins[i]), np.less(bins[keyword][:-1], stackBins[i+1]))
            selection = np.append(selection, stackingDictData[keyword][rFilter])
        if len(selection) == 0:
            continue
        else:
            mean[i] = np.nanmean(selection)
            stdDev[i] = np.nanstd(selection)
    
    return stackBins, mean, stdDev

def readHeader(snapPath):
    """
    INPUT:  snapPath:   path to the hdf5 file containing the simulation snapshot; example: h5py.File('/net/galaxy-data/export/galaxydata/EAGLE_snapshots/RecalL0025N0752/snapshot_012_z003p017/snap_012_z003p017.0.hdf5'

    RETURN: a:      the scale factor of the snapshot
            h:      hubble-factor of the snapshot
            boxsize: side length of the simulation cube
    """

    f = h5py.File(snapPath, 'r')

    a       = f['Header'].attrs.get('Time')         # Scale factor.
    h       = f['Header'].attrs.get('HubbleParam')  # h.
    boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].

    f.close()

    return a, h, boxsize

def stackProfiles_av(keywordList, stackingDictData, stackingDictBins, keywordToGnDict, binSize, cellLimit):
    """
    generates stacked profile (average) of the profiles stored in dictionary passed to the function
    INPUT:  keywordList: list holding the identifiers of the various profiles to be stacked (strings)
            stackingDictData: dictionary holding the profiles to be stacked, keywords are in the keywordList (strings)
            stackingDictBins: dictionary holding the bin edges of the profiles, same keywords as stackingDictData (strings)
            keywordToGnDict: dictionary linking groupnumber to the keywords (strings)
            binSize: binSize of the stacked profile in units of cells (float)
            cellLimit: max value to where the profile goes out to (float)
            len(stackingDictData[keyword]) = len(stackingDictBins[keyword]) - 1 
    RETURN: stackBins: the bin edges corresponding to the stacked profile (units of cells), 1D float numpy array
            average: the stacked, average profile; 1D float numpy array
            stdDev: standard deviation of the stacked profile; 1D float numpy array
            len(average) = len(stdDev) = len(stackBins) - 1 
    """
    # initialze the necessary variables
    bins = {}

    # cycle through the keywordList/List
    for keyword in keywordList:
        #gn = keywordToGnDict[keyword]
        bins[keyword] = stackingDictBins[keyword] 

    # define the bins for the stack
    stackBins = np.arange(0,cellLimit + 1, binSize, dtype=float)
    # initialize the median, perc25 and perc75 arrays
    average = np.zeros_like(stackBins[:-1], dtype=float)
    stdDev = np.zeros_like(stackBins[:-1], dtype=float)

    # cycle through the rescaledBins and stack the profiles of each halo in the appropriate bins
    for i in range(0, len(average)):
        # re-initialize selection and for each bin
        selection = np.array([])
        for keyword in keywordList:
            rFilter = np.logical_and(np.greater_equal(bins[keyword][:-1], stackBins[i]), np.less(bins[keyword][:-1], stackBins[i+1]))
            selection = np.append(selection, stackingDictData[keyword][rFilter])
        if len(selection) == 0:
            continue
        else:
            average[i] = np.average(selection)
            stdDev[i] = np.std(selection)
    
    return stackBins, average, stdDev

"""
def loadParticles(itype, snap_path, centre, load_region_length, attributes, bulkVel=np.array([0.,0.,0.])):
    ''' 
    modified version of the read_eagle function read_galaxy, modified by Marius
    loads the particles from the EAGLE snapshot in the region specified by centre and load_region_length
    INPUT:  itpye: integer indicating the particle type  
            snap_path: string containing the path to the snapshot
            centre: the coordinates of the galaxy's centre, directly from the SQL database (comoving Mpc)
            load_region_length : length of box in the snapshot's units (comoving Mpc/h)
            attributes: list holding the particle attributes to be loaded (see EAGLE paper for their names)

    RETURN: data:  the attributes of the itype particles, that have been loaded from the snapshot
    '''
    data = {}

    a, h, boxsize = readHeader(snap_path)

    #Intitialize read_eagle module
    eagle_data = EagleSnapshot(snap_path)

    #steph: centre comes from SQL and is stored in comoving coordinates, need to multiply by h so that the coordinates match those of the snapshot
    #steph: put centre into cMpc*h units
    #steph: centre_h (=centre*h) and load_region_length need to have the same units!
    centre_h = centre*h 

    #select a region to load, a 'load_ragion_length' cMpc*h cube centred on 'centre'
    region = np.array([(centre_h[0]-0.5*load_region_length), (centre_h[0]+0.5*load_region_length), (centre_h[1]-0.5*load_region_length), (centre_h[1]+0.5*load_region_length), (centre_h[2]-0.5*load_region_length), (centre_h[2]+0.5*load_region_length)])

    eagle_data.select_region(*region)

    #load data using read_eagle, load conversion factors manually
    f = h5py.File(snap_path, 'r')
    for att in attributes:
        tmp = eagle_data.read_dataset(itype, att)
        cgs = f['PartType%i/%s'%(itype, att)].attrs.get('CGSConversionFactor')
        aexp = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
        hexp = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')
        #this multiplication switches from comoving to proper
        data[att] = np.multiply(tmp, cgs * a**aexp * h**hexp, dtype='f8')
    
    if ('Coordinates' in attributes):
        #need to centre the coordinates on the galaxy's centre, i.e. subtract centre from coordinates
        #convert the centre into physical cgs units to do subtraction (same units as 'Coordinates')
        centre_phys = centre_h * h**-1 * a * f['PartType%i/Coordinates'%(itype)].attrs.get('CGSConversionFactor') 
        #correct for periodic box
        loadBox = load_region_length * h**-1 * a * f['PartType%i/Coordinates'%(itype)].attrs.get('CGSConversionFactor')
        boxSize_phys = boxsize * h**-1 *a * f['PartType%i/Coordinates'%(itype)].attrs.get('CGSConversionFactor')
        for i in range(3):
            if centre_phys[i] <= 0.5*loadBox:
                data['Coordinates'][:,i] = np.where(data['Coordinates'][:,i] > 0.5*(boxSize_phys), data['Coordinates'][:,i] - boxSize_phys, data['Coordinates'][:,i])
            elif centre_phys[i] >= (boxSize_phys - 0.5*loadBox):
                data['Coordinates'][:,i] = np.where(data['Coordinates'][:,i] < 0.5*(boxSize_phys), data['Coordinates'][:,i] + boxSize_phys, data['Coordinates'][:,i])


        data['cenCoordinates'] = data['Coordinates']-centre_phys
   
    if ('Velocity' in attributes):
        #need to subtract the bulk velocity from the velocity of the individual particles to centre on the halo
        #the units of bulkvel are expected to be km/s, need to be converted from km/s to cm/s
        bulkVel *= 1e5
        data['cenVelocity'] = data['Velocity'] - bulkVel

    f.close()

    return data
"""
