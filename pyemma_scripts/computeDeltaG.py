import numpy as np
import helper
import glob
from scipy.ndimage import filters

clusteringObject = helper.loadMSM('clustering_object.pkl')
r = clusteringObject.clustercenters

MSMObject = helper.loadMSM('MSM_object.pkl')
pi = MSMObject.stationary_distribution

#filename = "output.txt"
#data = np.loadtxt(filename)

#r = data[:,0:3]
#pi = data[:,3]

#d = 0.05
d = 0.25

allPointsFilename = "test/MSM1/allPoints_1.dat"
allPoints = np.loadtxt(allPointsFilename)
maxval = allPoints.max(axis=0)
minval = allPoints.min(axis=0)

#Rounded floor and ceiling in intervals of "d" (e.g., floor of 1.73 with d = 0.5, will be 1.5 instead of 1.0, in order to optimize box creation.
#An extra box is included in the ceiling, so that all the points are contained in the range given by arange
bins = [np.arange(np.floor(minval[i]) + d*int((minval[i] - np.floor(minval[i]))/d),
                    np.ceil(maxval[i]) + d*(int((maxval[i] - np.ceil(maxval[i]))/d) + 1),
        d) for i in range(3)]

"""
bins = [np.arange(np.floor(minval[i]),
                    np.ceil(maxval[i]),
        d) for i in range(3)]
"""


numberOfClusters = r.shape[0]
originalFilenames = glob.glob("test/MSM1/traj_*")
histogram = np.array([])
microstateVolume = np.zeros(numberOfClusters)

print "Number of clusters", numberOfClusters

originalCoordinates = []
for originalFilename in originalFilenames:
    trajOriginalCoordinates = np.loadtxt(originalFilename, usecols=(1,2,3))
    originalCoordinates.append(trajOriginalCoordinates)


for i in range(numberOfClusters):
#for i in range(2):
    allCoords = []
    for trajOriginalCoordinates, dtraj in zip(originalCoordinates, clusteringObject.dtrajs):
        belongingFrames = np.argwhere(dtraj==i)
        trajCoords = trajOriginalCoordinates[belongingFrames, :]
        trajCoords = trajCoords.flatten().tolist()


        """
        if allCoords.size == 0:
            allCoords = np.copy(trajCoords)
        else:
            allCoords = np.vstack((allCoords, trajCoords))
        """
        allCoords.extend(trajCoords)

    allCoords = np.reshape(allCoords, (-1,3))

    current_hist, edges = np.histogramdd(allCoords, bins=bins, normed=True)

    filtered_hist = filters.gaussian_filter(current_hist, sigma=1)

    microstateVolume[i] = len(np.argwhere(filtered_hist > 0)) * d**3
    print i, microstateVolume[i]

    if histogram.size == 0:
        histogram = pi[i]*current_hist
    else:
        histogram += pi[i]*current_hist


np.savetxt("volumeOfClusters.dat", microstateVolume)

microstateVolume = np.loadtxt("volumeOfClusters.dat")

kb = 0.0019872041
T = 300
beta = 1/(kb*T)

"""
gpmf = -kb*T*np.log(pi/microstateVolume)
gpmf -= gpmf.min()

bindingVolume = 0
for g, volume in zip(gpmf, microstateVolume):
   bindingVolume += np.exp(-beta * g) * volume 

deltaW = -gpmf.max()

print r.shape, gpmf.shape
np.savetxt("pmf_xyzg.dat", np.hstack((r, np.expand_dims(gpmf,axis=1))))

"""
gpmf = -kb*T*np.log(histogram)
gpmf -= gpmf.min()

#inExplorationRange = np.argwhere(gpmf != np.inf)
deltaW = -gpmf[gpmf != np.inf].max()

bindingVolume = 0
for i,g in np.ndenumerate(gpmf):
    bindingVolume += np.exp(-beta*g)
bindingVolume *= d**3

deltaG = deltaW - kb*T*np.log(bindingVolume/1661)

#np.savetxt("pmf_xyzg.dat", np.hstack((r, np.expand_dims(gpmf,axis=1))))

print "Delta G     Delta W     Binding Volume:     Binding Volume contribution"
print "%.3f\t%.3f\t%.3f\t%.3f" % (deltaG, deltaW, bindingVolume, -kb*T*np.log(bindingVolume/1661))
