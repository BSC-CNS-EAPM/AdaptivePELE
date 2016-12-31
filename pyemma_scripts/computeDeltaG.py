import numpy as np
import helper
import glob
import sys
import argparse
from scipy.ndimage import filters

def parseArgs():
    parser = argparse.ArgumentParser(description="Script that computes delta G")
    parser.add_argument('trajectory', type=str)
    args = parser.parse_args()
    return args.trajectory

def writePDB(pmf_xyzg, title="clusters.pdb"):
    templateLine = "HETATM%s  H1  CLT L 502    %s%s%s  0.75%s           H\n"

    content = ""
    for i, line in enumerate(pmf_xyzg):
        number = str(i).rjust(5)
        x = ("%.3f"%line[0]).rjust(8)
        y = ("%.3f"%line[1]).rjust(8)
        z = ("%.3f"%line[2]).rjust(8)
        g = ("%.3f"%line[3]).rjust(8)

        content += templateLine%(number, x, y, z, g)
    f = open(title, 'w')
    f.write(content)
    f.close()


clusteringObject = helper.loadMSM('clustering_object.pkl')
allClusters = clusteringObject.clustercenters

MSMObject = helper.loadMSM('MSM_object.pkl')
pi = MSMObject.stationary_distribution


r = allClusters[MSMObject.connected_sets[0]]


#filename = "output.txt"
#data = np.loadtxt(filename)

#r = data[:,0:3]
#pi = data[:,3]

d = 1.0

trajWildcard = parseArgs()
originalFilenames = glob.glob(trajWildcard)

originalCoordinates = []
for i, originalFilename in enumerate(originalFilenames):
    trajOriginalCoordinates = np.loadtxt(originalFilename, usecols=(1,2,3))
    originalCoordinates.append(trajOriginalCoordinates)

maxval = 3*[-np.inf]
minval = 3*[np.inf]
for coord in originalCoordinates:
    cmaxval = coord.max(axis=0)
    cminval = coord.min(axis=0)
    maxval = np.maximum(cmaxval, maxval)
    minval = np.minimum(cminval, minval)

print "Maximum bounds", maxval, "Minimum bounds", minval

#Rounded floor and ceiling in intervals of "d" (e.g., floor of 1.73 with d = 0.5, will be 1.5 instead of 1.0, in order to optimize box creation.
#An extra box is included in the ceiling, so that all the points are contained in the range given by arange
bins = np.array([np.arange(np.floor(minval[i]) + d*int((minval[i] - np.floor(minval[i]))/d),
                    np.ceil(maxval[i]) + d*(int((maxval[i] - np.ceil(maxval[i]))/d) + 1),
        d) for i in range(3)])

"""
bins = [np.arange(np.floor(minval[i]),
                    np.ceil(maxval[i]),
        d) for i in range(3)]
"""


numberOfClusters = r.shape[0]
histogram = np.array([])
histogramFreq = np.array([])
histograms = []

microstateVolume = np.zeros(numberOfClusters)

print "Number of clusters", numberOfClusters



for i in range(numberOfClusters):
    allCoords = []
    for j,(trajOriginalCoordinates, dtraj) in enumerate(zip(originalCoordinates, clusteringObject.dtrajs)):
        assert dtraj.shape[0] == trajOriginalCoordinates.shape[0]
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

    current_hist, edges = np.histogramdd(allCoords, bins=bins)
    histograms.append(current_hist)

    #filtered_hist = filters.gaussian_filter(current_hist, sigma=1)
    filtered_hist = current_hist/current_hist.sum()

    nonZeroIndices = np.argwhere(filtered_hist > 0)

    # microstateVolume[i] = len(nonZeroIndices) * d**3

    if histogram.size == 0:
        histogramFreq = pi[i]*filtered_hist
        histogram = np.copy(current_hist)
    else:
        histogramFreq += pi[i]*filtered_hist
        histogram += current_hist


for i in range(numberOfClusters):
    histogramCluster = histograms[i]
    histogramTotal = histogram[histogramCluster > 0]
    histogramCluster = histogramCluster[histogramCluster > 0]
    microstateVolume[i] = (histogramCluster/histogramTotal).sum() * d**3


np.savetxt("volumeOfClusters.dat", microstateVolume)

microstateVolume = np.loadtxt("volumeOfClusters.dat")

kb = 0.0019872041
T = 300
beta = 1/(kb*T)

gpmf = -kb*T*np.log(pi/microstateVolume)
gpmf -= gpmf.min()
#debug
#for i,j in enumerate(gpmf):
#    print i, j, pi[i], microstateVolume[i]

deltaW = -gpmf[gpmf != np.inf].max()
print "bound    Delta G     Delta W     Binding Volume:     Binding Volume contribution"

upperGpmfValues = np.arange(0,-deltaW,0.5)

for upperGpmfValue in upperGpmfValues:
    bindingVolume = 0
    for g, volume in zip(gpmf, microstateVolume):
        if g <= upperGpmfValue:
            bindingVolume += np.exp(-beta * g) * volume 
    deltaG = deltaW - kb*T*np.log(bindingVolume/1661)
    print "%.1f\t%.3f\t%.3f\t%.3f\t%.3f" % (upperGpmfValue, deltaG, deltaW, bindingVolume, -kb*T*np.log(bindingVolume/1661))


pmf_xyzg = np.hstack((r, np.expand_dims(gpmf,axis=1)))
np.savetxt("pmf_xyzg.dat", pmf_xyzg)

writePDB(pmf_xyzg)


sys.exit()
"""
"""

gpmf = -kb*T*np.log(histogram)
gpmf -= gpmf.min()

#inExplorationRange = np.argwhere(gpmf != np.inf)
deltaW = -gpmf[gpmf != np.inf].max()

for indices in np.argwhere(gpmf != np.inf):
    i = indices[0]
    j = indices[1]
    k = indices[2]
    #print bins[0][i], bins[1][j], bins[2][k], gpmf[i,j,k]

upperGpmfValues = np.arange(0,-deltaW,0.5)
bindingVolumes = []
deltaGs = []
print "Delta G     Delta W     Binding Volume:     Binding Volume contribution"
for upperGpmfValue in upperGpmfValues:
    bindingVolume = 0
    for i,g in np.ndenumerate(gpmf[gpmf <= upperGpmfValue]):
        bindingVolume += np.exp(-beta*g)
    bindingVolume *= d**3

    deltaG = deltaW - kb*T*np.log(bindingVolume/1661)

    deltaGs.append(deltaG)
    bindingVolumes.append(bindingVolume)

    #np.savetxt("pmf_xyzg.dat", np.hstack((r, np.expand_dims(gpmf,axis=1))))

    print "%.3f\t%.3f\t%.3f\t%.3f" % (deltaG, deltaW, bindingVolume, -kb*T*np.log(bindingVolume/1661))

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(upperGpmfValues, deltaGs)
plt.figure(2)
plt.plot(upperGpmfValues, bindingVolumes)
plt.show()
