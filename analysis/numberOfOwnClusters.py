import socket
machine = socket.gethostname()
import matplotlib
if machine == "bsccv03":
   matplotlib.use('wxagg')
import matplotlib.pyplot as plt
import numpy as np
import os
import collections

clusteringFileSizeColumn = 5
#clusteringFileSizeColumn = 4
clusteringFolder = "clustering"
summaryFile = "summary.txt"
clusteringSummaryFile = os.path.join(clusteringFolder, summaryFile)
templetizedClusteringSummaryFile = os.path.join("%d", clusteringSummaryFile)


def getClusteringSummaryContent(file):
    if os.path.isfile(file): 
        return np.genfromtxt(file)
    else: return []


allFolders = os.listdir('.')
epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]
numberOfEpochs=int(len(epochFolders))


numberOfClusters = []
clustersPerEpoch = [{}]
for epoch in range(numberOfEpochs):
    clusteringSummary = getClusteringSummaryContent(templetizedClusteringSummaryFile%epoch)

    if clusteringSummary != []:
        numberOfClusters.append(len(clusteringSummary))
        epochDictionary = {}
        for line in clusteringSummary:
            size = line[clusteringFileSizeColumn]
            if not size in epochDictionary:
                epochDictionary[size] = len(np.argwhere(clusteringSummary[:,clusteringFileSizeColumn] == size))
        clustersPerEpoch.append(epochDictionary)

"""
for epochSummary in clustersPerEpoch:
    for size, numClusters in epochSummary.iteritems():
        if not size in clustersPerSize:
            clustersPerSize[size] = [numClusters]
        else:
            clustersPerSize[size].append(numClusters)
"""
clustersPerSize = collections.defaultdict(list)
for epochSummary in clustersPerEpoch:
    for size, numClusters in epochSummary.iteritems():
            clustersPerSize[size].append(numClusters)

for size, numClusters in clustersPerSize.iteritems():
    thisSizeClusters = clustersPerSize[size]
    withPaddedValues  = np.lib.pad(thisSizeClusters, (numberOfEpochs-1-len(thisSizeClusters),0), 'constant', constant_values=(0,)) #-1 because it assumes that in the last epoch, no clustering was done
    clustersPerSize[size] = withPaddedValues


sizes = clustersPerSize.keys()
sortedSizes = np.sort(sizes)
plt.plot(numberOfClusters, label="All clusters")
for size in sortedSizes:
    plt.plot(clustersPerSize[size], label=str(size))
plt.legend(loc=2)
#plt.title("n=64, same threshold, variable density")
#plt.savefig("../3ptb_4_64_numberOfClusters_6.png")
plt.show()
