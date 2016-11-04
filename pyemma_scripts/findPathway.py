# import MSMblocks
# import helper
# import runMarkovChainModel as markov
import matplotlib.pyplot as plt
import numpy as np
import time
import os
import pickle
from atomset import atomset
from clustering import clustering, clusteringTypes, thresholdcalculator
from spawning import spawning, densitycalculator


class OrderedContactsClustering(clustering.Clustering):
    def __init__(self, thresholdCalculator, resname=None,
                 reportBaseFilename=None, columnOfReportFile=None,
                 contactThresholdDistance=8, symmetries={}):
        clustering.Clustering.__init__(self, resname, reportBaseFilename,
                                       columnOfReportFile,
                                       contactThresholdDistance)
        self.type = clusteringTypes.CLUSTERING_TYPES.contacts
        self.thresholdCalculator = thresholdCalculator
        self.symmetries = symmetries
        self.initialPDB = None
        self.maxThreshold = self.thresholdCalculator.getMaxThreshold()
        self.distancesList = []

    def getOptimalMetric(self):
        optimalMetric = 100
        optimalMetricIndex = 0
        for i, cluster in enumerate(self.clusters.clusters):
            if cluster.getMetric() < optimalMetric:
                optimalMetric = cluster.getMetric()
                optimalMetricIndex = i
        return optimalMetricIndex

    def addSnapshotToCluster(self, snapshot, metrics=[], col=None):
        pdb = atomset.PDB()
        pdb.initialise(snapshot, resname=self.resname)
        if not self.initialPDB:
            self.initialPDB = pdb

        contacts = pdb.countContacts(self.resname,
                                     self.contactThresholdDistance)
        numberOfLigandAtoms = pdb.getNumberOfAtoms()
        contactsPerAtom = float(contacts)/numberOfLigandAtoms

        threshold = self.thresholdCalculator.calculate(contactsPerAtom)
        threshold2 = threshold**2
        threshold2 = 10

        minimumCloseCluster = 0
        nclusters = len(self.clusters.clusters)
        initial_scd = atomset.computeSquaredCentroidDifference(self.initialPDB,
                                                               pdb)
        snapshotPosition = None
        for distance_initial in self.distancesList:
            if (snapshotPosition is None and distance_initial > initial_scd):
                snapshotPosition = minimumCloseCluster
            if (abs(distance_initial-initial_scd) < threshold2):
                break
            minimumCloseCluster += 1
        closerClusterInd = None
        closerClusterRMSD = 100

        while (minimumCloseCluster < nclusters):
            cluster = self.clusters.clusters[minimumCloseCluster]
            distance = self.distancesList[minimumCloseCluster]
            if (snapshotPosition is None and distance > initial_scd):
                snapshotPosition = minimumCloseCluster
            if (abs(distance-initial_scd) > threshold2):
                break
            clusterRMSD = atomset.computeRMSD(cluster.pdb, pdb, self.symmetries)
            if (clusterRMSD < closerClusterRMSD):
                closerClusterRMSD = clusterRMSD
                closerClusterInd = minimumCloseCluster
            minimumCloseCluster += 1

        try:
            clusterToAssign = self.clusters.clusters[closerClusterInd]
            if closerClusterRMSD < clusterToAssign.threshold:
                clusterToAssign.addElement(metrics)
                return
        except TypeError:
            # When the first snapshot is read the index of the closest cluster
            # center is set to None
            pass

        # if made it here, the snapshot was not added into any cluster

        cluster = clustering.Cluster(pdb, thresholdRadius=threshold,
                                     contacts=contactsPerAtom, metrics=metrics,
                                     metricCol=col)
        if snapshotPosition is not None:
            self.clusters.insertCluster(snapshotPosition, cluster)
            self.distancesList.insert(snapshotPosition, initial_scd)
        else:
            self.clusters.addCluster(cluster)
            self.distancesList.append(initial_scd)


ntrajs = 63
thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
thresholdCalculator = thresholdCalculatorBuilder.build({
        "thresholdCalculator": {
            "type": "constant",
            "params": {
                "value": 2
            }
        }
})
densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
densityCalculator = densityCalculatorBuilder.build({
       "density": {
            "type": "heaviside",
            "params": {
                 "conditions": [1.5, 1.0, 0.0],
                 "values": [8.0, 2.37, 1.0, 0.5]
             }
        }
   }
)
symmetries = {"3225:C3:AEN": "3227:C5:AEN", "3224:C2:AEN": "3228:C6:AEN",
              "3230:N1:AEN": "3231:N2:AEN"}
ClOrd = OrderedContactsClustering(thresholdCalculator, resname="AEN",
                                  reportBaseFilename="report",
                                  columnOfReportFile=4, symmetries=symmetries)

spawningObject = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
if os.path.exists("ClOrd_3PTB.pkl"):
    with open("ClOrd_3PTB.pkl", "r") as f:
        ClOrd = pickle.load(f)
else:
    for i in range(24):
        path = ["/gpfs/scratch/bsc72/bsc72021/simulation/3ptb_4_64/automate_inv_4/%d/traj*" % i]
        paths_report = ["/gpfs/scratch/bsc72/bsc72021/simulation/3ptb_4_64/automate_inv_4/%d/report*" % i]
        trajs = clustering.getAllTrajectories(paths_report)
        total_snapshots = 0
        for traj in trajs:
            for line in open(traj, "r"):
                total_snapshots += 1
            total_snapshots -= 1
        print "Total snapsthots for epoch %d: %d" % (i, total_snapshots)
        # startTimeCont = time.time()
        # ClCont.cluster(path)
        # endTimeCont = time.time()
        # print "Total time of clustering contacts, epoch %d: %.6f" % (i, endTimeCont-startTimeCont)
        # print "Number of clusters contacts epoch %d: %d" % (i, len(ClCont.clusters.clusters))
        # degeneraciesCont = spawningObject.calculate(ClCont.clusters.clusters, ntrajs, {})
        # ClCont.writeOutput("clsummary",degeneraciesCont,"ClCont_3PTB.pkl", False)
        # os.rename("clsummary/summary.txt", "results/summary_ClCont_3PTB.txt")
        startTimeOrd = time.time()
        ClOrd.cluster(path)
        endTimeOrd = time.time()
        print "Total time of clustering ordered contacts, epoch %d: %.6f" % (i, endTimeOrd-startTimeOrd)
        print "Number of clusters ordered contacts epoch %d: %d" % (i, len(ClOrd.clusters.clusters))
        degeneraciesOrd = spawningObject.calculate(ClOrd.clusters.clusters, ntrajs, {})
        ClOrd.writeOutput("clsummary",degeneraciesOrd,"ClOrd_3PTB.pkl", False)
        os.rename("clsummary/summary.txt", "results/summary_ClOrd_3PTB.txt")

# obtain a connectivity matrix from microstates
distanceThreshold = 10
initial_cluster = 0
final_cluster = ClOrd.getOptimalMetric()
nclusters = ClOrd.clusters.getNumberClusters()
distMatrix = np.zeros((nclusters, nclusters)) + 100
for rowind, clusterInit in enumerate(ClOrd.clusters.clusters):
    cluster_distance = ClOrd.distancesList[rowind]
    if rowind == final_cluster:
        continue
    # minimumCloseCluster = rowind+1
    minimumCloseCluster = 0
    for distance_initial in ClOrd.distancesList:
        if (abs(distance_initial-cluster_distance) < distanceThreshold):
            break
        minimumCloseCluster += 1

    while (minimumCloseCluster < nclusters):
        cluster = ClOrd.clusters.clusters[minimumCloseCluster]
        distance = ClOrd.distancesList[minimumCloseCluster]
        if (abs(distance-cluster_distance) > distanceThreshold):
            break
        clusterRMSD = atomset.computeRMSD(cluster.pdb, clusterInit.pdb,
                                          ClOrd.symmetries)
        distMatrix[rowind, minimumCloseCluster] = clusterRMSD
        # distMatrix[minimumCloseCluster, rowind] = -clusterRMSD
        minimumCloseCluster += 1
np.save("distMatrix.npy", distMatrix)
ClOrd.clusters.clusters[0].writePDB("cluster_0.pdb")
ClOrd.clusters.clusters[1].writePDB("cluster_1.pdb")
ClOrd.clusters.clusters[2].writePDB("cluster_2.pdb")
# plt.imshow(distMatrix)
# plt.colorbar()
# plt.savefig("connectivityMatrix.png")
# use graph algorithm to establish a path
