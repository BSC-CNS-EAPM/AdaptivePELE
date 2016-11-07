import time
import os
import pickle
from atomset import atomset, RMSDCalculator
from clustering import clustering, clusteringTypes, thresholdcalculator
from spawning import spawning, densitycalculator
from analysis import analyseClustering


class OrderedContactsClustering(clustering.Clustering):
    def __init__(self, thresholdCalculator, resname=None,
                 reportBaseFilename=None, columnOfReportFile=None,
                 contactThresholdDistance=8, symmetries=[]):
        clustering.Clustering.__init__(self, resname, reportBaseFilename,
                                       columnOfReportFile,
                                       contactThresholdDistance)
        self.type = clusteringTypes.CLUSTERING_TYPES.contacts
        self.thresholdCalculator = thresholdCalculator
        self.symmetries = symmetries
        self.initialPDB = None
        self.maxThreshold = self.thresholdCalculator.getMaxThreshold()
        self.distancesList = []
        self.RMSDCalculator = RMSDCalculator.RMSDCalculator(symmetries)

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

        initial_scd = atomset.computeSquaredCentroidDifference(self.initialPDB,
                                                               pdb)
        snapshotPosition, closerClusterInd, closerClusterRMSD = self.getCloserCluster(initial_scd, threshold2, pdb)
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

    def getCloserCluster(self, initial_scd, threshold2, pdb):
        minimumCloseCluster = 0
        nclusters = len(self.clusters.clusters)
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
            clusterRMSD = self.RMSDCalculator.computeRMSD(cluster.pdb, pdb)
            if (clusterRMSD < closerClusterRMSD):
                closerClusterRMSD = clusterRMSD
                closerClusterInd = minimumCloseCluster
            minimumCloseCluster += 1
        return snapshotPosition, closerClusterInd, closerClusterRMSD


def printNumberSnapshotsEpoch(paths_report, i):
    trajs = clustering.getAllTrajectories(paths_report)
    total_snapshots = 0
    for traj in trajs:
        for line in open(traj, "r"):
            total_snapshots += 1
        total_snapshots -= 1
    print "Total snapsthots for epoch %d: %d" % (i, total_snapshots)


def clusterSnapshotsEpoch(path, i, ntrajs, clusteringObject):
    startTimeOrd = time.time()
    ClOrd.cluster(path)
    endTimeOrd = time.time()
    print "Total time of clustering ordered contacts, epoch %d: %.6f" % (i, endTimeOrd-startTimeOrd)
    print "Number of clusters ordered contacts epoch %d: %d" % (i, len(ClOrd.clusters.clusters))
    degeneraciesOrd = spawningObject.calculate(ClOrd.clusters.clusters, ntrajs, {})
    ClOrd.writeOutput("clsummary", degeneraciesOrd, clusteringObject, False)
    os.remove("clsummary/summary.txt")


def createPathway(initial_cluster, final_cluster, ClOrd):
    pathway = [initial_cluster]
    rowind = final_cluster
    while (rowind > 0):
        pathway.insert(1, rowind)
        clusterInit = ClOrd.clusters.clusters[rowind]
        cluster_distance = ClOrd.distancesList[rowind]
        minimumCloseCluster = rowind-1
        while (minimumCloseCluster > 0):
            distance_initial = ClOrd.distancesList[minimumCloseCluster]
            if (abs(distance_initial-cluster_distance) < distanceThreshold):
                break
            minimumCloseCluster -= 1

        minimumRMSD = 1000
        closerCluster = rowind-1
        while (minimumCloseCluster > 0):
            cluster = ClOrd.clusters.clusters[minimumCloseCluster]
            distance = ClOrd.distancesList[minimumCloseCluster]
            if (abs(distance-cluster_distance) > distanceThreshold):
                break
            clusterRMSD = ClOrd.RMSDCalculator.computeRMSD(cluster.pdb,
                                                           clusterInit.pdb)
            if (clusterRMSD < minimumRMSD):
                minimumRMSD = clusterRMSD
                closerCluster = minimumCloseCluster
            minimumCloseCluster -= 1
        rowind = closerCluster
    return pathway


def writePathwayTrajectory(ClOrd, pathway, filename):
    pathwayFile = open(filename, "w")
    pathwayFile.write("REMARK 000 File created using PELE++\n")
    pathwayFile.write("REMARK 000 Pathway trajectory created using findPathway program\n")
    pathwayFile.write("REMARK 000 List of cluster belonging to the pathway %s\n" % ' '.join(map(str, pathway)))
    for i, step_cluster in enumerate(pathway):
        cluster = ClOrd.clusters.clusters[step_cluster]
        pathwayFile.write("MODEL %d\n" % (i+1))
        pdbStr = cluster.pdb.pdb
        pdbList = pdbStr.split("\n")
        for line in pdbList:
            line = line.strip()
            # Avoid writing previous REMARK block
            if line.startswith("REMARK ") or line.startswith("MODEL "):
                continue
            elif line:
                pathwayFile.write(line+"\n")
        pathwayFile.write("ENDMDL\n")
    pathwayFile.close()


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
# symmetries = [{"3225:C3:AEN": "3227:C5:AEN", "3224:C2:AEN": "3228:C6:AEN"},
#               {"3230:N1:AEN": "3231:N2:AEN"}]
symmetries = {}
resname = "DAJ"
ClOrd = OrderedContactsClustering(thresholdCalculator, resname=resname,
                                  reportBaseFilename="report",
                                  columnOfReportFile=4, symmetries=symmetries)

spawningObject = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
trajFolder = "/gpfs/scratch/bsc72/bsc72755/adaptiveSampling/simulations/4DAJ_4_64_epsilon"
clusteringObject = "ClOrd_4DAJ.pkl"
if os.path.exists(clusteringObject):
    with open(clusteringObject, "r") as f:
        ClOrd = pickle.load(f)
else:
    allFolders = os.listdir(trajFolder)
    Epochs = [epoch for epoch in allFolders if epoch.isdigit()]
    for i in range(len(Epochs)):
        path = [trajFolder+"/%d/traj*" % i]
        paths_report = [trajFolder+"/%d/report*" % i]
        printNumberSnapshotsEpoch(paths_report, i)
        clusterSnapshotsEpoch(path, i, ntrajs, clusteringObject)

# obtain a connectivity matrix from microstates
# use graph algorithm to establish a path
distanceThreshold = 10
initial_cluster = 0
final_cluster = ClOrd.getOptimalMetric()
nclusters = ClOrd.clusters.getNumberClusters()
pathway = createPathway(initial_cluster, final_cluster, ClOrd)

# write pathway into a single trajectory
filename = "pathway4DAJ.pdb"
writePathwayTrajectory(ClOrd, pathway, filename)
# create clustering object with only the pathway clusters
ClPath = OrderedContactsClustering(thresholdCalculator, resname=resname,
                                   reportBaseFilename="report",
                                   columnOfReportFile=4, symmetries=symmetries)
ClPath.clusters.clusters = map(lambda x: ClOrd.clusters.clusters[x], pathway)

# spawning along the trajectory
spawningParams = spawning.SpawningParams()
spawningParams.epsilon = 0.25
spawningObject = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
spawningPathway = spawning.EpsilonDegeneracyCalculator()
# TODO:Try different spawining params
densityCalculator = densityCalculatorBuilder.build({})
degeneracies = spawningObject.calculate(ClPath.clusters.clusters, ntrajs, {})

title = "Degeneracies along pathway"
title2 = "Metrics along pathway"
comCoord, metrics = analyseClustering.extractCOMMatrix(ClPath.clusters.clusters,
                                                       resname)[:2]
analyseClustering.plotClusters(comCoord, degeneracies, title)
analyseClustering.plotClusters(comCoord, metrics, title2)
import matplotlib.pyplot as plt
plt.show()
