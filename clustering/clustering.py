import sys
import atomset
import glob
import numpy as np
import os
import blockNames
from sklearn.cluster import AffinityPropagation

class Clusters:
    def __init__(self):
        self.clusters = []
    def addCluster(self, cluster):
        self.clusters.append(cluster)
    def printClusters(self, verbose=False):
        for i,cluster in enumerate(self.clusters):
            print "--------------"
            print "CLUSTER #%d" % i 
            print "--------------"
            cluster.printCluster(verbose)
            print ""
    def __eq__(self, other):
        return self.clusters == other.clusters

class Cluster:
    def __init__(self, pdb, thresholdRadius=0, contactMap=None, contacts=0, metric=0):
        self.pdb = pdb
        self.elements = 1
        self.threshold = thresholdRadius
        self.contacts = contacts
        self.contactMap = contactMap
        self.metric = metric
    def addElement(self, metric):
        self.elements += 1
        self.metric = min(metric, self.metric) #this for the moment
    def printCluster(self, verbose=False):
        if verbose:
            print self.pdb.printAtoms()
        print "Elements: ", self.elements
        print "Min energy: ", self.metric
        if self.threshold != 0:
            print "Radius threshold: ", self.threshold
        if not self.contactMap is None:
            print "Number of contacts in contact map", np.sum(self.contactMap)

    def writePDB(self, path):
        self.pdb.writePDB(path)
    def getMetric(self):
        return self.metric

    def __eq__(self, other):
        return self.pdb == other.pdb\
             and self.elements == other.elements\
             and self.threshold == other.threshold\
             and self.contacts == other.contacts\
             and abs(self.metric-other.metric) < 1e-7

from abc import ABCMeta, abstractmethod
class ThresholdCalculator():
    def __init__(self):
        self.type = "BaseClass" #change for abstract attribute
        pass

    @abstractmethod
    def calculate(self, contacts):
        pass

class ThresholdCalculatorConstant(ThresholdCalculator):
    def __init__(self, value):
        self.thresholdCalculator = ThresholdCalculator()
        self.type = "Constant"
        self.value

    def caclulate(self, contacts):
        return self.value

class ThresholdCalculatorHeaviside(ThresholdCalculator):
    def __init__(self, conditions=[15,10], values=[2,3,4]):
        self.thresholdCalculator = ThresholdCalculator()
        self.type = "Constant"
    
        if len(values) != len(conditions) and len(values) != len(conditions) + 1:
            raise ValueError('The number of values must be equal or one more, than the number of conditions')

        self.conditions = conditions
        self.values = values

    def calculate(self, contacts):
        for i in range(len(self.conditions)):
            #change, so that whole condition is in array
            if contacts > self.conditions[i]:
                return self.values[i]
        return self.values[-1]


class Clustering:
    def __init__(self, resname=None, reportBaseFilename=None, columnOfReportFile=None):
        self.clusters = Clusters()
        if reportBaseFilename:
            self.reportBaseFilename = reportBaseFilename + "_%d"
        else:
            self.reportBaseFilename = None
        self.resname = resname
        self.col = columnOfReportFile

    def __eq__(self, other):
        return self.clusters == other.clusters\
            and self.reportBaseFilename == other.reportBaseFilename\
            and self.resname == other.resname\
            and self.col == other.col


class ContactsClustering(Clustering):
    def __init__(self, resname=None, reportBaseFilename=None, columnOfReportFile=None):
        Clustering.__init__(self, resname, reportBaseFilename, columnOfReportFile)
        self.thresholdCalculator = ThresholdCalculatorHeaviside()

    def cluster(self, paths):
        trajectories = getAllTrajectories(paths)
        for trajectory in trajectories:
            trajNum = getTrajNum(trajectory)

            snapshots = getSnapshots(trajectory, True)
            if self.reportBaseFilename:
                reportFilename = os.path.join(os.path.split(trajectory)[0], self.reportBaseFilename%trajNum)
                metrics = np.loadtxt(reportFilename, usecols=(self.col,))
                if metrics.shape == ():
                    metrics = np.array([metrics])
            else:
                metrics = np.zeros(len(snapshots))

            for num, snapshot in enumerate(snapshots):
                self.addSnapshotToCluster(snapshot, metrics[num])

    def addSnapshotToCluster(self, snapshot, metric=0):
        pdb = atomset.PDB()
        pdb.initialise(snapshot, resname = self.resname)
        for clusterNum,cluster in enumerate(self.clusters.clusters):
            if atomset.computeRMSD(cluster.pdb, pdb) < cluster.threshold:
                cluster.addElement(metric)
                return

        #if made it here, the snapshot was not added into any cluster
        contactThresholdDistance = 8
        contacts = pdb.countContacts(self.resname, contactThresholdDistance)

        threshold = self.thresholdCalculator.calculate(contacts)
        cluster = Cluster (pdb, thresholdRadius = threshold, contacts=contacts, metric=metric)
        self.clusters.addCluster(cluster)
        return len(self.clusters.clusters)-1



class ContactMapClustering(Clustering):
    def cluster(self, paths):
        """Clusters the snapshots of the trajectories provided using the
        affinity propagation algorithm and the contactMaps similarity.

        The snapshots are clustered together with the previous found clusters.
        If and old cluster is found connected to a new one (with the new one
        being the exemplar) it is ignored and only the new snapshots are
        counted as members of the cluster. The minimum metric is used as the
        metric of the cluster"""
        trajectories = getAllTrajectories(paths)
        for trajectory in trajectories:
            trajNum = getTrajNum(trajectory)
            snapshots = getSnapshots(trajectory, True)
            metrics = np.zeros(len(snapshots)+len(self.clusters.clusters))
            if self.reportBaseFilename:
                reportFilename = os.path.join(os.path.split(trajectory)[0], self.reportBaseFilename%trajNum)
                metrics[:len(snapshots)] = np.loadtxt(reportFilename, usecols=(self.col,))
                if metrics.shape == ():
                    metrics = np.array([metrics])
            else:
                metrics = np.zeros(len(snapshots))
            #Prepare data for clustering
            contactmaps = []
            ids = []
            pdb_list = []
            # preferences = np.zeros(len(snapshots)+len(self.clusters.clusters))
            contactThresholdDistance = 8
            for num, snapshot in enumerate(snapshots):
                pdb = atomset.PDB()
                pdb.initialise(snapshot)
                pdb_list.append(pdb)
                contactmaps.append(pdb.createContactMap(self.resname,                    contactThresholdDistance))
                ids.append("snapshot:%d"%num)
            new_snapshot_limit = num
            for clusterNum,cluster in enumerate(self.clusters.clusters):
                contactmaps.append(cluster.contactMap)
                ids.append("cluster:%d"%clusterNum)
                #pdb_list.append(cluster.pdb)
                metrics[clusterNum+new_snapshot_limit+1] = cluster.metric
                # preferences[new_snapshot_limit+1+clusterNum] = cluster.elements
            cluster_center_indices, indices = clusterContactMaps(np.array(contactmaps))
            center_ind = 0
            for index in cluster_center_indices:
                cluster_index = int(ids[index].split(":")[-1])
                cluster_members, = np.where(indices[:new_snapshot_limit+1] == center_ind)
                elements_in_cluster = cluster_members.size
                if elements_in_cluster != 0:
                    best_metric_ind = cluster_members[metrics[cluster_members].argmin()]
                    # snapshot identified as exemplar by the algortihm
                    #best_pdb = pdb_list[best_metric_ind]
                    # Not update center pdb structure of old clusters
                    best_metric = metrics[best_metric_ind]
                    #best_contactMap = contactmaps[best_metric_ind]
                else:
                    #best_pdb = pdb_list[index]
                    best_metric = metrics[index]
                    best_contactMap = contactmaps[index]
                if index > new_snapshot_limit:
                    cluster = self.clusters.clusters[cluster_index]
                    #cluster.pdb = best_pdb
                    cluster.metric = best_metric
                    #cluster.contactMap = best_contactMap
                    cluster.elements += elements_in_cluster
                else:
                    cluster = Cluster(pdb_list[index],
                                      contactMap=contactmaps[index], metric=best_metric)
                    self.clusters.addCluster(cluster)
                    cluster.elements += elements_in_cluster-1
                center_ind += 1

class ClusteringBuilder:
    def buildClustering(self, method, resname=None, reportBaseFilename=None, columnOfReportFile=None):
        if method == blockNames.ClusteringTypes.contacts:
            return ContactsClustering(resname, reportBaseFilename, columnOfReportFile)
        elif method == blockNames.ClusteringTypes.contactMap:
            return ContactMapClustering(resname, reportBaseFilename, columnOfReportFile)
        else:
            sys.exit("Unknown clustering method! Choices are: " +
                     str(clusteringTypes.CLUSTERING_TYPE_TO_STRING_DICTIONARY.values()))

def clusterContactMaps(contactmaps, preferences=None):
    contactmaps = contactmaps.reshape((contactmaps.shape[0],-1))
    affinitypropagation = AffinityPropagation(damping=0.9,preference=preferences,verbose=False).fit(contactmaps)
    cluster_center_indices = affinitypropagation.cluster_centers_indices_
    labels = affinitypropagation.labels_
    return cluster_center_indices, np.array(labels)

def getAllTrajectories(paths):
    files = []
    for path in paths:
        files += glob.glob(path)
    return files

def getSnapshots(trajectoryFile, verbose=False):
    return atomset.getPDBSnapshots(trajectoryFile, verbose)

def getTrajNum(trajFilename):
    return int(trajFilename.split("_")[-1][:-4])

