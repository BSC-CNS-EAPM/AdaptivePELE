import sys
import atomset
import glob
import numpy as np
import os
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

class contactsClustering(Clustering):
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
        threshold, contacts = self.thresholdCalculator(pdb, self.resname)
        cluster = Cluster (pdb, thresholdRadius = threshold, contacts=contacts, metric=metric)
        self.clusters.addCluster(cluster)
        return len(self.clusters.clusters)-1

    def thresholdCalculator(self, pdb, ligandResname):
        contactThresholdDistance = 8

        numberOfContactsThresholdCompletelyBuried = 15
        numberOfContactsThresholdSemiBuried = 10
        numberOfContactsThresholdBulk = 1

        contacts = pdb.countContacts(ligandResname, contactThresholdDistance)

        """
        if contacts > numberOfContactsThresholdCompletelyBuried:
            return 1, contacts
        elif contacts > numberOfContactsThresholdSemiBuried:
            return 3, contacts
        elif contacts > numberOfContactsThresholdBulk:
            return 5, contacts
        else:
            return 8, contacts

        #with PR, works diamonds
        if contacts > numberOfContactsThresholdCompletelyBuried:
            return 3, contacts
        elif contacts > numberOfContactsThresholdSemiBuried:
            return 5, contacts
        elif contacts > numberOfContactsThresholdBulk:
            return 7, contacts
        else:
            return 10, contacts

        if contacts > numberOfContactsThresholdCompletelyBuried:
            return 2, contacts
        elif contacts > numberOfContactsThresholdSemiBuried:
            return 4, contacts
        elif contacts > numberOfContactsThresholdBulk:
            return 6, contacts
        else:
            return 8, contacts
        """

        #return 2, contacts
        if contacts > numberOfContactsThresholdCompletelyBuried:
            return 2, contacts
        elif contacts > numberOfContactsThresholdSemiBuried:
            return 3, contacts
        elif contacts > numberOfContactsThresholdBulk:
            return 4, contacts
        else:
            return 4, contacts



class contactMapClustering(Clustering):
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
            #Prepare data for clustering
            contactmaps = []
            ids = []
            pdb_list = []
            preferences = np.zeros(len(snapshots)+len(self.clusters.clusters))
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
                preferences[new_snapshot_limit+clusterNum] = cluster.elements
            cluster_center_indices, degeneracies, indices = clusterContactMaps(np.array(contactmaps), preferences)
            for index in cluster_center_indices:
                cluster_index = int(ids[index].split(":")[-1])
                if index > new_snapshot_limit:
                    cluster_members, = np.where(indices == index)
                    cluster_members = cluster_members[cluster_members < new_snapshot_limit]
                    if metrics.size == 0:
                        best_metric = metrics[cluster_members].min()
                    else:
                        best_metric = 1e4
                    self.clusters.clusters[cluster_index].addElement(best_metric)
                    cluster.elements += degeneracies[indices[index]]-1
                else:
                    pdb = pdb_list[cluster_index]
                    contactMap = contactmaps[cluster_index]
                    metric = metrics[cluster_index]
                    cluster = Cluster (pdb, contactMap=contactMap, metric=metric)
                    cluster.elements += degeneracies[indices[index]]-1
                    self.clusters.addCluster(cluster)

class ClusteringBuilder:
    #TODO: add proper parameter handling for the builder(no hardcoded strings)
    def buildClustering(self, method, resname=None, reportBaseFilename=None, columnOfReportFile=None):
        if method == "contacts":
            return contactsClustering(resname, reportBaseFilename, columnOfReportFile)
        elif method == "contactmap":
            return contactMapClustering(resname, reportBaseFilename, columnOfReportFile)

def clusterContactMaps(contactmaps, preferences):
    contactmaps = contactmaps.reshape((contactmaps.shape[0],-1))
    affinitypropagation = AffinityPropagation(preference=preferences).fit(contactmaps)
    cluster_center_indices = affinitypropagation.cluster_centers_indices_
    labels = affinitypropagation.labels_
    labels,indices,degeneracies = np.unique(labels, return_counts=True, return_inverse=True)
    assert degeneracies.shape == cluster_center_indices.shape
    return cluster_center_indices, degeneracies, indices

def getAllTrajectories(paths):
    files = []
    for path in paths:
        files += glob.glob(path)
    return files

def getSnapshots(trajectoryFile, verbose=False):
    return atomset.getPDBSnapshots(trajectoryFile, verbose)

def getTrajNum(trajFilename):
    return int(trajFilename.split("_")[-1][:-4])

