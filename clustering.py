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



class ContactMapClustering(Clustering):
    def cluster(self, paths):
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
                metrics[clusterNum+new_snapshot_limit+1] = cluster.metric
                # preferences[new_snapshot_limit+1+clusterNum] = cluster.elements
            cluster_center_indices, degeneracies, indices = clusterContactMaps(np.array(contactmaps))
            center_ind = 0
            clusters_to_pop=[]
            for index in cluster_center_indices:
                cluster_index = int(ids[index].split(":")[-1])
                cluster_members, = np.where(indices == center_ind)
                best_metric_ind = metrics[cluster_members].argmin()
                if best_metric_ind > new_snapshot_limit:
                    cluster_best_ind = int(ids[best_metric_ind].split(":")[-1])
                    cluster_best = self.clusters.clusters[cluster_best_ind]
                    if best_metric_ind != index:
                        clusters_to_pop.append(cluster_best_ind)
                    best_metric = cluster_best.metric
                    best_pdb = cluster_best.pdb
                    best_contactMap = cluster_best.contactMap
                    best_elements = cluster_best.elements
                else:
                    # snapshot identified as exemplar by the algortihm
                    best_metric = metrics[best_metric_ind]
                    best_pdb = pdb_list[best_metric_ind]
                    best_contactMap = contactmaps[best_metric_ind]
                    best_elements = 0 
                if index > new_snapshot_limit:
                    cluster = self.clusters.clusters[cluster_index]
                    cluster.addElement(best_metric)
                    cluster.elements += best_elements+degeneracies[indices[index]]-2
                    cluster.pdb = best_pdb
                    cluster.contactMap = best_contactMap
                    #One of the counts will account for the own cluster and the other
                    # for the addElement method
                else:
                    cluster = Cluster(best_pdb, contactMap=best_contactMap, metric=best_metric)
                    cluster.elements += best_elements+degeneracies[indices[index]]-1
                    self.clusters.addCluster(cluster)
                center_ind += 1
            #delete clusters that have been merged
            for pop_index in sorted(clusters_to_pop, reverse=True):
                del self.clusters.clusters[pop_index]

class ClusteringBuilder:
    def buildClustering(self, method, resname=None, reportBaseFilename=None, columnOfReportFile=None):
        if method == blockNames.ClusteringTypes.contacts:
            return ContactsClustering(resname, reportBaseFilename, columnOfReportFile)
        elif method == blockNames.ClusteringTypes.contactMap:
            return ContactMapClustering(resname, reportBaseFilename, columnOfReportFile)
        else:
            sys.exit("Unknown clustering method! Choices are: " +
                     str(blockNames.CLUSTERING_TYPE_TO_STRING_DICTIONARY.values()))

def clusterContactMaps(contactmaps, preferences=None):
    contactmaps = contactmaps.reshape((contactmaps.shape[0],-1))
    affinitypropagation = AffinityPropagation(damping=0.9,preference=preferences,verbose=False).fit(contactmaps)
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

