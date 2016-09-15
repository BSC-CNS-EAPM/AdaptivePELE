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
        contactThresholdDistance = 8
        contacts = pdb.countContacts(self.resname, contactThresholdDistance)

        threshold = self.thresholdCalculator(contacts)
        cluster = Cluster (pdb, thresholdRadius = threshold, contacts=contacts, metric=metric)
        self.clusters.addCluster(cluster)
        return len(self.clusters.clusters)-1

    def thresholdCalculator(self, contacts):

        numberOfContactsThresholdCompletelyBuried = 15
        numberOfContactsThresholdSemiBuried = 10
        numberOfContactsThresholdBulk = 1

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

        #return 2
        if contacts > numberOfContactsThresholdCompletelyBuried:
            return 2
        elif contacts > numberOfContactsThresholdSemiBuried:
            return 3
        elif contacts > numberOfContactsThresholdBulk:
            return 4
        else:
            return 4


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
        pdb_list = []
        metrics = []
        contactmaps = []
        preferences = []
        ids = []
        for trajectory in trajectories:
            trajNum = getTrajNum(trajectory)
            snapshots = getSnapshots(trajectory, True)
            if self.reportBaseFilename:
                reportFilename = os.path.join(os.path.split(trajectory)[0], self.reportBaseFilename%trajNum)
                metricstraj = np.loadtxt(reportFilename, usecols=(self.col,))
                if metricstraj.shape == ():
                    metricstraj = np.array([metricstraj])
            else:
                metricstraj = np.zeros(len(snapshots))
            #Prepare data for clustering
            metrics.extend(metricstraj.tolist())
            contactThresholdDistance = 8
            for snapshot in snapshots:
                pdb = atomset.PDB()
                pdb.initialise(snapshot)
                pdb_list.append(pdb)
                contactMap = pdb.createContactMap(self.resname,                    contactThresholdDistance)
                contactmaps.append(contactMap)
                preferences.append(contactMap.sum())

        preferences = float((min(preferences)-max(preferences))/2)
        cluster_center_indices, indices = clusterContactMaps(np.array(contactmaps),
                           preferences=preferences)
        center_ind = 0
        new_clusters = Clusters()
        new_ids = []
        new_contactmaps = []
        new_metrics = []
        elements = []
        for index in cluster_center_indices:
            cluster_members, = np.where(indices == center_ind)
            elements_in_cluster = cluster_members.size
            best_metric_ind = cluster_members[np.array(metrics)[cluster_members].argmin()]
            best_metric = metrics[best_metric_ind]
            cluster = Cluster(pdb_list[index],
                                    contactMap=contactmaps[index], metric=best_metric)
            new_clusters.addCluster(cluster)
            cluster.elements += elements_in_cluster-1
            new_ids.append('new:%d'%center_ind)
            new_contactmaps.append(cluster.contactMap)
            new_metrics.append(cluster.metric)
            elements.append(cluster.elements)
            center_ind += 1
        if len(self.clusters.clusters) > 0:
            new_clusters_lim = center_ind-1
            for clusterNum,cluster in enumerate(self.clusters.clusters):
                new_contactmaps.append(cluster.contactMap)
                new_ids.append("old:%d"%clusterNum)
                new_metrics.append(cluster.metric)
                elements.append(cluster.elements)

            preferences = -np.array(elements)
            cluster_center_indices, indices = clusterContactMaps(np.array(new_contactmaps),
                            preferences=preferences)
            center_ind = 0
            final_clusters = Clusters()
#            if cluster_center_indices is None:
                # if there are no contacts between ligand and protein the
                # clustering returns NaNs, it should not happen in a normal
                # run
#                import pdb as debug
#                debug.set_trace()
            for index in cluster_center_indices:
                cluster_index = int(new_ids[index].split(":")[-1])
                cluster_members, = np.where(indices == center_ind)
                elements_in_cluster = np.array(elements)[cluster_members].sum()
                best_metric_ind = cluster_members[np.array(new_metrics)[cluster_members].argmin()]
                # snapshot identified as exemplar by the algortihm
                # Not update center pdb structure of old clusters
                best_metric = new_metrics[best_metric_ind]
                if index <= new_clusters_lim:
                    cluster = new_clusters.clusters[cluster_index]
                else:
                    cluster = self.clusters.clusters[cluster_index]
                cluster.metric = best_metric
                cluster.elements += (elements_in_cluster-cluster.elements)
                final_clusters.addCluster(cluster)
                center_ind += 1
            self.clusters = final_clusters
        else:
            self.clusters = new_clusters

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

