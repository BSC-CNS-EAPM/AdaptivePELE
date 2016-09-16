import sys
import atomset
import glob
import numpy as np
import os
import blockNames
import clusteringTypes
import thresholdcalculator
import utilities
import pickle
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans

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

    def writeOutput(self, outputPath, degeneracy, outputObject):
        """
            Writes all the clustering information in outputPath

            outputPath [In] Folder that will contain all the clustering information
            degeneracy [In] Degeneracy of each cluster. It must be in the same order as in the self.clusters list
            outputObject [In] Output name for the pickle object
        """
        utilities.cleanup(outputPath) 
        utilities.makeFolder(outputPath)

        summaryFilename = os.path.join(outputPath, "summary.txt")
        summaryFile = open(summaryFilename, 'w')
        summaryFile.write("#cluster size degeneracy threshold contacts metric\n")

        for i,cluster in enumerate(self.clusters.clusters):
            outputFilename = "cluster_%d.pdb"%i
            outputFilename = os.path.join(outputPath, outputFilename)
            cluster.writePDB(outputFilename)
            summaryFile.write("%d %d %d %.1f %d %.3f\n"%(i, cluster.elements,
                                                         degeneracy[i],
                                                         cluster.threshold,
                                                         cluster.contacts,
                                                         cluster.metric))
        summaryFile.close()

        with open(outputObject, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)


class ContactsClustering(Clustering):
    def __init__(self, thresholdCalculator, resname=None, reportBaseFilename=None, columnOfReportFile=None):
        Clustering.__init__(self, resname, reportBaseFilename, columnOfReportFile)
        self.thresholdCalculator = thresholdCalculator

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


#TODO: refactor
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


class ContactMapAgglomerativeClustering(Clustering):
    def __init__(self, nclusters, resname=None, reportBaseFilename=None, columnOfReportFile=None):
        Clustering.__init__(self, resname, reportBaseFilename, columnOfReportFile)
        self.nclusters = nclusters

    def cluster(self, paths):
        """Clusters the snapshots of the trajectories provided using 
         a hierarchical clustering algorithm and the contactMaps similarity.

        The snapshots are clustered together with the previous found clusters.
        If and old cluster is found connected to a new one (with the new one
        being the exemplar) it is ignored and only the new snapshots are
        counted as members of the cluster. The minimum metric is used as the
        metric of the cluster"""
        trajectories = getAllTrajectories(paths)
        pdb_list = []
        metrics = []
        contactmaps = []
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
            for num, snapshot in enumerate(snapshots):
                pdb = atomset.PDB()
                pdb.initialise(snapshot, atomname="CA")
                pdb_list.append(pdb)
                contactMap = pdb.createContactMap(self.resname,                    contactThresholdDistance)
                contactmaps.append(contactMap)
        clusters, labels = clusterAgglomerativeContactMaps(np.array(contactmaps), self.nclusters)
        new_clusters = Clusters()
        new_ids = []
        new_contactmaps = []
        new_metrics = []
        elements = []
        for index in clusters:
            cluster_members, = np.where(labels == index)
            elements_in_cluster = cluster_members.size
            cluster_index = clusterKmeans(np.array(contactmaps)[cluster_members],1)
            cluster_center = cluster_members[cluster_index]
            best_metric_ind = cluster_members[np.array(metrics)[cluster_members].argmin()]
            best_metric = metrics[best_metric_ind]
            #debug.set_trace()
            cluster = Cluster(pdb_list[cluster_center],
                                contactMap=contactmaps[cluster_center], metric=best_metric)
            new_clusters.addCluster(cluster)
            cluster.elements += elements_in_cluster-1
            new_ids.append('new:%d'%index)
            new_contactmaps.append(cluster.contactMap)
            new_metrics.append(cluster.metric)
            elements.append(cluster.elements)

        cluster_limit = clusters.size
        if len(self.clusters.clusters) > 0:
            #recluster with previous clusters
            for clusterNum,cluster in enumerate(self.clusters.clusters):
                new_contactmaps.append(cluster.contactMap)
                new_ids.append("cluster:%d"%clusterNum)
                new_metrics.append(cluster.metric)
                elements.append(cluster.elements)

            final_clusters = Clusters()
            clusters, labels = clusterAgglomerativeContactMaps(np.array(new_contactmaps), self.nclusters)
            for index in clusters:
                cluster_members, = np.where(labels == index)
                elements_in_cluster = np.array(elements)[cluster_members].sum()
                print cluster_members
                cluster_index = clusterKmeans(np.array(new_contactmaps)[cluster_members],1)
                cluster_center = cluster_members[cluster_index]
                best_metric_ind = cluster_members[np.array(new_metrics)[cluster_members].argmin()]
                best_metric = new_metrics[best_metric_ind]
                cluster_index = int(new_ids[cluster_center].split(":")[-1])
                if cluster_center < cluster_limit:
                    cluster = new_clusters.clusters[cluster_index]
                else:
                    cluster = self.clusters.clusters[cluster_index]
                cluster.metric = best_metric
                cluster.elements += (elements_in_cluster-cluster.elements)
                final_clusters.addCluster(cluster)
                #debug.set_trace()
            self.clusters = final_clusters
        else:
            self.clusters = new_clusters

class ClusteringBuilder:
    def buildClustering(self, clusteringBlock, reportBaseFilename=None, columnOfReportFile=None):
        try:
            resname = clusteringBlock[blockNames.ClusteringTypes.ligandResname].upper()
            clusteringType = clusteringBlock[blockNames.ClusteringTypes.type]
        except KeyError as err:
            err.message=err.message + ": Need to provide mandatory parameter in clustering block"
            raise KeyError(err.message)

        if clusteringType == blockNames.ClusteringTypes.contacts:
            thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
            thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
            return ContactsClustering(thresholdCalculator, resname, reportBaseFilename, columnOfReportFile)
        elif clusteringType == blockNames.ClusteringTypes.contactMap:
            return ContactMapClustering(resname, reportBaseFilename, columnOfReportFile)
        elif clusteringType == blockNames.ClusteringTypes.agglomerative:
            nclusters = clusteringBlock[blockNames.ClusteringTypes.nclusters]
            return  ContactMapAgglomerativeClustering(nclusters, resname, reportBaseFilename, columnOfReportFile)
        else:
            sys.exit("Unknown clustering method! Choices are: " +
                     str(clusteringTypes.CLUSTERING_TYPE_TO_STRING_DICTIONARY.values()))

#TODO: should it be a class method?
def clusterContactMaps(contactmaps, preferences=None):
    contactmaps = contactmaps.reshape((contactmaps.shape[0],-1))
    affinitypropagation = AffinityPropagation(damping=0.9,preference=preferences,verbose=False).fit(contactmaps)
    cluster_center_indices = affinitypropagation.cluster_centers_indices_
    labels = affinitypropagation.labels_
    return cluster_center_indices, np.array(labels)

def clusterAgglomerativeContactMaps(contactmaps, n_clusters):
    contactmaps = contactmaps.reshape((contactmaps.shape[0],-1))
    agglomerative = AgglomerativeClustering(n_clusters=n_clusters,
                                            linkage='complete').fit(contactmaps)
    labels = agglomerative.labels_
    clusters = np.unique(labels)
    return clusters, labels

def clusterKmeans(contactmaps, n_clusters):
    contactmaps = contactmaps.reshape((contactmaps.shape[0],-1))
    kmeans = KMeans(n_clusters=n_clusters).fit(contactmaps)
    center = kmeans.cluster_centers_[0]
    contactmaps -= center
    distances = contactmaps.sum(axis=1)
    cluster_center = abs(distances).argmin()
    return cluster_center

def getAllTrajectories(paths):
    files = []
    for path in paths:
        files += glob.glob(path)
    return files

def getSnapshots(trajectoryFile, verbose=False):
    return atomset.getPDBSnapshots(trajectoryFile, verbose)

def getTrajNum(trajFilename):
    return int(trajFilename.split("_")[-1][:-4])

