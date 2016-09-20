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
import pdb as debug

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
    def getContacts(self):
        return self.contacts

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


class ContactMapClustering(Clustering):
    def cluster(self, paths):
        """Cluster the snapshots of the trajectories provided using the
        affinity propagation algorithm and the contactMaps similarity.

        A double clustering process is going on. First, all snapshots from an
        epoch are clustered together, and then these clusters are clustered
        together with the previous clusters, if any
        
        Paths [in] list with the path to the trajectories to cluster"""
        trajectories = getAllTrajectories(paths)
        contactThresholdDistance = 8

        pdb_list, metrics, contactmaps = processSnapshots(trajectories, self.reportBaseFilename,
                         self.col, contactThresholdDistance, self.resname)
        preferences = map(np.sum,contactmaps)

        preferences = float((min(preferences)-max(preferences))/2)
        self.firstClusteringParams = clusteringResultsParameters(pdb_list=pdb_list,
                                                            metrics=metrics,
                                                            contactmaps=contactmaps)
        self.secondClusteringParams = clusteringResultsParameters(contactmaps=[])
        contactmapsArray = np.array(self.firstClusteringParams.contactmaps)
        cluster_center_indices, indices = clusterContactMaps(contactmapsArray, preferences=preferences)

        self.processClusterResults("first", cluster_center_indices, indices)

        if len(self.clusters.clusters) > 0:
            for clusterNum,cluster in enumerate(self.clusters.clusters):
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.ids.append("old:%d"%clusterNum)
                self.secondClusteringParams.metrics.append(cluster.metric)
                self.secondClusteringParams.elements.append(cluster.elements)

            preferences = float((min(self.secondClusteringParams.elements)-max(self.secondClusteringParams.elements))/2)
#            preferences = -np.array(self.secondClusteringParams.elements)
            #prefereneces = None
            secondcontactmapsArray = np.array(self.secondClusteringParams.contactmaps)
            cluster_center_indices, indices = clusterContactMaps(secondcontactmapsArray, preferences=preferences)
            debug.set_trace()
            self.processClusterResults("second", cluster_center_indices, indices)
            self.clusters = self.secondClusteringParams.newClusters
        else:
            self.clusters = self.firstClusteringParams.newClusters

    def processClusterResults(self, process, cluster_center_indices, indices):
        """ Analyse the results obtained from the affinity propagation
        algorithm.
        
        Process [In] String to decide wether it is the results of the first or
        second clustering
        cluster_center_indices [In] Numpy array with the snapshots taht are
        cluster centers
        Indices [In] Numpy array with the cluster each snapshot belongs to"""
        center_ind = 0
        for index in cluster_center_indices:

            cluster_members, = np.where(indices == center_ind)

            if process == "first":
                elements_in_cluster = cluster_members.size
                metricsArray = np.array(self.firstClusteringParams.metrics)
                best_metric_ind = cluster_members[metricsArray[cluster_members].argmin()]
                best_metric = self.firstClusteringParams.metrics[best_metric_ind]
                cluster = Cluster(self.firstClusteringParams.pdb_list[best_metric_ind],
                                        contactMap=self.firstClusteringParams.contactmaps[best_metric_ind], metric=best_metric)
                cluster.elements += elements_in_cluster-1
                self.secondClusteringParams.ids.append('new:%d'%center_ind)
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.metrics.append(cluster.metric)
                self.secondClusteringParams.elements.append(cluster.elements)
                self.firstClusteringParams.newClusters.addCluster(cluster)
            else:
                cluster_index = int(self.secondClusteringParams.ids[index].split(":")[-1])
                elementsArray = np.array(self.secondClusteringParams.elements)
                metricsArray = np.array(self.secondClusteringParams.metrics)
                elements_in_cluster = elementsArray[cluster_members].sum()
                best_metric_ind = cluster_members[metricsArray[cluster_members].argmin()]
                best_metric = self.secondClusteringParams.metrics[best_metric_ind]
                if index < len(self.firstClusteringParams.newClusters.clusters):
                    cluster = self.firstClusteringParams.newClusters.clusters[cluster_index]
                else:
                    cluster = self.clusters.clusters[cluster_index]
                cluster.metric = best_metric
                cluster.elements += (elements_in_cluster-cluster.elements)
                self.secondClusteringParams.newClusters.addCluster(cluster)

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
        contactThresholdDistance = 8

        pdb_list, metrics, contactmaps = processSnapshots(trajectories, self.reportBaseFilename,
                         self.col, contactThresholdDistance, self.resname)

        clusters, labels = clusterAgglomerativeContactMaps(np.array(contactmaps), self.nclusters)
        new_clusters = Clusters()
        new_ids = []
        new_contactmaps = []
        new_metrics = []
        elements = []
        for index in clusters:
            cluster_members, = np.where(labels == index)
            elements_in_cluster = cluster_members.size
            #cluster_index = clusterKmeans(np.array(contactmaps)[cluster_members],1)
            # Instead of using Kmeans to obtain a "center" or representative
            # cluster, select one of the members randomly
            cluster_index = np.random.randint(0, elements_in_cluster)
            cluster_center = cluster_members[cluster_index]
            best_metric_ind = cluster_members[np.array(metrics)[cluster_members].argmin()]
            best_metric = metrics[best_metric_ind]
            #debug.set_trace()
            cluster = Cluster(pdb_list[cluster_center],
                                contactMap=contactmaps[cluster_center], metric=best_metric)
            #cluster = Cluster(pdb_list[best_metric_ind],
            #                    contactMap=contactmaps[best_metric_ind], metric=best_metric)
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
                #cluster_index = clusterKmeans(np.array(new_contactmaps)[cluster_members],1)
                # Instead of using Kmeans to obtain a "center" or representative
                # cluster, select one of the members randomly
                cluster_index = np.random.randint(0, cluster_members.size)
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

class clusteringResultsParameters:
    def __init__(self, pdb_list=[], metrics=[], contactmaps=[]):
            self.pdb_list = pdb_list
            self.metrics = metrics
            self.contactmaps = contactmaps
            self.elements = []
            self.ids = []
            self.newClusters = Clusters()

#TODO: should it be a class method?
def clusterContactMaps(contactmaps, preferences=None):
    contactmaps = contactmaps.reshape((contactmaps.shape[0],-1))
    affinitypropagation = AffinityPropagation(damping=0.9,preference=preferences,verbose=False).fit(contactmaps)
    cluster_center_indices = affinitypropagation.cluster_centers_indices_
    labels = affinitypropagation.labels_
    debug.set_trace()
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

def processSnapshots(trajectories, reportBaseFilename, col, contactThresholdDistance, resname):
    pdb_list = []
    metrics = []
    contactmaps = []
    for trajectory in trajectories:
        trajNum = getTrajNum(trajectory)
        snapshots = getSnapshots(trajectory, True)
        if reportBaseFilename:
            reportFilename = os.path.join(os.path.split(trajectory)[0], reportBaseFilename%trajNum)
            metricstraj = np.loadtxt(reportFilename, usecols=(col,))
            if metricstraj.shape == ():
                metricstraj = np.array([metricstraj])
        else:
            metricstraj = np.zeros(len(snapshots))
        #Prepare data for clustering
        metrics.extend(metricstraj.tolist())
        for num, snapshot in enumerate(snapshots):
            pdb = atomset.PDB()
            pdb.initialise(snapshot, atomname="CA")
            pdb_list.append(pdb)
            contactMap = pdb.createContactMap(resname, contactThresholdDistance)
            contactmaps.append(contactMap)
    return pdb_list, metrics, contactmaps
