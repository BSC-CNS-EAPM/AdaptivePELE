import sys
import atomset
import glob
import numpy as np
import os
import blockNames
import clusteringTypes
import thresholdcalculator
from utilities import utilities
import pickle
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
# TODO: to be removed when not used
import pdb as debug


class Clusters:
    def __init__(self):
        self.clusters = []

    def addCluster(self, cluster):
        self.clusters.append(cluster)

    def printClusters(self, verbose=False):
        for i, cluster in enumerate(self.clusters):
            print "--------------"
            print "CLUSTER #%d" % i
            print "--------------"
            cluster.printCluster(verbose)
            print ""

    def __eq__(self, other):
        return self.clusters == other.clusters


class Cluster:
    def __init__(self, pdb, thresholdRadius=0, contactMap=None, contacts=0,
                 metric=0, density=None):
        """
            contacts stands for contacts/ligandAtom
        """
        self.pdb = pdb
        self.elements = 1
        self.threshold = thresholdRadius
        self.density = density
        self.contacts = contacts
        self.contactMap = contactMap
        self.metric = metric

    def addElement(self, metric):
        self.elements += 1
        self.metric = min(metric, self.metric)  # this for the moment

    def printCluster(self, verbose=False):
        if verbose:
            print self.pdb.printAtoms()
        print "Elements: ", self.elements
        print "Min energy: ", self.metric
        if self.threshold != 0:
            print "Radius threshold: ", self.threshold
        if self.contactMap is not None:
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
    def __init__(self, resname=None, reportBaseFilename=None,
                 columnOfReportFile=None, contactThresholdDistance=8):
        self.clusters = Clusters()
        if reportBaseFilename:
            self.reportBaseFilename = reportBaseFilename + "_%d"
        else:
            self.reportBaseFilename = None
        self.resname = resname
        self.col = columnOfReportFile
        self.contactThresholdDistance = contactThresholdDistance

    def __eq__(self, other):
        return self.clusters == other.clusters\
            and self.reportBaseFilename == other.reportBaseFilename\
            and self.resname == other.resname\
            and self.col == other.col

    # Move the cluster methods of contactsClustering to the Clustering
    # superclass in order to make accessible to contactMapAccumaltiveClustering
    # and avoid duplicate code
    def cluster(self, paths):
        trajectories = getAllTrajectories(paths)
        for trajectory in trajectories:
            trajNum = utilities.getTrajNum(trajectory)

            snapshots = utilities.getSnapshots(trajectory, True)
            if self.reportBaseFilename:
                reportFilename = os.path.join(os.path.split(trajectory)[0],
                                              self.reportBaseFilename % trajNum)
                metrics = np.loadtxt(reportFilename, usecols=(self.col,))
                if metrics.shape == ():
                    metrics = np.array([metrics])
            else:
                metrics = np.zeros(len(snapshots))

            for num, snapshot in enumerate(snapshots):
                self.addSnapshotToCluster(snapshot, metrics[num])

    def writeOutput(self, outputPath, degeneracy, outputObject, writeAll):
        """
            Writes all the clustering information in outputPath

            outputPath [In] Folder that will contain all the clustering information
            degeneracy [In] Degeneracy of each cluster. It must be in the same order
            as in the self.clusters list
            outputObject [In] Output name for the pickle object
            writeAll [In] Boolean, wether to write pdb files for all cluster in addition
            of the summary
        """
        utilities.cleanup(outputPath)
        utilities.makeFolder(outputPath)

        summaryFilename = os.path.join(outputPath, "summary.txt")
        summaryFile = open(summaryFilename, 'w')
        summaryFile.write("#cluster size degeneracy contacts threshold density metric\n")

        for i, cluster in enumerate(self.clusters.clusters):
            if writeAll:
                outputFilename = "cluster_%d.pdb" % i
                outputFilename = os.path.join(outputPath, outputFilename)
                cluster.writePDB(outputFilename)

            if cluster.metric:
                writeString = "%d %d %d %.2f %.1f %.1f %.3f\n" % (i, cluster.elements,
                                                                degeneracy[i],
                                                                cluster.contacts,
                                                                cluster.threshold,
                                                                cluster.density,
                                                                cluster.metric)
            else:
                writeString = "%d %d %d %.2f %.1f %.1f -\n" % (i, cluster.elements,
                                                             degeneracy[i],
                                                             cluster.contacts,
                                                             cluster.threshold,
                                                             cluster.density)
            summaryFile.write(writeString)
        summaryFile.close()

        with open(outputObject, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)


class ContactsClustering(Clustering):
    def __init__(self, thresholdCalculator, resname=None,
                 reportBaseFilename=None, columnOfReportFile=None,
                 contactThresholdDistance=8):
        Clustering.__init__(self, resname, reportBaseFilename,
                            columnOfReportFile, contactThresholdDistance)
        self.thresholdCalculator = thresholdCalculator


    def addSnapshotToCluster(self, snapshot, metric=0):
        pdb = atomset.PDB()
        pdb.initialise(snapshot, resname=self.resname)
        for clusterNum, cluster in enumerate(self.clusters.clusters):
            if atomset.computeRMSD(cluster.pdb, pdb) < cluster.threshold:
                cluster.addElement(metric)
                return

        # if made it here, the snapshot was not added into any cluster
        contacts = pdb.countContacts(self.resname, self.contactThresholdDistance)
        numberOfLigandAtoms = pdb.getNumberOfAtoms()
        contactsPerAtom = float(contacts)/numberOfLigandAtoms

        threshold = self.thresholdCalculator.calculate(contactsPerAtom)
        cluster = Cluster(pdb, thresholdRadius=threshold, contacts=contactsPerAtom,
                          metric=metric)
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

        pdb_list, metrics, contactmaps = processSnapshots(trajectories, self.reportBaseFilename,
                         self.col, self.contactThresholdDistance, self.resname)
        preferences = map(np.sum,contactmaps)

        preferences = float((min(preferences)-max(preferences))/2)
        preferences = None
        self.firstClusteringParams = clusteringResultsParameters(pdb_list=pdb_list,
                                                            metrics=metrics,
                                                            contactmaps=contactmaps)
        self.secondClusteringParams = clusteringResultsParameters(contactmaps=[])
        contactmapsArray = np.array(self.firstClusteringParams.contactmaps)
        cluster_center_indices, indices = clusterContactMaps(contactmapsArray, preferences=preferences)

        self.processClusterResults("first", cluster_center_indices, indices)
        # print "Number of clusters after first clustering", len(self.firstClusteringParams.newClusters.clusters)
        if len(self.clusters.clusters) > 0:
            for clusterNum, cluster in enumerate(self.clusters.clusters):
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.ids.append("old:%d" % clusterNum)
                self.secondClusteringParams.metrics.append(cluster.metric)
                self.secondClusteringParams.elements.append(cluster.elements)

            preferences = float((min(self.secondClusteringParams.elements)-max(self.secondClusteringParams.elements))/2)
            preferences = 0.5
            # preferences = -np.array(self.secondClusteringParams.elements)
            # prefereneces = None
            secondcontactmapsArray = np.array(self.secondClusteringParams.contactmaps)
            cluster_center_indices, indices = clusterContactMaps(secondcontactmapsArray, preferences=preferences)
            self.processClusterResults("second", cluster_center_indices, indices)
            # print "Number of clusters after second clustering", len(self.secondClusteringParams.newClusters.clusters)
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

                cluster_index = selectRandomCenter(cluster_members,
                                                   metricsArray[cluster_members])

                contacts = self.firstClusteringParams.pdb_list[cluster_index].countContacts(self.resname,
                                                                 self.contactThresholdDistance)

                numberOfLigandAtoms = self.firstClusteringParams.pdb_list[cluster_index].getNumberOfAtoms()
                contactsPerAtom = float(contacts)/numberOfLigandAtoms

                cluster = Cluster(self.firstClusteringParams.pdb_list[cluster_index],
                                  contactMap=self.firstClusteringParams.contactmaps[cluster_index],
                                  contacts=contactsPerAtom, metric=best_metric)

                cluster.elements += elements_in_cluster-1
                self.secondClusteringParams.ids.append('new:%d' % center_ind)
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.metrics.append(cluster.metric)
                self.secondClusteringParams.elements.append(cluster.elements)
                self.firstClusteringParams.newClusters.addCluster(cluster)
            else:
                elementsArray = np.array(self.secondClusteringParams.elements)
                metricsArray = np.array(self.secondClusteringParams.metrics)
                elements_in_cluster = elementsArray[cluster_members].sum()
                index = selectRandomCenter(cluster_members,
                                           metricsArray[cluster_members])
                cluster_index = int(self.secondClusteringParams.ids[index].split(":")[-1])
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
    def __init__(self, nclusters, resname=None, reportBaseFilename=None,
                 columnOfReportFile=None, contactThresholdDistance=8):
        Clustering.__init__(self, resname, reportBaseFilename,
                            columnOfReportFile, contactThresholdDistance)
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

        pdb_list, metrics, contactmaps = processSnapshots(trajectories, self.reportBaseFilename,
                         self.col, self.contactThresholdDistance, self.resname)


        self.firstClusteringParams = clusteringResultsParameters(pdb_list=pdb_list,
                                     metrics=metrics, contactmaps=contactmaps)
        self.secondClusteringParams = clusteringResultsParameters(contactmaps=[])
        clusters, labels = clusterAgglomerativeContactMaps(np.array(contactmaps), self.nclusters)
        self.processClusterResults("first", clusters, labels)

        if len(self.clusters.clusters) > 0:
            for clusterNum, cluster in enumerate(self.clusters.clusters):
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.ids.append("old:%d" % clusterNum)
                self.secondClusteringParams.metrics.append(cluster.metric)
                self.secondClusteringParams.elements.append(cluster.elements)

            new_contactmaps = np.array(self.secondClusteringParams.contactmaps)
            clusters, labels = clusterAgglomerativeContactMaps(new_contactmaps,
                                                               self.nclusters)
            self.processClusterResults("second", clusters, labels)
            self.clusters = self.secondClusteringParams.newClusters
        else:
            self.clusters = self.firstClusteringParams.newClusters

    def processClusterResults(self, process, clusters, labels):
        """ Analyse the results obtained from the affinity propagation
        algorithm.

        Process [In] String to decide wether it is the results of the first or
        second clustering
        clusters [In] Numpy array with the indexes the clusters
        Labels [In] Numpy array with the cluster each snapshot belongs to"""
        for index in clusters:
            cluster_members, = np.where(labels == index)

            if process == "first":
                elements_in_cluster = cluster_members.size
                # contactmapsArray = np.array(self.firstClusteringParams.contactmaps)
                # cluster_index = clusterKmeans(contactmapsArray[cluster_members],1)
                # Instead of using Kmeans to obtain a "center" or representative
                # cluster, select one of the members randomly
                metricsArray = np.array(self.firstClusteringParams.metrics)
                best_metric_ind = cluster_members[metricsArray[cluster_members].argmin()]
                best_metric = self.firstClusteringParams.metrics[best_metric_ind]

                cluster_center = selectRandomCenter(cluster_members,
                                                    metricsArray[cluster_members])
                # cluster = Cluster(self.firstClusteringParams.pdb_list[best_metric_ind],
                #                   contactMap=self.firstClusteringParams.contactmaps[best_metric_ind], metric=best_metric)

                contacts = self.firstClusteringParams.pdb_list[cluster_center].countContacts(self.resname,
                                                                 self.contactThresholdDistance)

                numberOfLigandAtoms = self.firstClusteringParams.pdb_list[cluster_center].getNumberOfAtoms()
                contactsPerAtom = float(contacts)/numberOfLigandAtoms

                cluster = Cluster(self.firstClusteringParams.pdb_list[cluster_center],
                                  contactMap=self.firstClusteringParams.contactmaps[cluster_center],
                                  contacts=contactsPerAtom, metric=best_metric)

                cluster.elements += elements_in_cluster-1
                self.secondClusteringParams.ids.append('new:%d' % index)
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.metrics.append(cluster.metric)
                self.secondClusteringParams.elements.append(cluster.elements)
                self.firstClusteringParams.newClusters.addCluster(cluster)
            else:
                elementsArray = np.array(self.secondClusteringParams.elements)
                metricsArray = np.array(self.secondClusteringParams.metrics)
                elements_in_cluster = elementsArray[cluster_members].sum()
                best_metric_ind = cluster_members[metricsArray[cluster_members].argmin()]
                best_metric = self.secondClusteringParams.metrics[best_metric_ind]
                # contactmapsArray = np.array(self.firstClusteringParams.contactmaps)
                # cluster_index = clusterKmeans(contactmapsArray[cluster_members],1)
                # Instead of using Kmeans to obtain a "center" or representative
                # cluster, select one of the members randomly
                cluster_center = selectRandomCenter(cluster_members,
                                                    metricsArray[cluster_members])
                cluster_index = int(self.secondClusteringParams.ids[cluster_center].split(":")[-1])
                if cluster_center < self.nclusters:
                    cluster = self.firstClusteringParams.newClusters.clusters[cluster_index]
                else:
                    cluster = self.clusters.clusters[cluster_index]
                cluster.metric = best_metric
                cluster.elements += (elements_in_cluster-cluster.elements)
                self.secondClusteringParams.newClusters.addCluster(cluster)


class ContactMapAccumulativeClustering(Clustering):
    def __init__(self, thresholdCalculator, resname=None,
                 reportBaseFilename=None, columnOfReportFile=None,
                 contactThresholdDistance=8):
        Clustering.__init__(self, resname, reportBaseFilename,
                            columnOfReportFile, contactThresholdDistance)
        self.thresholdCalculator = thresholdCalculator

    def addSnapshotToCluster(self, snapshot, metric=0):
        pdb = atomset.PDB()
        pdb.initialise(snapshot, resname=self.resname)
        contactMap, contacts = pdb.createContactMap(self.resname, self.contactThresholdDistance)
        for clusterNum, cluster in enumerate(self.clusters.clusters):
            # correlation
            similarity = calculateCorrelationContactMaps(contactMap, cluster.contactMap)
            similarity += 1  # Necessary to omit negative correlations
            similarity /= 2.0  # Correlation values need to be higher now
            distance = 1-similarity

            # Jaccard index
            # intersectContactMaps = (contactMap == cluster.contactMap).sum()
            # averageContacts = contactMap.size + cluster.contactMap.size - intersectContactMaps
            # similarity = float(intersectContactMaps)/averageContacts
            # distance = 1-similarity

            # simple index
            # differenceContactMaps = np.abs(contactMap-cluster.contactMap).sum()
            # averageContacts = (0.5*(contactMap.sum()+cluster.contactMap.sum()))
            # distance = differenceContactMaps/averageContacts
            if distance == 0 or distance < cluster.threshold:
                cluster.addElement(metric)
                return

        # if made it here, the snapshot was not added into any cluster
        numberOfLigandAtoms = pdb.getNumberOfAtoms()
        contactsPerAtom = float(contacts)/numberOfLigandAtoms

        threshold = self.thresholdCalculator.calculate(contactsPerAtom)
        cluster = Cluster(pdb, thresholdRadius=threshold, metric=metric,
                          contacts=contactsPerAtom, contactMap=contactMap)
        self.clusters.addCluster(cluster)
        return len(self.clusters.clusters)-1


class ClusteringBuilder:
    def buildClustering(self, clusteringBlock, reportBaseFilename=None, columnOfReportFile=None):
        paramsBlock = clusteringBlock[blockNames.ClusteringTypes.params]
        try:
            resname = paramsBlock[blockNames.ClusteringTypes.ligandResname].upper()
            clusteringType = clusteringBlock[blockNames.ClusteringTypes.type]
            contactThresholdDistance = paramsBlock[blockNames.ClusteringTypes.contactThresholdDistance]
        except KeyError as err:
            err.message += ": Need to provide mandatory parameter in clustering block"
            raise KeyError(err.message)

        if clusteringType == blockNames.ClusteringTypes.contacts:
            thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
            thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
            return ContactsClustering(thresholdCalculator, resname,
                                      reportBaseFilename, columnOfReportFile)
        elif clusteringType == blockNames.ClusteringTypes.contactMapAffinity:
            return ContactMapClustering(resname, reportBaseFilename,
                                        columnOfReportFile,
                                        contactThresholdDistance)
        elif clusteringType == blockNames.ClusteringTypes.contactMapAgglomerative:
            nclusters = paramsBlock[blockNames.ClusteringTypes.nclusters]
            return ContactMapAgglomerativeClustering(nclusters, resname,
                                                     reportBaseFilename,
                                                     columnOfReportFile,
                                                     contactThresholdDistance)
        elif clusteringType == blockNames.ClusteringTypes.contactMapAccumulative:
            thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
            thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
            return ContactMapAccumulativeClustering(thresholdCalculator, resname,
                                      reportBaseFilename, columnOfReportFile)
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


# TODO: should it be a class method?
def clusterContactMaps(contactmaps, preferences=None):
    contactmaps = contactmaps.reshape((contactmaps.shape[0], -1))
    affinitypropagation = AffinityPropagation(damping=0.9,
                                              preference=preferences,
                                              verbose=False).fit(contactmaps)
    cluster_center_indices = affinitypropagation.cluster_centers_indices_
    labels = affinitypropagation.labels_
    return cluster_center_indices, np.array(labels)


def clusterAgglomerativeContactMaps(contactmaps, n_clusters):
    contactmaps = contactmaps.reshape((contactmaps.shape[0], -1))
    agglomerative = AgglomerativeClustering(n_clusters=n_clusters,
                                            linkage='complete').fit(contactmaps)
    labels = agglomerative.labels_
    clusters = np.unique(labels)
    return clusters, labels


def clusterKmeans(contactmaps, n_clusters):
    contactmaps = contactmaps.reshape((contactmaps.shape[0], -1))
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


def selectRandomCenter(cluster_members, metrics_weights):
    metrics_weights -= metrics_weights.max()
    if abs(metrics_weights.sum()) < 1e-8:
        metrics_weights = None
    else:
        T = 500
        kbT = 0.001987*T
        metrics_weights = np.exp(-metrics_weights/kbT)
        metrics_weights /= metrics_weights.sum()
    cluster_index = np.random.choice(cluster_members,
                                     p=metrics_weights)
    return cluster_index


def calculateCorrelationContactMaps(contactMap, clusterContactMap):
    contactMap1 = contactMap.reshape((1, -1))
    contactMap2 = clusterContactMap.reshape((1, -1))
    return np.corrcoef(contactMap1, contactMap2)[0, 1]


def processSnapshots(trajectories, reportBaseFilename, col,
                     contactThresholdDistance, resname):
    pdb_list = []
    metrics = []
    contactmaps = []
    for trajectory in trajectories:
        trajNum = utilities.getTrajNum(trajectory)
        snapshots = utilities.getSnapshots(trajectory, True)
        if reportBaseFilename:
            reportFilename = os.path.join(os.path.split(trajectory)[0],
                                          reportBaseFilename % trajNum)
            metricstraj = np.loadtxt(reportFilename, usecols=(col,))
            if metricstraj.shape == ():
                metricstraj = np.array([metricstraj])
        else:
            metricstraj = np.zeros(len(snapshots))
        # Prepare data for clustering
        metrics.extend(metricstraj.tolist())
        for num, snapshot in enumerate(snapshots):
            pdb = atomset.PDB()
            pdb.initialise(snapshot, resname=resname)
            pdb_list.append(pdb)
            contactMap, contacts = pdb.createContactMap(resname, contactThresholdDistance)
            contactmaps.append(contactMap)
    return pdb_list, metrics, contactmaps
