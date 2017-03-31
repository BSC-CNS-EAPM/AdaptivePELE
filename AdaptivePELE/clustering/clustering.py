import sys
import glob
import numpy as np
import os
import pickle
# import cPickle as pickle
import clusteringTypes
import thresholdcalculator
from AdaptivePELE.constants import blockNames
import AdaptivePELE.atomset.atomset as atomset
from AdaptivePELE.utilities import utilities
from AdaptivePELE.atomset import SymmetryContactMapEvaluator as sym
from AdaptivePELE.atomset import RMSDCalculator
from scipy import stats
import socket
import networkx as nx
import heapq
# if "bsccv" not in socket.gethostname():
#     from sklearn.cluster import AffinityPropagation
#     from sklearn.cluster import AgglomerativeClustering
#     from sklearn.cluster import KMeans


class Clusters:
    def __init__(self):
        self.clusters = []

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"clusters": self.clusters}
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.clusters = state['clusters']

    def addCluster(self, cluster):
        self.clusters.append(cluster)

    def insertCluster(self, index, cluster):
        """ Insert a cluster in a specified index"""
        self.clusters.insert(index, cluster)

    def getNumberClusters(self):
        return len(self.clusters)

    def getCluster(self, clusterNum):
        return self.clusters[clusterNum]

    def printClusters(self, verbose=False):
        for i, cluster in enumerate(self.clusters):
            print "--------------"
            print "CLUSTER #%d" % i
            print "--------------"
            cluster.printCluster(verbose)
            print ""

    def __eq__(self, other):
        return self.clusters == other.clusters


class AltStructures:
    """
        Helper class, each cluster will have an instance of AltStructures that
        will maintain a priority queue (pq) of alternative structures to spawn
        from encoded as tuples (priority, PDB).
    """
    def __init__(self):
        self.altStructPQ = []
        self.limitSize = 10

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"altStructPQ": self.altStructPQ, "limitSize": self.limitSize}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.limitSize = state['limitSize']
        self.altStructPQ = state['altStructPQ']

    def altSpawnSelection(self, centerPair):
        """
            Select an alternative PDB from the cluster center to spawn from
        """
        subpopulations = [i[0] for i in self.altStructPQ]
        totalSubpopulation = sum(subpopulations)
        # Create a list of the population distributed between the cluster
        # center and the alternative structures
        weights = 1.0/np.array([centerPair[0]-totalSubpopulation]+subpopulations)
        weights /= weights.sum()
        # This function only works on numpy >= 1.7, on life we have 1.6
        # ind = np.random.choice(range(len(self.altStructPQ)), p=weights)
        r = stats.rv_discrete(values=(range(self.sizePQ()), weights))
        ind = r.rvs(size=10)[0]
        # The first value of the distribution is always the cluster center
        if ind == 0:
            print "cluster center"
            return centerPair[1]
        else:
            # pick an alternative structure from the priority queue
            print "alternative structure"
            return self.altStructPQ[ind][1].pdb

    def cleanPQ(self):
        """
            Ensure that the alternative structures priority queue has no more
            elements than the limit in order to ensure efficiency
        """
        if len(self.altStructPQ) < self.limitSize:
            return
        limit = len(self.altStructPQ)
        del self.altStructPQ[self.limitSize-limit:]

    def addStructure(self, PDB, threshold, resname, contactThreshold, similarityEvaluator):
        i = 0
        for priority, subCluster in self.altStructPQ:
            isSimilar, distance = similarityEvaluator.isElement(PDB, subCluster, resname, contactThreshold)
            if distance < subCluster.threshold/2.0:
                subCluster.addElement([])
                del self.altStructPQ[i]
                heapq.heappush(self.altStructPQ, (subCluster.elements, subCluster))
                if len(self.altStructPQ) > 2*self.limitSize:
                    self.cleanPQ()
                return
            i += 1
        newCluster = Cluster(PDB, thresholdRadius=threshold, contactThreshold=contactThreshold, contactMap=similarityEvaluator.contactMap)
        heapq.heappush(self.altStructPQ, (1, newCluster))
        if len(self.altStructPQ) > 2*self.limitSize:
            self.cleanPQ()

    def sizePQ(self):
        return len(self.altStructPQ)


class Cluster:
    """
        A cluster contains a representative structure(pdb), the number of
        elements, its density, threshold, number of contacts,
        a contactMap(sometimes) and a metric
    """
    def __init__(self, pdb, thresholdRadius=None, contactMap=None,
                 contacts=None, metrics=[], metricCol=None, density=None,
                 contactThreshold=8, altSelection=False):
        """
            contacts stands for contacts/ligandAtom
        """
        self.pdb = pdb
        self.altStructure = AltStructures()
        self.elements = 1
        self.threshold = thresholdRadius
        self.density = density
        self.contacts = contacts
        self.contactMap = contactMap
        self.metrics = metrics
        self.metricCol = metricCol
        self.contactThreshold = contactThreshold
        self.altSelection = altSelection

        if self.threshold is None:
            self.threshold2 = None
        else:
            self.threshold2 = thresholdRadius*thresholdRadius

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"pdb": self.pdb, "altStructure": self.altStructure,
                 "elements": self.elements, "threshold": self.threshold,
                 "densitiy": self.density, "contacts": self.contacts,
                 "contactMap": self.contactMap,  "metrics": self.metrics,
                 "metricCol": self.metricCol, "threshold2": self.threshold2,
                 "contactThreshold": self.contactThreshold,
                 "altSelection": self.altSelection}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.pdb = state['pdb']
        self.altStructure = state.get('altStructPQ', AltStructures())
        self.elements = state['elements']
        self.threshold = state.get('threshold')
        self.density = state.get('density')
        self.contacts = state.get('contacts')
        self.contactMap = state.get('contactMap')
        self.metrics = state.get('metrics', [])
        self.metricCol = state.get('metricCol')
        self.threshold2 = state.get('threshold2')
        self.contactThreshold = state.get('contactThreshold', 8)
        self.altSelection = state.get('altSelection', False)

    def getMetric(self):
        if len(self.metrics):
            return self.metrics[self.metricCol]
        else:
            return None

    def getMetricFromColumn(self, numcol):
        if len(self.metrics):
            return self.metrics[numcol]
        else:
            return None

    def addElement(self, metrics):
        self.elements += 1
        if self.metrics is None:
            # Special case where cluster in created during clustering of
            # initial structures
            self.metrics = metrics
            return
        if len(metrics) and len(self.metrics):
            # Set all metrics to the minimum value
            self.metrics = np.minimum(self.metrics, metrics)

    def printCluster(self, verbose=False):
        if verbose:
            print self.pdb.printAtoms()
        print "Elements: ", self.elements
        print "Metrics: ", self.metrics
        if self.threshold != 0:
            print "Radius threshold: ", self.threshold
        print "Number of contacts: %.2f" % self.contacts

    def writePDB(self, path):
        self.pdb.writePDB(str(path))

    def getContacts(self):
        return self.contacts

    def writeSpawningStructure(self, path):
        """
            With 50 % probability select the cluster center to spawn in
            the next epoch
        """
        if not self.altSelection or self.altStructure.sizePQ() == 0:
            print "cluster center"
            self.pdb.writePDB(str(path))
        else:
            self.altStructure.altSpawnSelection((self.elements, self.pdb)).writePDB(str(path))

    def __eq__(self, other):
        return self.pdb == other.pdb\
             and self.elements == other.elements\
             and self.threshold == other.threshold\
             and self.contacts == other.contacts\
             and np.allclose(self.metrics, other.metrics)


class ContactsClusteringEvaluator:
    def __init__(self, RMSDCalculator):
        self.RMSDCalculator = RMSDCalculator
        self.contacts = None
        # Only here for compatibility purpose
        self.contactMap = None

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"RMSDCalculator": self.RMSDCalculator,
                 "contacts": self.contacts, "contactMap": self.contactMap}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.RMSDCalculator = state.get('RMSDCalculator', RMSDCalculator.RMSDCalculator())
        self.contacts = state.get('contacts')
        self.contactMap = state.get('contactMap')

    def isElement(self, pdb, cluster, resname, contactThresholdDistance):
        dist = self.RMSDCalculator.computeRMSD(cluster.pdb, pdb)
        return dist < cluster.threshold, dist

    def cleanContactMap(self):
        self.contactMap = None
        self.contacts = None

    def checkAttributes(self, pdb, resname, contactThresholdDistance):
        if self.contacts is None:
            self.contacts = pdb.countContacts(resname, contactThresholdDistance)

    def getOuterLimit(self, node):
        # return max(node.cluster.threshold2, node.cluster.threshold2 + 10 - node.depth)
        return node.cluster.threshold2
        if node.depth > 5:
            return node.cluster.threshold2
        else:
            return node.cluster.threshold2 + 2 ** (6-node.depth)

    def getInnerLimit(self, cluster):
        return cluster.threshold2


class CMClusteringEvaluator:
    limitSlope = {8: 6, 6: 15, 4: 60, 10: 3}
    limitMax = {8: 2, 6: 0.8, 4: 0.2, 10: 4}

    def __init__(self, similarityEvaluator, symmetryEvaluator):
        self.similarityEvaluator = similarityEvaluator
        self.symmetryEvaluator = symmetryEvaluator
        self.contacts = None
        self.contactMap = None

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"similarityEvaluator": self.similarityEvaluator,
                 "symmetryEvaluator": self.symmetryEvaluator,
                 "contacts": self.contacts, "contactMap": self.contactMap}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.similarityEvaluator = state.get('similarityEvaluator')
        self.symmetryEvaluator = state.get('symmetryEvaluator')
        self.contacts = state.get('contacts')
        self.contactMap = state.get('contactMap')

    def isElement(self, pdb, cluster, resname, contactThresholdDistance):
        if self.contactMap is None:
            self.contactMap, self.contacts = self.symmetryEvaluator.createContactMap(pdb, resname, contactThresholdDistance)
            # self.contactMap, foo = self.symmetryEvaluator.createContactMap(pdb, resname, contactThresholdDistance)
            # self.contacts = pdb.countContacts(resname, 8)  # contactThresholdDistance)
        distance = self.similarityEvaluator.isSimilarCluster(self.contactMap, cluster.contactMap, self.symmetryEvaluator)
        return distance < cluster.threshold, distance

    def cleanContactMap(self):
        self.contactMap = None
        self.contacts = None

    def checkAttributes(self, pdb, resname, contactThresholdDistance):
        if self.contactMap is None:
            self.contactMap, self.contacts = self.symmetryEvaluator.createContactMap(pdb, resname, contactThresholdDistance)
            # self.contactMap, foo = self.symmetryEvaluator.createContactMap(pdb, resname, contactThresholdDistance)
            # self.contacts = pdb.countContacts(resname, 8)  # contactThresholdDistance)

    def getInnerLimit(self, cluster):
        # if cluster.contacts > self.limitMax[cluster.contactThreshold]:
        #     return 4.0
        # else:
        #     return 16-self.limitSlope[cluster.contactThreshold]*cluster.contacts

        # if cluster.contacts > 2.0:
        #     return 4.0
        # elif cluster.contacts <= 0.5:
        #     return 25.0
        # else:
        #     return 25-14*(cluster.contacts-0.5)

        if cluster.contacts > 2.0:
            return 4.0
        elif cluster.contacts < 0.5:
            return 25.0
        else:
            return 16-8*(cluster.contacts-0.5)

        # if cluster.contacts > 1.0:
        #     return 4.0
        # elif cluster.contacts > 0.75:
        #     return 9.0
        # elif cluster.contacts > 0.5:
        #     return 16.0
        # else:
        #      return 25

    def getOuterLimit(self, node):
        # return max(16, 16 * 2 - node.depth)
        return 16


class Clustering:
    """
        Base class for clustering methods, it defines a cluster method that
        contacts and accumulative inherit and use

        resname [In] String containing the three letter name of the ligan in the pdb
        reportBaseFilename [In] Name of the file that contains the metrics of the snapshots to cluster
        columnOfReportFile [In] Column of the report file that contain the metric of interest
        contactThresholdDistance [In] Distance at wich a ligand atom and a protein atom are
        considered in contact(default 8)
    """
    def __init__(self, resname=None, reportBaseFilename=None,
                 columnOfReportFile=None, contactThresholdDistance=8,
                 altSelection=False):
        self.type = "BaseClass"

        self.clusters = Clusters()
        if reportBaseFilename:
            self.reportBaseFilename = reportBaseFilename + "_%d"
        else:
            self.reportBaseFilename = None
        self.resname = resname
        self.col = columnOfReportFile
        self.contactThresholdDistance = contactThresholdDistance
        self.symmetries = []
        self.altSelection = altSelection
        self.conformationNetwork = nx.DiGraph()
        self.epoch = -1

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"type": self.type, "clusters": self.clusters,
                 "reportBaseFilename": self.reportBaseFilename,
                 "resname": self.resname, "col": self.col, "epoch": self.epoch,
                 "symmetries": self.symmetries,
                 "conformationNetwork": self.conformationNetwork,
                 "contactThresholdDistance": self.contactThresholdDistance,
                 "altSelection": self.altSelection}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.type = state['type']
        self.clusters = state['clusters']
        self.reportBaseFilename = state.get('reportBaseFilename')
        self.resname = state.get('resname')
        self.col = state.get('col')
        self.contactThresholdDistance = state.get('contactThresholdDistance', 8)
        self.symmetries = state.get('symmetries', [])
        self.altSelection = state.get('altSelection', False)
        self.conformationNetwork = state.get('conformationNetwork', nx.DiGraph())
        self.epoch = state.get('metricCol', -1)

    def setCol(self, col):
        self.col = col

        for cluster in self.clusters.clusters:
            cluster.metricCol = col

    def getCluster(self, clusterNum):
        return self.clusters.getCluster(clusterNum)


    def getNumberClusters(self):
        return self.clusters.getNumberClusters()

    def __eq__(self, other):
        return self.clusters == other.clusters\
            and self.reportBaseFilename == other.reportBaseFilename\
            and self.resname == other.resname\
            and self.col == other.col

    def clusterInitialStructures(self, initialStructures):
        """
            Cluster the initial structures. This will allow to obtain a working
            processor to cluster mapping (see SimulationRunner docs for more info)
            for simulation with multiple initial structures
            :param initialStructures: List of the initial structures of the simulation
            :type initialStructures: list
        """
        clusterInitial = []
        for i, structurePath in enumerate(initialStructures):
            pdb = atomset.PDB()
            pdb.initialise(str(structurePath), resname=self.resname)
            for clusterNum, cluster in enumerate(self.clusters.clusters):
                scd = atomset.computeSquaredCentroidDifference(cluster.pdb, pdb)
                if scd > self.clusteringEvaluator.getInnerLimit(cluster):
                    continue

                isSimilar, dist = self.clusteringEvaluator.isElement(pdb, cluster,
                                                                    self.resname, self.contactThresholdDistance)
                if isSimilar:
                    cluster.addElement([])
                    clusterInitial.append(clusterNum)
                    self.clusters.clusters[clusterNum].elements = 0
                    break

            if len(clusterInitial) == i+1:
                # If an initial structure is similar enough to be added as
                # element in a previous cluster, move to the next structure,
                # this function is practically a duplicate of the
                # addSnapshotToCluster function and this blocks substitutes the
                # return statement after adding an element
                continue
            # if made it here, the snapshot was not added into any cluster
            # Check if contacts and contactMap are set (depending on which kind
            # of clustering)
            self.clusteringEvaluator.checkAttributes(pdb, self.resname, self.contactThresholdDistance)
            contacts = self.clusteringEvaluator.contacts
            numberOfLigandAtoms = pdb.getNumberOfAtoms()
            contactsPerAtom = float(contacts)/numberOfLigandAtoms

            threshold = self.thresholdCalculator.calculate(contactsPerAtom)
            cluster = Cluster(pdb, thresholdRadius=threshold,
                            contacts=contactsPerAtom,
                            contactMap=self.clusteringEvaluator.contactMap,
                            metrics=None, metricCol=self.col,
                            contactThreshold=self.contactThresholdDistance,
                            altSelection=self.altSelection)
            self.clusters.addCluster(cluster)
            clusterNum = self.getNumberClusters()-1
            clusterInitial.append(clusterNum)
            self.clusters.clusters[clusterNum].elements = 0
            self.conformationNetwork.add_node(clusterNum, parent='root', epoch=0)
        return clusterInitial

    def cluster(self, paths, processorsToClusterMapping):
        """
            Cluster the snaptshots contained in the pahts folder
            :param paths: List of folders with the snapshots
            :type paths: list
        """
        self.epoch += 1
        trajectories = getAllTrajectories(paths)
        for trajectory in trajectories:
            trajNum = utilities.getTrajNum(trajectory)
            origCluster = processorsToClusterMapping[trajNum-1]
            snapshots = utilities.getSnapshots(trajectory, True)
            if self.reportBaseFilename:
                reportFilename = os.path.join(os.path.split(trajectory)[0],
                                              self.reportBaseFilename % trajNum)
                metrics = np.loadtxt(reportFilename, ndmin=2)

                for num, snapshot in enumerate(snapshots):
                    origCluster = self.addSnapshotToCluster(snapshot, origCluster, metrics[num], self.col)
            else:
                for num, snapshot in enumerate(snapshots):
                    origCluster = self.addSnapshotToCluster(snapshot, origCluster)
        for cluster in self.clusters.clusters:
            cluster.altStructure.cleanPQ()

    def writeOutput(self, outputPath, degeneracy, outputObject, writeAll):
        """
            Writes all the clustering information in outputPath

            :param outputPath: Folder that will contain all the clustering information
            :type outputPath: str
            :param degeneracy: Degeneracy of each cluster. It must be in the same order
            as in the self.clusters list
            :type degeneracy: list
            :param outputObject: Output name for the pickle object
            :type outputObject: str
            :param writeAll: Wether to write pdb files for all cluster in addition
            of the summary
            :type writeAll: bool
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

            metric = cluster.getMetric()
            if metric is None:
                writeString = "%d %d %d %.2f %.4f %.1f %.3f\n" % (i, cluster.elements,
                                                                  degeneracy[i],
                                                                  cluster.contacts,
                                                                  cluster.threshold,
                                                                  cluster.density,
                                                                  metric)
            else:
                writeString = "%d %d %d %.2f %.4f %.1f -\n" % (i, cluster.elements,
                                                               degeneracy[i],
                                                               cluster.contacts,
                                                               cluster.threshold,
                                                               cluster.density)
            summaryFile.write(writeString)
        summaryFile.close()

        with open(outputObject, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

    def addSnapshotToCluster(self, snapshot, origCluster, metrics=[], col=None):
        pdb = atomset.PDB()
        pdb.initialise(snapshot, resname=self.resname)
        self.clusteringEvaluator.cleanContactMap()
        for clusterNum, cluster in enumerate(self.clusters.clusters):
            scd = atomset.computeSquaredCentroidDifference(cluster.pdb, pdb)
            if scd > self.clusteringEvaluator.getInnerLimit(cluster):
                continue

            isSimilar, dist = self.clusteringEvaluator.isElement(pdb, cluster,
                                                                 self.resname, self.contactThresholdDistance)
            if isSimilar:
                if dist > cluster.threshold/2:
                    cluster.altStructure.addStructure(pdb, cluster.threshold, self.resname, self.contactThresholdDistance, self.clusteringEvaluator)
                cluster.addElement(metrics)
                if self.conformationNetwork.has_edge(origCluster, clusterNum):
                    self.conformationNetwork[origCluster][clusterNum]['transition'] += 1
                else:
                    self.conformationNetwork.add_edge(origCluster, clusterNum, transition=1)
                return origCluster

        # if made it here, the snapshot was not added into any cluster
        # Check if contacts and contactMap are set (depending on which kind
        # of clustering)
        self.clusteringEvaluator.checkAttributes(pdb, self.resname, self.contactThresholdDistance)
        contacts = self.clusteringEvaluator.contacts
        numberOfLigandAtoms = pdb.getNumberOfAtoms()
        contactsPerAtom = float(contacts)/numberOfLigandAtoms

        threshold = self.thresholdCalculator.calculate(contactsPerAtom)
        cluster = Cluster(pdb, thresholdRadius=threshold,
                          contacts=contactsPerAtom,
                          contactMap=self.clusteringEvaluator.contactMap,
                          metrics=metrics, metricCol=col,
                          contactThreshold=self.contactThresholdDistance,
                          altSelection=self.altSelection)
        self.clusters.addCluster(cluster)
        clusterNum = self.clusters.getNumberClusters()-1
        if clusterNum == origCluster:
            # The clusterNum should only be equal to origCluster when the first
            # cluster is created and the clusterInitialStructures function has
            # not been called, i.e. when usind the compareClustering script
            self.conformationNetwork.add_node(clusterNum, parent='root', epoch=self.epoch)
        else:
            self.conformationNetwork.add_node(clusterNum, parent=origCluster, epoch=self.epoch)
        self.conformationNetwork.add_edge(origCluster, clusterNum, transition=1)
        # If a new cluster is discovered during a trajectory, the next step in
        # the same trajectory will be considered to start from these new
        # cluster, thus resulting in a more precise conformation network and
        # smoother pathways
        return clusterNum

    def writeConformationNetwork(self, path):
        """
            Write the conformational network to file to visualize it
            :param path: Path where to write the network
            :type path: str
        """
        nx.write_edgelist(self.conformationNetwork, path)

    def writeFDT(self, path):
        """
            Write the first discovery tree to file in edgelist format to
            visualize it
            :param path: Path where to write the network
            :type path: str
        """
        with open(path, "w") as fw:
            for node, data in self.conformationNetwork.nodes_iter(data=True):
                if data['parent'] != 'root':
                    fw.write("%d\t%d\n" % (data['parent'], node))

    def writeConformationNodeMetric(self, path, metricCol):
        """
            Write the metric of each node in the conformation network in a
            tab-separated file
            :param path: Path where to write the network
            :type path: str
            :param metricCol: Column of the metric of interest
            :type metricCol: int
        """
        with open(path, "w") as f:
            for i, cluster in enumerate(self.clusters.clusters):
                metric = cluster.getMetricFromColumn(metricCol)
                if metric is None:
                    f.write("%d\t-\n" % i)
                else:
                    f.write("%d\t%.4f\n" % (i, metric))

    def writeConformationNodePopulation(self, path):
        """
            Write the population of each node in the conformation network in a
            tab-separated file
            :param path: Path where to write the network
            :type path: str
        """
        with open(path, "w") as f:
            for i, cluster in enumerate(self.clusters.clusters):
                f.write("%d\t%d\n" % (i, cluster.elements))


    def createPathwayToCluster(self, clusterLeave):
        """
            Retrace the FDT from a specific cluster to the root where it was
            discovered
            :param clusterLeave: End point of the pathway to reconstruct
            :type clusterLeave: int
        """
        pathway = []
        nodeLabel = clusterLeave
        while nodeLabel != "root":
            pathway.append(nodeLabel)
            nodeLabel = self.conformationNetwork.node[nodeLabel]['parent']
        return pathway[::-1]

    def getOptimalMetric(self, column=None):
        """
            Find the cluster with the best metric
        """
        optimalMetric = 100
        optimalMetricIndex = 0
        for i, cluster in enumerate(self.clusters.clusters):
            if column is None:
                metric = cluster.getMetric()
            else:
                metric = cluster.getMetricFromColumn(column)
            if metric < optimalMetric:
                optimalMetric = metric
                optimalMetricIndex = i
        return optimalMetricIndex

    def writePathwayTrajectory(self, pathway, filename):
        """
            Write a list of cluster forming a pathway into a trajectory pdb file
            :param pathway: List of clusters that form the pathway
            :type pathway: list
            :param filename: Path where to write the trajectory
            :type filename: str
        """
        pathwayFile = open(filename, "w")
        pathwayFile.write("REMARK 000 File created using PELE++\n")
        pathwayFile.write("REMARK 000 Pathway trajectory created using the FDT\n")
        pathwayFile.write("REMARK 000 List of cluster belonging to the pathway %s\n" % ' '.join(map(str, pathway)))
        for i, step_cluster in enumerate(pathway):
            cluster = self.clusters.clusters[step_cluster]
            pathwayFile.write("MODEL %d\n" % (i+1))
            pdbStr = cluster.pdb.pdb
            pdbList = pdbStr.split("\n")
            for line in pdbList:
                line = line.strip()
                # Avoid writing previous REMARK block
                if line.startswith("REMARK ") or line.startswith("MODEL ") or line == "END":
                    continue
                elif line:
                    pathwayFile.write(line+"\n")
            pathwayFile.write("ENDMDL\n")
        pathwayFile.close()

    def writePathwayOptimalCluster(self, filename):
        """
            Extracte the pathway to the cluster with the best metric as a
            trajectory and  write it to a PDB file
            :param filename: Path where to write the trajectory
            :type filename: str
        """
        optimalCluster = self.getOptimalMetric()
        pathway = self.createPathwayToCluster(optimalCluster)
        self.writePathwayTrajectory(pathway, filename)

    def calculateMetastabilityIndex(self):
        """
            Calculate the metastablity index, mI. mI is the ratio of transitions
            from a given cluster that remains within the same cluster
        """
        metInd = {}
        for node in self.conformationNetwork.nodes_iter():
            totalOut = 0
            selfTrans = 0
            for n, edge, data in self.conformationNetwork.out_edges_iter(node, data=True):
                if edge == node:
                    selfTrans = data['transition']
                totalOut += data['transition']
            if totalOut+selfTrans:
                metInd[node] = selfTrans/float(totalOut)
            else:
                metInd[node] = 0
        return metInd

    def writeMetastabilityIndex(self, filename, metInd=None):
        """
            Write the metastability index of each node to file.
            :param filename: Path where to write the trajectory
            :type filename: str
        """
        if metInd is None:
            metInd = self.calculateMetastabilityIndex()
        with open(filename, "w") as f:
            for node, met in metInd.iteritems():
                f.write("%d\t%.4f\n" % (node, met))


class ContactsClustering(Clustering):
    """
        Cluster together all snapshots that are closer to the cluster center
        than certain threshold. This threshold is assigned according to the
        ratio of number of contacts over the number of heavy atoms of the ligand

        thresholdCalculator [In] ThresholdCalculator object that calculate the
        threshold according to the contacts ratio
        resname [In] String containing the three letter name of the ligand
        in the pdb
        reportBaseFilename [In] Name of the file that contains the metrics of
        the snapshots to cluster
        columnOfReportFile [In] Column of the report file that contain the
        metric of interest
        contactThresholdDistance [In] Distance at wich a ligand atom and a protein atom are
        considered in contact(default 8)
    """
    def __init__(self, thresholdCalculator, resname=None,
                 reportBaseFilename=None, columnOfReportFile=None,
                 contactThresholdDistance=8, symmetries=[], altSelection=False):
        Clustering.__init__(self, resname, reportBaseFilename,
                           columnOfReportFile, contactThresholdDistance,
                            altSelection=altSelection)
        self.type = clusteringTypes.CLUSTERING_TYPES.contacts
        self.thresholdCalculator = thresholdCalculator
        self.symmetries = symmetries
        self.clusteringEvaluator = ContactsClusteringEvaluator(RMSDCalculator.RMSDCalculator(symmetries))

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"type": self.type, "clusters": self.clusters,
                 "reportBaseFilename": self.reportBaseFilename,
                 "resname": self.resname, "col": self.col, "epoch": self.epoch,
                 "symmetries": self.symmetries,
                 "conformationNetwork": self.conformationNetwork,
                 "contactThresholdDistance": self.contactThresholdDistance,
                 "altSelection": self.altSelection,
                 "thresholdCalculator": self.thresholdCalculator,
                 "clusteringEvaluator": self.clusteringEvaluator}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.type = state['type']
        self.clusters = state['clusters']
        self.reportBaseFilename = state.get('reportBaseFilename')
        self.resname = state.get('resname')
        self.col = state.get('col')
        self.contactThresholdDistance = state.get('contactThresholdDistance', 8)
        self.symmetries = state.get('symmetries', [])
        self.altSelection = state.get('altSelection', False)
        self.conformationNetwork = state.get('conformationNetwork', nx.DiGraph())
        self.epoch = state.get('metricCol', -1)
        self.thresholdCalculator = state.get('thresholdCalculator', thresholdcalculator.ThresholdCalculatorConstant())
        self.clusteringEvaluator = state.get('clusteringEvaluator', ContactsClusteringEvaluator(RMSDCalculator.RMSDCalculator(self.symmetries)))


class ContactMapAccumulativeClustering(Clustering):
    """ Cluster together all snapshots that have similar enough contactMaps.
        This similarity can be calculated with different methods (see similariyEvaluator documentation)

        thresholdCalculator [In] ThresholdCalculator object that calculate the threshold
        similarityEvaluator [In] SimilarityEvaluator object that will determine wether two snapshots
        are similar enough to belong to the same cluster
        resname [In] String containing the three letter name of the ligan in the pdb
        reportBaseFilename [In] Name of the file that contains the metrics of the snapshots to cluster
        columnOfReportFile [In] Column of the report file that contain the metric of interest
        contactThresholdDistance [In] Distance at wich a ligand atom and a protein
        atom are considered in contact(default 8)
    """
    def __init__(self, thresholdCalculator, similarityEvaluator, resname=None,
                 reportBaseFilename=None, columnOfReportFile=None,
                 contactThresholdDistance=8, symmetries=[], altSelection=False):
        Clustering.__init__(self, resname, reportBaseFilename,
                            columnOfReportFile, contactThresholdDistance,
                            altSelection=altSelection)
        self.type = clusteringTypes.CLUSTERING_TYPES.contactMapAccumulative
        self.thresholdCalculator = thresholdCalculator
        self.similarityEvaluator = similarityEvaluator
        self.symmetryEvaluator = sym.SymmetryContactMapEvaluator(symmetries)
        self.clusteringEvaluator = CMClusteringEvaluator(similarityEvaluator, self.symmetryEvaluator)

    def __getstate__(self):
        # Defining pickling interface to avoid problems when working with old
        # simulations if the properties of the clustering-related classes have
        # changed
        state = {"type": self.type, "clusters": self.clusters,
                 "reportBaseFilename": self.reportBaseFilename,
                 "resname": self.resname, "col": self.col, "epoch": self.epoch,
                 "symmetries": self.symmetries,
                 "conformationNetwork": self.conformationNetwork,
                 "contactThresholdDistance": self.contactThresholdDistance,
                 "altSelection": self.altSelection,
                 "thresholdCalculator": self.thresholdCalculator,
                 "similariyEvaluator": self.similarityEvaluator,
                 "symmetryEvaluator": self.symmetryEvaluator,
                 "clusteringEvaluator": self.clusteringEvaluator}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.type = state['type']
        self.clusters = state['clusters']
        self.reportBaseFilename = state.get('reportBaseFilename')
        self.resname = state.get('resname')
        self.col = state.get('col')
        self.contactThresholdDistance = state.get('contactThresholdDistance', 8)
        self.symmetries = state.get('symmetries', [])
        self.altSelection = state.get('altSelection', False)
        self.conformationNetwork = state.get('conformationNetwork', nx.DiGraph())
        self.epoch = state.get('metricCol', -1)
        self.thresholdCalculator = state.get('thresholdCalculator', thresholdcalculator.ThresholdCalculatorConstant(value=0.3))
        self.similarityEvaluator = state.get('similariyEvaluator', JaccardEvaluator())
        self.symmetryEvaluator = state.get('symmetryEvaluator', sym.SymmetryContactMapEvaluator(self.symmetries))
        self.clusteringEvaluator = state.get('clusteringEvaluator', CMClusteringEvaluator(self.similarityEvaluator, self.symmetryEvaluator))


class SequentialLastSnapshotClustering(Clustering):
    """
        Assigned  the last snapshot of the trajectory to a cluster.
        Only useful for PELE sequential runs
    """
    def cluster(self, paths):
        """
            Cluster the snaptshots contained in the pahts folder
            paths [In] List of folders with the snapshots
        """
        # Clean clusters at every step, so we only have the last snapshot of
        # each trajectory as clusters
        self.clusters = Clusters()
        trajectories = getAllTrajectories(paths)
        for trajectory in trajectories:
            trajNum = utilities.getTrajNum(trajectory)

            snapshots = utilities.getSnapshots(trajectory, True)
            if self.reportBaseFilename:
                reportFilename = os.path.join(os.path.split(trajectory)[0],
                                              self.reportBaseFilename % trajNum)
                metrics = np.loadtxt(reportFilename, ndmin=2)
                # Pass as cluster metrics the minimum value for each metric,
                # thus the metrics are not valid to do any spawning, only to
                # check the exit condition
                metrics = metrics.min(axis=0)

                self.addSnapshotToCluster(snapshots[-1], metrics, self.col)
            else:
                self.addSnapshotToCluster(snapshots[-1])

    def addSnapshotToCluster(self, snapshot, metrics=[], col=None):
        pdb = atomset.PDB()
        pdb.initialise(snapshot, resname=self.resname)
        contacts = pdb.countContacts(self.resname,
                                     self.contactThresholdDistance)
        numberOfLigandAtoms = pdb.getNumberOfAtoms()
        contactsPerAtom = float(contacts)/numberOfLigandAtoms

        cluster = Cluster(pdb, thresholdRadius=0,
                          contacts=contactsPerAtom, metrics=metrics,
                          metricCol=col)
        self.clusters.addCluster(cluster)


class ContactMapClustering(Clustering):
    def __init__(self, resname=None, reportBaseFilename=None,
                 columnOfReportFile=None, contactThresholdDistance=8,
                 symmetries=[]):
        Clustering.__init__(self, resname, reportBaseFilename,
                            columnOfReportFile, contactThresholdDistance)
        self.type = clusteringTypes.CLUSTERING_TYPES.contactMapAffinity
        self.symmetryEvaluator = sym.SymmetryContactMapEvaluator(symmetries)

    def cluster(self, paths):
        """Cluster the snapshots of the trajectories provided using the
        affinity propagation algorithm and the contactMaps similarity.

        A double clustering process is going on. First, all snapshots from an
        epoch are clustered together, and then these clusters are clustered
        together with the previous clusters, if any

        Paths [in] list with the path to the trajectories to cluster"""
        trajectories = getAllTrajectories(paths)

        pdb_list, metrics, contactmaps, contacts = processSnapshots(trajectories, self.reportBaseFilename,
                         self.col, self.contactThresholdDistance, self.resname, self.symmetryEvaluator)
        preferences = map(np.sum,contactmaps)

        preferences = float((min(preferences)-max(preferences))/2)
        preferences = None
        self.firstClusteringParams = clusteringResultsParameters(pdb_list=pdb_list,
                                                            metrics=metrics,
                                                            contactmaps=contactmaps,contacts=contacts)
        self.secondClusteringParams = clusteringResultsParameters(contactmaps=[])
        contactmapsArray = np.array(self.firstClusteringParams.contactmaps)
        cluster_center_indices, indices = clusterContactMaps(contactmapsArray, preferences=preferences)

        self.processClusterResults("first", cluster_center_indices, indices)
        # print "Number of clusters after first clustering", len(self.firstClusteringParams.newClusters.clusters)
        if len(self.clusters.clusters) > 0:
            for clusterNum, cluster in enumerate(self.clusters.clusters):
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.ids.append("old:%d" % clusterNum)
                self.secondClusteringParams.metrics.append(cluster.metrics)
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
                best_metric_ind = cluster_members[metricsArray[cluster_members, self.col].argmin()]
                best_metrics = self.firstClusteringParams.metrics[best_metric_ind]

                cluster_index = selectRandomCenter(cluster_members,
                                                   metricsArray[cluster_members, self.col])

                contacts = self.firstClusteringParams.contacts[cluster_index]
                numberOfLigandAtoms = self.firstClusteringParams.pdb_list[cluster_index].getNumberOfAtoms()
                contactsPerAtom = float(contacts)/numberOfLigandAtoms

                cluster = Cluster(self.firstClusteringParams.pdb_list[cluster_index],
                                  contactMap=self.firstClusteringParams.contactmaps[cluster_index],
                                  contacts=contactsPerAtom, metrics=best_metrics, metricCol=self.col)

                cluster.elements += elements_in_cluster-1
                self.secondClusteringParams.ids.append('new:%d' % center_ind)
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.metrics.append(cluster.metrics)
                self.secondClusteringParams.elements.append(cluster.elements)
                self.firstClusteringParams.newClusters.addCluster(cluster)
            else:
                elementsArray = np.array(self.secondClusteringParams.elements)
                metricsArray = np.array(self.secondClusteringParams.metrics)
                elements_in_cluster = elementsArray[cluster_members].sum()
                index = selectRandomCenter(cluster_members,
                                           metricsArray[cluster_members, self.col])
                cluster_index = int(self.secondClusteringParams.ids[index].split(":")[-1])
                best_metric_ind = cluster_members[metricsArray[cluster_members, self.col].argmin()]
                best_metrics = self.secondClusteringParams.metrics[best_metric_ind]
                if index < len(self.firstClusteringParams.newClusters.clusters):
                    cluster = self.firstClusteringParams.newClusters.clusters[cluster_index]
                else:
                    cluster = self.clusters.clusters[cluster_index]
                cluster.metrics = best_metrics
                cluster.elements += (elements_in_cluster-cluster.elements)
                self.secondClusteringParams.newClusters.addCluster(cluster)

            center_ind += 1


class ContactMapAgglomerativeClustering(Clustering):
    """nclusters[In] Number of the cluster to generate with the hierarchical clustering algorithm
        resname [In] String containing the three letter name of the ligan in the pdb
        reportBaseFilename [In] Name of the file that contains the metrics of the snapshots to cluster
        columnOfReportFile [In] Column of the report file that contain the metric of interest
        contactThresholdDistance [In] Distance at wich a ligand atom and a protein atom are
        considered in contact(default 8)
    """
    def __init__(self, nclusters, resname=None, reportBaseFilename=None,
                 columnOfReportFile=None, contactThresholdDistance=8,
                 symmetries=[]):
        Clustering.__init__(self, resname, reportBaseFilename,
                            columnOfReportFile, contactThresholdDistance)
        self.type = clusteringTypes.CLUSTERING_TYPES.contactMapAgglomerative
        self.nclusters = nclusters
        self.symmetryEvaluator = sym.SymmetryContactMapEvaluator(symmetries)

    def cluster(self, paths):
        """Clusters the snapshots of the trajectories provided using
         a hierarchical clustering algorithm and the contactMaps similarity.

        The snapshots are clustered together with the previous found clusters.
        If and old cluster is found connected to a new one (with the new one
        being the exemplar) it is ignored and only the new snapshots are
        counted as members of the cluster. The minimum metric is used as the
        metric of the cluster"""
        trajectories = getAllTrajectories(paths)

        pdb_list, metrics, contactmaps, contacts = processSnapshots(trajectories, self.reportBaseFilename,
                         self.col, self.contactThresholdDistance, self.resname, self.symmetryEvaluator)


        self.firstClusteringParams = clusteringResultsParameters(pdb_list=pdb_list,
                                     metrics=metrics, contactmaps=contactmaps, contacts=contacts)
        self.secondClusteringParams = clusteringResultsParameters(contactmaps=[])
        clusters, labels = clusterAgglomerativeContactMaps(np.array(contactmaps), self.nclusters)
        self.processClusterResults("first", clusters, labels)

        if len(self.clusters.clusters) > 0:
            for clusterNum, cluster in enumerate(self.clusters.clusters):
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.ids.append("old:%d" % clusterNum)
                self.secondClusteringParams.metrics.append(cluster.metrics)
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
                best_metric_ind = cluster_members[metricsArray[cluster_members, self.col].argmin()]
                best_metrics = self.firstClusteringParams.metrics[best_metric_ind]

                cluster_center = selectRandomCenter(cluster_members,
                                                    metricsArray[cluster_members, self.col])
                # cluster = Cluster(self.firstClusteringParams.pdb_list[best_metric_ind],
                #                   contactMap=self.firstClusteringParams.contactmaps[best_metric_ind], metric=best_metric)

                contacts = self.firstClusteringParams.contacts[cluster_center]
                numberOfLigandAtoms = self.firstClusteringParams.pdb_list[cluster_center].getNumberOfAtoms()
                contactsPerAtom = float(contacts)/numberOfLigandAtoms

                cluster = Cluster(self.firstClusteringParams.pdb_list[cluster_center],
                                  contactMap=self.firstClusteringParams.contactmaps[cluster_center],
                                  contacts=contactsPerAtom, metrics=best_metrics, metricCol=self.col)

                cluster.elements += elements_in_cluster-1
                self.secondClusteringParams.ids.append('new:%d' % index)
                self.secondClusteringParams.contactmaps.append(cluster.contactMap)
                self.secondClusteringParams.metrics.append(cluster.metrics)
                self.secondClusteringParams.elements.append(cluster.elements)
                self.firstClusteringParams.newClusters.addCluster(cluster)
            else:
                elementsArray = np.array(self.secondClusteringParams.elements)
                metricsArray = np.array(self.secondClusteringParams.metrics)
                elements_in_cluster = elementsArray[cluster_members].sum()
                best_metric_ind = cluster_members[metricsArray[cluster_members, self.col].argmin()]
                best_metrics = self.secondClusteringParams.metrics[best_metric_ind]
                # contactmapsArray = np.array(self.firstClusteringParams.contactmaps)
                # cluster_index = clusterKmeans(contactmapsArray[cluster_members],1)
                # Instead of using Kmeans to obtain a "center" or representative
                # cluster, select one of the members randomly
                cluster_center = selectRandomCenter(cluster_members,
                                                    metricsArray[cluster_members, self.col])
                cluster_index = int(self.secondClusteringParams.ids[cluster_center].split(":")[-1])
                if cluster_center < self.nclusters:
                    cluster = self.firstClusteringParams.newClusters.clusters[cluster_index]
                else:
                    cluster = self.clusters.clusters[cluster_index]
                cluster.metrics = best_metrics
                cluster.elements += (elements_in_cluster-cluster.elements)
                self.secondClusteringParams.newClusters.addCluster(cluster)


class ClusteringBuilder:
    def buildClustering(self, clusteringBlock, reportBaseFilename=None, columnOfReportFile=None):
        paramsBlock = clusteringBlock[blockNames.ClusteringTypes.params]
        try:
            resname = str(paramsBlock[blockNames.ClusteringTypes.ligandResname].upper())
            clusteringType = clusteringBlock[blockNames.ClusteringTypes.type]
            contactThresholdDistance = paramsBlock[blockNames.ClusteringTypes.contactThresholdDistance]
            altSelection = paramsBlock.get(blockNames.ClusteringTypes.alternativeStructure, False)
        except KeyError as err:
            err.message += ": Need to provide mandatory parameter in clustering block"
            raise KeyError(err.message)
        if clusteringType == blockNames.ClusteringTypes.contacts:
            symmetries = paramsBlock.get(blockNames.ClusteringTypes.symmetries, [])

            thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
            thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
            return ContactsClustering(thresholdCalculator, resname,
                                      reportBaseFilename, columnOfReportFile,
                                      contactThresholdDistance, symmetries,
                                      altSelection=altSelection)
        elif clusteringType == blockNames.ClusteringTypes.lastSnapshot:

            return SequentialLastSnapshotClustering(resname, reportBaseFilename,
                                                    columnOfReportFile,
                                                    contactThresholdDistance)
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
            symmetries = paramsBlock.get(blockNames.ClusteringTypes.symmetries,[])
            thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
            thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
            try:
                similarityEvaluatorType = paramsBlock[blockNames.ClusteringTypes.similarityEvaluator]
            except KeyError:
                raise ValueError("No similarity Evaluator specified!!")
            similarityBuilder = similarityEvaluatorBuilder()
            similarityEvaluator = similarityBuilder.build(similarityEvaluatorType)
            return ContactMapAccumulativeClustering(thresholdCalculator, similarityEvaluator, resname,
                                                    reportBaseFilename, columnOfReportFile,
                                                    contactThresholdDistance, symmetries, altSelection)
        else:
            sys.exit("Unknown clustering method! Choices are: " +
                     str(clusteringTypes.CLUSTERING_TYPE_TO_STRING_DICTIONARY.values()))


class clusteringResultsParameters:
    """
        Helper object to pass parameters in the ContactMap affinity and agglomerative clustering
    """
    def __init__(self, pdb_list=[], metrics=[], contactmaps=[], contacts=[]):
            self.pdb_list = pdb_list
            self.metrics = metrics
            self.contactmaps = contactmaps
            self.elements = []
            self.ids = []
            self.newClusters = Clusters()
            self.contacts = contacts


class similarityEvaluatorBuilder:
    def build(self, similarityEvaluatorType):
        if similarityEvaluatorType == blockNames.ClusteringTypes.differenceDistance:
            return differenceDistanceEvaluator()
        elif similarityEvaluatorType == blockNames.ClusteringTypes.Jaccard:
            return JaccardEvaluator()
        elif similarityEvaluatorType == blockNames.ClusteringTypes.correlation:
            return correlationEvaluator()
        else:
            sys.exit("Unknown threshold calculator type! Choices are: " + str(clusteringTypes.SIMILARITY_TYPES_TO_STRING_DICTIONARY.values()))


class differenceDistanceEvaluator:
    """
        Evaluate the similarity of two contactMaps by calculating the ratio of
        the number of differences over the average of elements in the contacts
        maps
    """
    def isSimilarCluster(self, contactMap, clusterContactMap, symContactMapEvaluator):
        """
            Evaluate if two contactMaps are similar or not, return True if yes,
            False otherwise
        """
        return symContactMapEvaluator.evaluateDifferenceDistance(contactMap, clusterContactMap)


class JaccardEvaluator:
    """
        Evaluate the similarity of two contactMaps by calculating the Jaccard
        index, that is, the ratio between the intersection of the two contact
        maps and their union
    """
    def isSimilarCluster(self, contactMap, clusterContactMap, symContactMapEvaluator):
        """
            Evaluate if two contactMaps are similar or not, return True if yes,
            False otherwise
        """
        return symContactMapEvaluator.evaluateJaccard(contactMap, clusterContactMap)


class correlationEvaluator:
    """
        Evaluate the similarity of two contact maps by calculating their
        correlation
    """
    def isSimilarCluster(self, contactMap, clusterContactMap, symContactMapEvaluator):
        """
            Evaluate if two contactMaps are similar or not, return True if yes,
            False otherwise
        """
        return symContactMapEvaluator.evaluateCorrelation(contactMap, clusterContactMap)


# TODO: should it be a class method?
def clusterContactMaps(contactmaps, preferences=None):
    """
        Cluster multiple contact maps with the affinity propagation algorithm
        contactmaps [In] Array of contact maps
        preferences [In] Input value to the affinity propagation algorithm, it
        affects the number of clusters that will be generated
    """
    contactmaps = contactmaps.reshape((contactmaps.shape[0], -1))
    affinitypropagation = AffinityPropagation(damping=0.9,
                                              preference=preferences,
                                              verbose=False).fit(contactmaps)
    cluster_center_indices = affinitypropagation.cluster_centers_indices_
    labels = affinitypropagation.labels_
    return cluster_center_indices, np.array(labels)


def clusterAgglomerativeContactMaps(contactmaps, n_clusters):
    """
        Cluster multiple contact maps with the hierarchical clustering algorithm
        contactmaps [In] Array of contact maps
        n_clusters [In] Number of clusters that the clustering should generate
    """
    contactmaps = contactmaps.reshape((contactmaps.shape[0], -1))
    agglomerative = AgglomerativeClustering(n_clusters=n_clusters,
                                            linkage='complete').fit(contactmaps)
    labels = agglomerative.labels_
    clusters = np.unique(labels)
    return clusters, labels


def clusterKmeans(contactmaps, n_clusters):
    """
        Cluster  multiple contact maps with the kmeans clustering algorithm
        contactmaps [In] Array of contact maps
        n_clusters [In] Number of clusters that the clustering should generate
    """
    contactmaps = contactmaps.reshape((contactmaps.shape[0], -1))
    kmeans = KMeans(n_clusters=n_clusters).fit(contactmaps)
    center = kmeans.cluster_centers_[0]
    contactmaps -= center
    distances = contactmaps.sum(axis=1)
    cluster_center = abs(distances).argmin()
    return cluster_center


def getAllTrajectories(paths):
    """
        Find all the trajectory files in the paths specified

        :param paths: The path where to find the trajectories
        :type paths: str
        :returns: list -- A list with the names of all th trajectories in paths
    """
    files = []
    for path in paths:
        files += glob.glob(path)
    # sort the files obtained by glob by name, so that the results will be the
    # same on all computers
    return sorted(files)


def selectRandomCenter(cluster_members, metrics_weights):
    """
        Extract a center randomly from a cluster, weighted according to their
        metrics
    """
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


def processSnapshots(trajectories, reportBaseFilename, col,
                     contactThresholdDistance, resname, symmetryEvaluator):
    """
        Create list of pdb, contactMaps, metrics and contacts from a series of
        snapshots
    """
    pdb_list = []
    metrics = []
    contactmaps = []
    contacts = []
    for trajectory in trajectories:
        trajNum = utilities.getTrajNum(trajectory)
        snapshots = utilities.getSnapshots(trajectory, True)
        if reportBaseFilename:
            reportFilename = os.path.join(os.path.split(trajectory)[0],
                                          reportBaseFilename % trajNum)
            metricstraj = np.loadtxt(reportFilename, ndmin=2)
        else:
            metricstraj = np.zeros(len(snapshots))
        # Prepare data for clustering
        metrics.extend(metricstraj.tolist())
        for num, snapshot in enumerate(snapshots):
            pdb = atomset.PDB()
            pdb.initialise(snapshot, resname=str(resname))
            pdb_list.append(pdb)
            contactMap, contactnum = symmetryEvaluator.createContactMap(pdb, resname, contactThresholdDistance)
            contactmaps.append(contactMap)
            contacts.append(contactnum)
    return pdb_list, metrics, contactmaps, contacts
