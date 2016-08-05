import sys
import atomset
import glob
import numpy as np
import os

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
    def __init__(self, pdb, thresholdRadius, contacts=0, metric=0):
        self.pdb = pdb
        self.elements = 1
        self.threshold = thresholdRadius
        self.contacts = contacts
        self.metric = metric
    def addElement(self, metric):
        self.elements += 1
        self.metric = min(metric, self.metric) #this for the moment
    def printCluster(self, verbose=False):
        if verbose:
            print self.pdb.printAtoms()
        print "Elements: ", self.elements
        print "Radius threshold: ", self.threshold
        print "Min energy: ", self.metric
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
        cluster = Cluster (pdb, threshold, contacts, metric)
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

def getAllTrajectories(paths):
    files = []
    for path in paths:
        files += glob.glob(path)
    return files

def getSnapshots(trajectoryFile, verbose=False):
    return atomset.getPDBSnapshots(trajectoryFile, verbose)

def getTrajNum(trajFilename):
    return int(trajFilename.split("_")[-1][:-4])

