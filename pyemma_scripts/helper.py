import os
import numpy as np
import cPickle

def writeClusterCenters(cl, outputFilename):
    np.savetxt(outputFilename, cl.clustercenters)

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

def saveMSM(MSM_object):
    """Save the MSM object to avoid having to run again
    the more computationally expensive part"""
    with open("MSM_object.pkl","w") as MSMfile:
        cPickle.dump(MSM_object, MSMfile, -1)

def saveClustering(clustering_object, clusteringFile):
    """Save the clustering object to avoid having to run again
    the more computationally expensive part"""
    with open(clusteringFile,"w") as clustering_file:
        cPickle.dump(clustering_object, clustering_file, -1)

def loadMSM(MSMFile):
    with open(MSMFile) as MSMfile:
        MSM_object = cPickle.load(MSMfile)
    return MSM_object

def loadClustering(clusteringFile):
    with open(clusteringFile,"r") as clustering_file:
        cl = cPickle.load(clustering_file)
    return cl
