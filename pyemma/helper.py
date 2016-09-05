import os
import numpy as np
import cPickle

def writeClusterCenters(cl, outputFilename):
    np.savetxt(outputFilename, cl.clustercenters)

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

def saveMSM(MSM_object, clustering_object):
    """Save the MSM object and clustering objects to avoid having to run again
    the more computationally expensive part"""
    with open("MSM_object.pkl","w") as MSMfile:
        cPickle.dump(MSM_object, MSMfile, -1)
    with open("clustering_object.pkl","w") as clustering_file:
        cPickle.dump(clustering_object, clustering_file, -1)

def loadMSM():
    with open("MSM_object.pkl","r") as MSMfile:
        MSM_object = cPickle.load(MSMfile)
    with open("clustering_object.pkl","r") as clustering_file:
        cl = cPickle.load(clustering_file)
    return MSM_object,cl
