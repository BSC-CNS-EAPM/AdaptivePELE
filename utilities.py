import os
import shutil
import pickle
from clustering import clustering

def cleanup(tmpFolder):
    if os.path.exists(tmpFolder):
        shutil.rmtree(tmpFolder)

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

def writeStructures(clusteringObject, listStructures, outputPath="cluster.pdb"):
    with open(clusteringObject, "rb") as f:
        clObject = pickle.load(f)
    nameStructure = os.path.splitext(outputPath)
    outputName = nameStructure[0]+'_%d'+nameStructure[1]
    path = os.path.split(outputName)
    if path[0]:
        makeFolder(path[0])
    for element in listStructures:
        clObject.clusters.clusters[element].writePDB(path[1] % element)
