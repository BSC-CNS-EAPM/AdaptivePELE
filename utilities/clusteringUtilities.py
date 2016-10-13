import os
import pickle

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
