import os
import pickle
import utilities


def writeStructures(clusteringObject, listStructures, outputPath="cluster.pdb"):
    with open(clusteringObject, "rb") as f:
        clObject = pickle.load(f)
    nameStructure = os.path.splitext(outputPath)
    outputName = nameStructure[0]+'_%d'+nameStructure[1]
    path = os.path.split(outputName)
    pathToWrite = path[1]
    if path[0]:
        utilities.makeFolder(path[0])
        pathToWrite = os.path.join(path[0], path[1])

    if listStructures is None or len(listStructures) == 0: #If no listStructures, write all
        listStructures = range(len(clObject.clusters.clusters))

    for element in listStructures:
        print pathToWrite, pathToWrite%element
        clObject.clusters.clusters[element].writePDB(pathToWrite % element)
