import argparse
from spawning import spawning, densitycalculator
from clustering import clustering
from atomset import atomset
from utilities import utilities
from scipy.sparse.csgraph import shortest_path
import numpy as np


class PathwayError(Exception):
    pass


def parseArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('clusteringObj', type=str)
    parser.add_argument('col', type=int)
    parser.add_argument('ntrajs', type=int)
    parser.add_argument('threshold', type=int)
    parser.add_argument('pathwayFilename', type=str, default="pathway.pdb", nargs='?')
    args = parser.parse_args()
    return args


def getOptimalCluster(clusteringObj, col):
    optimalMetric = 1000000
    optimalIndex = 0
    for i, cluster in enumerate(clusteringObj.clusters.clusters):
        if cluster.getMetricFromColumn(col) < optimalMetric:
            optimalMetric = cluster.getMetricFromColumn(col)
            optimalIndex = i
    return optimalIndex


def createNetworkMatrix(clusteringObj, threshold):
    n = clusteringObj.clusters.getNumberClusters()
    matrix = np.zeros((n, n))
    for i, cluster in enumerate(clusteringObj.clusters.clusters):
        for j in range(i, n):
            if i == j:
                continue
            cluster2 = clusteringObj.getCluster(j)
            if atomset.computeSquaredCentroidDifference(cluster.pdb, cluster2.pdb) > threshold:
                matrix[i, j] = matrix[j, i] = np.inf
            else:
                matrix[i, j] = matrix[j, i] = clusteringObj.RMSDCalculator.computeRMSD(cluster.pdb, cluster2.pdb)
    return matrix


def obtainShortestPath(distanceMatrix):
    foo, predecessors = shortest_path(distanceMatrix, return_predecessors=True)
    return predecessors


def createPathway(initial_cluster, final_cluster, predecessors):
    pathway = [final_cluster]
    i, j = initial_cluster, final_cluster
    while (i != j) and j >= 0:
        j = predecessors[i, j]
        pathway.insert(0, j)
    if pathway[0] < 0:
        # The network defined by the distanceMatrix is not connected
        raise PathwayError("Unable to find a path since network is not "
                           "connected, use a bigger threshold")

    return pathway


def writePathwayTrajectory(ClOrd, pathway, filename):
    pathwayFile = open(filename, "w")
    pathwayFile.write("REMARK 000 File created using PELE++\n")
    pathwayFile.write("REMARK 000 Pathway trajectory created using findPathway program\n")
    pathwayFile.write("REMARK 000 List of cluster belonging to the pathway %s\n" % ' '.join(map(str, pathway)))
    for i, step_cluster in enumerate(pathway):
        cluster = ClOrd.clusters.clusters[step_cluster]
        pathwayFile.write("MODEL %d\n" % (i+1))
        pdbStr = cluster.pdb.pdb
        pdbList = pdbStr.split("\n")
        for line in pdbList:
            line = line.strip()
            # Avoid writing previous REMARK block
            if line.startswith("REMARK ") or line.startswith("MODEL "):
                continue
            elif line:
                pathwayFile.write(line+"\n")
        pathwayFile.write("ENDMDL\n")
    pathwayFile.close()


def main(args):

    # Parameters
    clusteringObj = utilities.readClusteringObject(args.clusteringObj)
    col = args.col
    pathwayFilename = args.pathwayFilename
    ntrajs = args.ntrajs
    threshold = args.threshold

    # use graph algorithm to establish a path
    initial_cluster = 0
    final_cluster = getOptimalCluster(clusteringObj, col)
    distanceMatrix = createNetworkMatrix(clusteringObj, threshold)
    predecessors = obtainShortestPath(distanceMatrix)
    pathway = createPathway(initial_cluster, final_cluster, predecessors)
    print "Pathway clusters:"
    print pathway

    # write pathway into a single trajectory
    writePathwayTrajectory(clusteringObj, pathway, pathwayFilename)

    # create clustering object with only the pathway clusters
    ClPath = clustering.Clustering()
    ClPath.clusters.clusters = map(lambda x: clusteringObj.clusters.clusters[x],
                                   pathway)

    # spawning along the trajectory
    spawningParams = spawning.SpawningParams()
    densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
    densityCalculator = densityCalculatorBuilder.build({})
    spawningPathway = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
    # Set a least 1 processors from the extrems of the path
    degeneracies = spawningPathway.calculate(ClPath.clusters.clusters,
                                              ntrajs-2, spawningParams)
    degeneracies[0] += 1
    degeneracies[-1] += 1
    print "degeneracies over pathway:"
    print degeneracies
    print ""

if __name__ == "__main__":
    args = parseArgs()
    main(args)
