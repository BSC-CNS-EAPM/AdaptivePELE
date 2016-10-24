import sys
import MSMblocks
import helper
import runMarkovChainModel as markov
import trajectories
from pyemma.coordinates.clustering import AssignCenters
import msm
import numpy as np
import matplotlib.pyplot as plt

import pdb as debug


def readParams(control_file):
    params = MSMblocks.readParams(control_file)
    trajectoryFolder = params["trajectoryFolder"]
    trajectoryFolder2 = params["trajectoryFolder2"]
    trajectoryBasename = params["trajectoryBasename"]
    numClusters = params["numClusters"]
    lagtimes = params["lagtimes"]
    itsOutput = params["itsOutput"]
    numberOfITS = params["numberOfITS"]
    itsErrors = params["itsErrors"]
    lagtime = params.get("lagtime", 0)
    sampleSize = params.get("sampleSize", None)
    numRuns = params.get("numRuns", 1)
    return trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, sampleSize, numRuns


def getTransitionMatrix(trajectoryFolder, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime):
    prepareMSM = MSMblocks.PrepareMSM(numClusters, trajectoryFolder, trajectoryBasename)
    cl = prepareMSM.getClusteringObject()
    calculateMSM = MSMblocks.MSM(cl, lagtimes, 0, itsOutput, numberOfITS, itsErrors, None)
    if not lagtime:
        lagtime = calculateMSM.calculateITS()
    calculateMSM.createMSM(lagtime)
    calculateMSM.check_connectivity()
    MSM_object = calculateMSM.getMSM_object()
    helper.saveMSM(MSM_object)
    counts = MSM_object.count_matrix_full
    counts += 1.0/numClusters
    transition = markov.buildTransitionMatrix(counts)

    return transition, MSM_object.stationary_distribution, lagtime


def assignNewTrajecories(X, lagtime, sampleSize=None):

    if sampleSize:
        indices = np.array(np.random.choice(range(len(X)), sampleSize))
        Xsample = map(lambda x: X[x], indices)
    else:
        Xsample = X
    assign = AssignCenters("discretized/clusterCenters.dat")
    dTrajs = assign.assign(Xsample)
    secondMSM = msm.estimateMSM(dTrajs, lagtime)
    counts = secondMSM.count_matrix_full
    counts += 1.0/counts.shape[0]
    transition = markov.buildTransitionMatrix(counts)
    return transition


def main(controlFile):
    trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, sampleSize, numRuns = readParams(controlFile)
    refTransition, refStationaryDist, lagtime = getTransitionMatrix(trajectoryFolder, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime)
    X = trajectories.loadCOMFiles(trajectoryFolder2, trajectoryBasename)
    entropies = []
    try:
        allSizes = range(5, sampleSize, 5)
    except TypeError:
        allSizes = [None]
    for i in allSizes:
        if not i % 100:
            print "Starting loop for sample of size %d" % i
        relativeEntropy = 0
        for j in range(numRuns):
            transitionMatrix = assignNewTrajecories(X, lagtime, i)
            try:
                relativeEntropy += markov.getRelativeEntropy(refStationaryDist, refTransition, transitionMatrix)
            except:
                j -= 1
        entropies.append(relativeEntropy/float(numRuns))
    if i:
        plt.plot(allSizes, entropies)
        plt.show()
    else:
        print entropies

if __name__ == "__main__":
    controlFile = sys.argv[1]
    MSM_object = main(controlFile)
