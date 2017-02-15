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


def assignNewTrajecories(X, goldenMSMClusterCenters, lagtime, ntrajs=None, length=None):

    if ntrajs:
        indices = np.array(np.random.choice(range(len(X)), ntrajs))
    else:
        indices = range(len(X))

    Xsample = map(lambda x: X[x][:length], indices)

    assign = AssignCenters(goldenMSMClusterCenters)
    dTrajs = assign.assign(Xsample)
    counts = markov.estimateCountMatrix(dTrajs, goldenMSMClusterCenters.shape[0], lagtime)
    #secondMSM = msm.estimateMSM(dTrajs, lagtime)
    #counts = secondMSM.count_matrix_full
    counts += 1.0/counts.shape[0]
    transition = markov.buildTransitionMatrix(counts)
    return transition


def main(controlFile):
    trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, sampleSize, numRuns = readParams(controlFile)
    refTransition, refStationaryDist, lagtime = getTransitionMatrix(trajectoryFolder, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime)
    X = trajectories.loadCOMFiles(trajectoryFolder2, trajectoryBasename)

    goldenMSMClusterCenters = np.loadtxt("discretized/clusterCenters.dat")

    #np.random.seed(250793)

    entropies = []
    try:
        numberOfTrajs = range(5,20,5) + range(20, sampleSize, 20)

        #only trying different traj lengths if sampleSize is defined in control file
        shortestTrajSize = min([len(i) for i in X])
        lowerLimit = 2*lagtime
        upperLimit = shortestTrajSize
        #allTrajLengths = range(lowerLimit, upperLimit, 200)
        allTrajLengths = range(lowerLimit, 550, 50) + range(750, upperLimit, 250)
        #allTrajLengths = [None]
    except TypeError:
        numberOfTrajs = [None]
        allTrajLengths = [None]

    for length in allTrajLengths:
        if length: print "Working with trajectories of length: %d" % length
        lengthEntropies = []
        for ntrajs in numberOfTrajs:
            if not ntrajs % 25:
                print "Starting loop for sample of %d trajs" % ntrajs
            relativeEntropy = 0
            for j in range(numRuns):
                transitionMatrix = assignNewTrajecories(X, goldenMSMClusterCenters, lagtime, ntrajs, length)
                try:
                    s = markov.getRelativeEntropy(refStationaryDist, refTransition, transitionMatrix)
                    relativeEntropy += s
                except ValueError:
                    j -= 1
            lengthEntropies.append(relativeEntropy/float(numRuns))
        entropies.append(lengthEntropies)

    np.save("matrix_1.npy", entropies)

    for i, length in enumerate(allTrajLengths):
        for j, ntrajs  in enumerate(numberOfTrajs):
            print length, ntrajs, entropies[i][j]
        print ""

    if ntrajs and length:
        plt.figure(1)
        plt.imshow(entropies, interpolation="nearest", origin="lower", aspect="auto", extent=[numberOfTrajs[0], numberOfTrajs[-1], allTrajLengths[0], allTrajLengths[-1]])
        #plt.imshow(entropies, interpolation="nearest", extent=[numberOfTrajs[0], numberOfTrajs[-1], 800, 801])
        plt.colorbar()
        plt.figure(2)
        plt.imshow(entropies, interpolation="bilinear", origin="lower", aspect="auto",  extent=[numberOfTrajs[0], numberOfTrajs[-1], allTrajLengths[0], allTrajLengths[-1]])
        plt.colorbar()
        plt.show()
    elif ntrajs:
        print numberOfTrajs, entropies[0]
        plt.plot(numberOfTrajs, entropies[0])
        plt.show()
    else:
        print entropies

if __name__ == "__main__":
    controlFile = sys.argv[1]
    MSM_object = main(controlFile)
