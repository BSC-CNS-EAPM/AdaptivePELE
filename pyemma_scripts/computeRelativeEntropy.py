import os
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
    stride = params.get("stride", 1)
    return trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, sampleSize, numRuns, stride


def getTransitionMatrix(trajectoryFolder, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, stride):
    prepareMSM = MSMblocks.PrepareMSM(numClusters, trajectoryFolder, trajectoryBasename, stride)
    cl = prepareMSM.getClusteringObject()
    calculateMSM = MSMblocks.MSM(cl, lagtimes, 0, itsOutput, numberOfITS, itsErrors, None, dtrajs=prepareMSM.dtrajs)
    if not lagtime:
        lagtime = calculateMSM.calculateITS()
    calculateMSM.createMSM(lagtime)
    calculateMSM.check_connectivity()
    MSM_object = calculateMSM.getMSM_object()
    helper.saveMSM(MSM_object)
    counts = MSM_object.count_matrix_full
    counts += 1.0/numClusters
    transition = markov.buildTransitionMatrix(counts)

    fullStationaryDistribution = np.zeros(MSM_object.nstates_full)
    active = MSM_object.active_set
    fullStationaryDistribution[active] = MSM_object.stationary_distribution
    return transition, fullStationaryDistribution, lagtime

def makeRandomSampleOfNtrajs(X, ntrajs=None, length=None):
    if ntrajs:
        indices = np.array(np.random.choice(range(len(X)), ntrajs))
    else:
        indices = range(len(X))

    Xsample = map(lambda x: X[x][:length,:], indices)
    return Xsample

def assignNewTrajecories(Xsample, goldenMSMClusterCenters, lagtime):
    assign = AssignCenters(goldenMSMClusterCenters)
    dTrajs = assign.assign(Xsample)
    counts = markov.estimateCountMatrix(dTrajs, goldenMSMClusterCenters.shape[0], lagtime)
    #secondMSM = msm.estimateMSM(dTrajs, lagtime)
    #counts = secondMSM.count_matrix_full
    counts += 1.0/counts.shape[0]
    transition = markov.buildTransitionMatrix(counts)
    return transition

def plotIsocostLines(extent, allTrajLengths, numberOfTrajs, steps=10):
    minCost = allTrajLengths[0]*numberOfTrajs[0]
    maxCost = allTrajLengths[-1]*numberOfTrajs[-1]
    d = (maxCost - minCost) / steps
    for cost in np.arange(minCost, maxCost, d):
        x = np.arange(extent[0], extent[1], 1)
        y = cost / x
        plt.plot(x,y, color="black")


def main(controlFile):
    trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, sampleSize, numRuns, stride = readParams(controlFile)
    refTransition, refStationaryDist, lagtime = getTransitionMatrix(trajectoryFolder, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, stride)


    goldenMSMClusterCenters = np.loadtxt("discretized/clusterCenters.dat")

    #np.random.seed(250793)

    seq = True
    entropies = []

    if seq:
        try:
            X = trajectories.loadCOMFiles(trajectoryFolder2, trajectoryBasename)

            dTrajs = 50
            numberOfTrajs = range(63, 512, dTrajs)
            numberOfTrajs = range(50, 502, dTrajs)
            # numberOfTrajs = range(50, sampleSize, 50)

            #only trying different traj lengths if sampleSize is defined in control file
            shortestTrajSize = min([len(i) for i in X])
            lowerLimit = 2*lagtime
            upperLimit = shortestTrajSize
            upperLimit = 1001
            dTrajLengths = 100
            allTrajLengths = range(lowerLimit, upperLimit, dTrajLengths)
            #allTrajLengths = range(lowerLimit, 401, 100) 
            #allTrajLengths = [None]
        except TypeError:
            numberOfTrajs = [None]
            allTrajLengths = [None]
    else:
        # epochFolders = [int(i) for i in os.listdir(trajectoryFolder2) if i.isdigit()]
        dTrajs = 1
        epochFolders = range(3)
        epochFolders.sort()
        # import pdb as debug
        # debug.set_trace()
        lowerLimit = 200+50
        numberOfTrajs = epochFolders
        upperLimit = 401
        dTrajLengths = 50
        allTrajLengths = range(lowerLimit, upperLimit, dTrajLengths)


    for length in allTrajLengths:
        if length: print "Working with trajectories of length: %d" % length
        lengthEntropies = []
        if not seq:
            X = []
        for ntrajs in numberOfTrajs:

            if not seq:
                currentEpochFolder = os.path.join(trajectoryFolder2, str(ntrajs))
                currentX = trajectories.loadCOMFiles(currentEpochFolder, trajectoryBasename)
                X.extend(currentX) 

            if not ntrajs % 25:
                print "Starting loop for sample of %d trajs" % ntrajs
            relativeEntropy = 0
            for j in range(numRuns):
                if seq:
                    Xsample = makeRandomSampleOfNtrajs(X, ntrajs, length)
                else:
                    Xsample = map(lambda x: x[:length,:], X)

                transitionMatrix = assignNewTrajecories(Xsample, goldenMSMClusterCenters, lagtime)
                try:
                    s = markov.getRelativeEntropy(refStationaryDist, refTransition, transitionMatrix)
                    print s
                    relativeEntropy += s
                except ValueError:
                    j -= 1
            lengthEntropies.append(relativeEntropy/float(numRuns))
            print relativeEntropy, lengthEntropies
        entropies.append(lengthEntropies)

    np.save("matrix_adaptive.npy", entropies)
    #entropies = np.load("matrix_adaptive.npy")
    entropies = np.log10(entropies)

    for i, length in enumerate(allTrajLengths):
        for j, ntrajs  in enumerate(numberOfTrajs):
            print length, ntrajs, entropies[i][j]
        print ""

    if ntrajs and length:
        plt.figure(1)
        if seq:
            extent = [numberOfTrajs[0] - dTrajs/2, numberOfTrajs[-1] + dTrajs/2,
                    allTrajLengths[0] - dTrajLengths/2, allTrajLengths[-1] + dTrajLengths/2]
            #plot isocost lines
            plotIsocostLines(extent, allTrajLengths, numberOfTrajs, 11)
            plt.imshow(entropies, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
        else:
            plt.imshow(entropies, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
        import os
        cwd = os.getcwd()
        cwd = cwd.replace("/", "_")
        #plt.save(cwd + ".eps")
        #plt.show()
        #plt.imshow(entropies, interpolation="nearest", extent=[numberOfTrajs[0], numberOfTrajs[-1], 800, 801])
        plt.colorbar()
        plt.savefig(cwd + ".eps")
        plt.figure(2)
        if seq:
            plotIsocostLines(extent, allTrajLengths, numberOfTrajs, 11)
            plt.imshow(entropies, interpolation="bilinear", origin="lower", aspect="auto",  extent=extent)
        else:
            plt.imshow(entropies, interpolation="bilinear", origin="lower", aspect="auto",  extent=extent)
        plt.colorbar()
        import os
        cwd = os.getcwd()
        cwd = cwd.replace("/", "_")
        #plt.save(cwd + "2.eps")
        plt.savefig(cwd + "2.eps")
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
