import buildMSM
import sys
import MSMblocks

def main(controlFile):
    trajectoryFolder, trajectoryBasename, numClusters, lagtimes, numPCCA, itsOutput, numberOfITS, itsErrors, error_estimationCK, state_labels, outfile_fluxTPT = buildMSM.readParams(controlFile)

    prepareMSM = MSMblocks.PrepareMSM(numClusters, trajectoryFolder, trajectoryBasename)
    cl = prepareMSM.getClusteringObject()
    calculateMSM = MSMblocks.MSM(cl, lagtimes, 0, itsOutput, numberOfITS, itsErrors, error_estimationCK)
    lagtime = calculateMSM.calculateITS()
    calculateMSM.createMSM(lagtime)
    calculateMSM.check_connectivity()
    MSM_object = calculateMSM.getMSM_object()

if __name__ == "__main__":
    controlFile = sys.argv[1]
    main(controlFile)
