import os
import shutil
import numpy as np


def cleanup(tmpFolder):
    if os.path.exists(tmpFolder):
        shutil.rmtree(tmpFolder)


def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)


def getSnapshots(trajectoryFile, verbose=False):
    inputFile = open(trajectoryFile, "r")
    inputFileContent = inputFile.read()

    snapshots = inputFileContent.split("ENDMDL")[:-1]
    if not verbose:
        return snapshots

    remarkInfo = "REMARK 000 File created using PELE++\nREMARK source            : %s\nREMARK original model nr : %d\nREMARK First snapshot is 1, not 0 (as opposed to report)\n"
    snapshotsWithInfo = [remarkInfo % (file, i+1)+snapshot for i, snapshot in enumerate(snapshots)]
    return snapshotsWithInfo


def getTrajNum(trajFilename):
    return int(trajFilename.split("_")[-1][:-4])


def calculateContactMapEigen(contactMap):
    nLig, nCA = contactMap.shape
    extendedCM = np.zeros((nLig+nCA, nLig+nCA))
    extendedCM[nCA:, :nCA] = contactMap
    extendedCM[:nCA, nCA:] = contactMap.T
    assert (extendedCM == extendedCM.T).all(), "Extended ContactMap not symmetric"
    eiv, eic = np.linalg.eigh(extendedCM)
    return eiv, eic


def generateReciprocalAtoms(pairsDictionary):
    tmpDict = {}
    for key,val in pairsDictionary.iteritems():
        tmpDict[val] = key
    return pairsDictionary.update(tmpDict)

def assertSymmetriesDict(pairsDictionary, PDB):
    for key in pairsDictionary:
        assert key in PDB.atoms, "Symmetry atom %s not found in initial structure" % key
