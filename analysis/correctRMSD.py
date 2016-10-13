import atomset
import os
import glob
import numpy as np
from utilities import utilities


def getRMSD(traj, nativePDB, resname, symmetries):
    snapshots = utilities.getSnapshots(traj)

    rmsds = np.zeros(len(snapshots))

    for i, snapshot in enumerate(snapshots):
        snapshotPDB = atomset.PDB()
        snapshotPDB.initialise(snapshot, resname=resname)

        rmsds[i] = atomset.computeRMSD(nativePDB, snapshotPDB, symmetries)

    return rmsds




def extendReportWithRmsd(reportFile, rmsds):
    newShape = reportFile.shape
    newShape[1] += 1
    fixedReport = np.zeros(newShape)
    fixedReport[:,:-1] = reportFile
    fixedReport[:,-1] = rmsds
    return fixedReport


def main():
    folder = "."
    outputFilename = "fixedReport_%d"
    trajName = "*traj*.pdb"
    reportName = "*report_%d"
    resname = "K5Y"
    nativeFilename = "/gpfs/scratch/bsc72/bsc72755/adaptiveSampling/data/4K5Y/4K5Y_native.pdb"
    symmetries = {"C12:K5Y":"C13:K5Y", "C10:K5Y":"C17:K5Y", "C14:K5Y":"C16:K5Y", "C19:K5Y":"C20:K5Y", "C18:K5Y":"C21:K5Y"}
    rmsdColInReport = 5


    nativePDB = atomset.PDB()
    nativePDB.initialise(nativeFilename, resname=resname)

    utilities.generateReciprocalAtoms(symmetries)


    allFolders = os.listdir(folder)
    epochs = [epoch for epoch in allFolders if epoch.isdigit()]

    for epoch in epochs:
        os.chdir(epoch)
        allTrajs = glob.glob(trajName)

        for traj in allTrajs:
            rmsds = getRMSD(traj, nativePDB, resname, symmetries)
            trajNum = utilities.getTrajNum(traj)
            try:
                reportFilename = glob.glob(reportName%trajNum)[0]
            except IndexError:
                raise IndexError("File %s not found in folder %s"%(reportName%trajNum, epoch))

            reportFile = np.loadtxt(reportFilename, ndmin=2)

            if rmsdColInReport < reportFile.shape[1]:
                reportFile[:,rmsdColInReport] = rmsds
                fixedReport = reportFile
            else:
                fixedReport = extendReportWithRmsd(reportFile, rmsds)

            np.savetxt(outputFilename%trajNum, fixedReport, fmt='%.4f')

        os.chdir("..")

if __name__ == "__main__":
    main()
