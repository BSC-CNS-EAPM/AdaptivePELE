import atomset
import os
import glob
import numpy as np
from utilities import utilities
import argparse
import json

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

def parseArguments():
    desc = "Program that fixes RMSD symmetries of a PELE report file."\
            "Control file is a JSON file that contains \"resname\", \"native\", "\
            "symmetries, and, optionally, the column to substitute in report. "\
            "Example of content:"\
            "{"\
            "\"resname\" : \"K5Y\","\
            "\"native\" : \"native.pdb\","\
            "\"symmetries\" : {\"4122:C12:K5Y\":\"4123:C13:K5Y\", \"4120:C10:K5Y\":\"4127:C17:K5Y\"},"\
            "\"column\" = 5"\
            "}"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("controlFile", type=str, help="Control File name")
    args = parser.parse_args()

    return  args.controlFile, args.column, args.threshold, args.stepsPerEpoch, args.seq

def readControlFile(controlFile):
    jsonFile = open(controlFile, 'r').read()
    parsedJSON = json.loads(jsonFile)
    resname = parsedJSON["resname"]
    nativeFilename = parsedJSON["native"]
    symmetries = parsedJSON["symmetries"]
    rmsdColInReport = parsedJSON.get("column")
    if not rmsdColInReport:
        #append to the end
        rmsdColInReport = -1

    return resname, nativeFilename, symmetries, rmsdColInReport

def main(controlFile):
    #Constants
    folder = "."
    outputFilename = "fixedReport_%d"
    trajName = "*traj*.pdb"
    reportName = "*report_%d"
    #end constants

    resname, nativeFilename, symmetries, rmsdColInReport = readControlFile(controlFile)


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

            if rmsdColInReport > 0 and rmsdColInReport < reportFile.shape[1]:
                reportFile[:,rmsdColInReport] = rmsds
                fixedReport = reportFile
            else:
                fixedReport = extendReportWithRmsd(reportFile, rmsds)

            np.savetxt(outputFilename%trajNum, fixedReport, fmt='%.4f')

        os.chdir("..")

if __name__ == "__main__":
    controlFile = parseArguments()
    main()
