# coding: utf-8
import sys
reload(sys)
sys.setdefaultencoding('utf-8')
import os
import argparse
import glob
from AdaptivePELE.atomset import atomset
import re
import numpy as np

DIRECTORY = '%s/extractedCoordinates'
BASEOUTPUTFILENAME = 'coord'
LIGANDTRAJECTORYBASENAME = "traj_ligand_%s.pdb"


def parseArguments():
    desc = "Extracts coordinates in <currentFolder>/extractedCoordinates/coord*.\
            It either extracts the resname COM coordinates or those of an atomId, depending on the input.\
            It automatically detects whether it is an adaptive or a sequential PELE run looking for folders\
            with numeric names."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-f", "--folderWithTrajs", default = ".",
                        help="Folder with trajectories (or epochs)")
    parser.add_argument("-atomId", default="", help="serial:atomName:resname, e.g. 2048:C1:AIN")
    parser.add_argument("-resname", default="", help="Ligand resname")
    parser.add_argument("-s", "--enforceSequential", action="store_true", help="Enforce the consideration as sequential run")
    parser.add_argument("-w", "--writeLigandTrajectory", action="store_true", help="It writes a traj_ligand_XXX.pdb file with the ligand coordinates. The user must delete the original trajectory (if wanted)")
    # parser.add_argument("-f", nargs='+', help="Files to get coordinates")
    args = parser.parse_args()
    return args.folderWithTrajs, args.atomId, args.resname, args.enforceSequential, args.writeLigandTrajectory


def loadAllResnameAtomsInPdb(filename, resname):
    fileContent = open(filename).read()
    fileContent = fileContent.split('ENDMDL')
    prunedFileContent = []
    for snapshot in fileContent:
        snapshot = snapshot.split('\n')
        prunedSnapshot =  [line for line in snapshot if line[17:20] == resname]
        prunedFileContent.append("\n".join(prunedSnapshot))
    return prunedFileContent

def extractFilenumber(filename):
    last = filename.rfind('.')
    first = filename.rfind('_')
    number = re.sub("[^0-9]", "", filename[first+1:last])
    return number

def getOutputFilename(directory, filename, baseOutputFilename):
    filenumber = extractFilenumber(filename)
    return os.path.join(directory, baseOutputFilename+"_"+filenumber+".dat")

def getPDBCOM(allCoordinates, resname):
    COMs = []
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=resname, heavyAtoms=True)
        COMs.append(pdb.extractCOM())
    return COMs

def getAtomCoord(allCoordinates, resname, atomId):
    coords = []
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=resname, heavyAtoms=True)
        coords.append(pdb.getAtom(atomId).getAtomCoords())
    return coords

def writeToFile(COMs, outputFilename):
    with open(outputFilename, 'w') as f:
        for i, line in enumerate(COMs):
            f.write(str(i) + ' ')
            for i in range(len(line) - 1):
                f.write(str(line[i]) + ' ')
            f.write(str(line[-1]) + '\n')

def writeFilenameExtractedCoordinates(filename, resname, atomId, pathFolder, writeLigandTrajectory):
    print filename
    allCoordinates = loadAllResnameAtomsInPdb(filename, resname)
    if writeLigandTrajectory:
        outputFilename = os.path.join(pathFolder, LIGANDTRAJECTORYBASENAME%extractFilenumber(filename))
        with open(outputFilename, 'w') as f:
            f.write("\nENDMDL\n".join(allCoordinates))

    # because of the way it's split, the last element is empty
    if atomId != "":
        coords = getAtomCoord(allCoordinates[:-1], resname, atomId)
    else:
        coords = getPDBCOM(allCoordinates[:-1], resname)

    outputFilename = getOutputFilename(DIRECTORY, filename,
                                       BASEOUTPUTFILENAME)
    writeToFile(coords, outputFilename % pathFolder)

def writeFilenamesExtractedCoordinates(pathFolder, resname, atomId, writeLigandTrajectory):
    if not os.path.exists(DIRECTORY % pathFolder):
        os.makedirs(DIRECTORY % pathFolder)

    originalPDBfiles = glob.glob(pathFolder+'/*trajec*.pdb')
    for filename in originalPDBfiles:
        writeFilenameExtractedCoordinates(filename, resname, atomId, pathFolder, writeLigandTrajectory)

def parseResname(atomId, resname):
    if atomId == "" and resname == "":
        sys.exit("Either resname or atomId should be provided")
    elif resname == "":
        resname = atomId.split(":")[-1] #the atom Id last element is the resname
    elif atomId != "":
        resnameInAtomId = atomId.split(":")[-1] #the atom Id last element is the resname
        if resnameInAtomId != resname:
            sys.exit("Residue name in resname and atomId do not match!")
    return resname

def main(folder=".", atomId="", resname="", enforceSequential=False, writeLigandTrajectory=True):
    resname = parseResname(atomId, resname)
    folderWithTrajs = folder

    if enforceSequential:
        folders = [""]
    else:
        allFolders = os.listdir(folderWithTrajs)
        folders = [epoch for epoch in allFolders if epoch.isdigit()]
        if len(folders) == 0:
            folders = [""]

    for folder in folders:
        pathFolder = os.path.join(folderWithTrajs, folder)
        writeFilenamesExtractedCoordinates(pathFolder, resname, atomId, writeLigandTrajectory)


if __name__ == "__main__":
    folder, atomId, resname, enforceSequential, writeLigandTrajectory = parseArguments()
    main(folder, atomId, resname, enforceSequential, writeLigandTrajectory)
