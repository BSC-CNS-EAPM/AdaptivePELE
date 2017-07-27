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
import shutil

class Constants:
    def __init__(self):
        self.extractedTrajectoryFolder = "%s/extractedCoordinates"
        self.baseExtractedTrajectoryName = "coord_"
        self.reportName = '*report_'
        self.outputTrajectoryFolder = "%s/repeatedExtractedCoordinates"
        self.ligandTrajectoryBasename = "traj_ligand_%s.pdb"
        self.gatherTrajsFolder = "allTrajs"
        self.gatherTrajsFilename = os.path.join(self.gatherTrajsFolder, "traj_%s_%s.dat")


def parseArguments():
    desc = "Program that extracts residue coordinates for a posterior MSM analysis.\
            It either extracts the resname COM coordinates or those of an atomId, depending on the input.\
            It then fills the rejected steps, which is not done by PELE.\
            Finally, trajectories are gathered together in the same allTrajs folder.\
            It automatically detects whether it is an adaptive or a sequential PELE run by looking for folders\
            with numeric names."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-f", "--folderWithTrajs", default = ".",
                        help="Folder with trajectories (or epochs)")
    #parser.add_argument("-atomId", action=AtomIdAction, help="serial:atomName:resname, e.g. 2048:C1:AIN")
    parser.add_argument("-atomIds", nargs='*', help="serial:atomName:resname, e.g. 2048:C1:AIN. May contain more than one atomId")
    parser.add_argument("-resname", default="", help="Ligand resname")
    parser.add_argument("-s", "--enforceSequential", action="store_true", help="Force the consideration as sequential run (non-adaptive)")
    parser.add_argument("--setNum", type=int, default=0, help="Sets the number to appear in gathered trajectory in order to avoid clashes between different sequential runs. Ignored in adaptive runs.")
    parser.add_argument("-w", "--writeLigandTrajectory", action="store_true", help="It writes a traj_ligand_XXX.pdb file with the ligand coordinates. The user must delete the original trajectory (if wanted)")
    parser.add_argument("-t", "--totalSteps", type=int, default=0, help="Total number of steps in traj. Equivalent to epoch length in adaptive runs")
    # parser.add_argument("-f", nargs='+', help="Files to get coordinates")
    args = parser.parse_args()

    return args.folderWithTrajs, args.atomIds, args.resname, args.enforceSequential, args.writeLigandTrajectory, args.totalSteps, args.setNum

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
    return os.path.join(directory, baseOutputFilename+filenumber+".dat")

def getPDBCOM(allCoordinates, resname):
    COMs = []
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=resname, heavyAtoms=True)
        COMs.append(pdb.extractCOM())
    return COMs

def getAtomCoord(allCoordinates, resname, atomIds):
    coords = []
    #If ever need to speed this up, build a Trajectory class that inherits from PDB
    #and loads the atom according to the position in the snapshot, rather than looking
    #for the atom
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=resname, heavyAtoms=True)
        snapshotcoords = []
        for atomId in atomIds:
            snapshotcoords.extend(pdb.getAtom(atomId).getAtomCoords())
        coords.append(snapshotcoords)
    return coords

def writeToFile(COMs, outputFilename):
    with open(outputFilename, 'w') as f:
        for i, line in enumerate(COMs):
            f.write(str(i) + ' ')
            for i in range(len(line) - 1):
                f.write(str(line[i]) + ' ')
            f.write(str(line[-1]) + '\n')

def writeFilenameExtractedCoordinates(filename, resname, atomIds, pathFolder, writeLigandTrajectory, constants):
    allCoordinates = loadAllResnameAtomsInPdb(filename, resname)
    if writeLigandTrajectory:
        outputFilename = os.path.join(pathFolder, constants.ligandTrajectoryBasename%extractFilenumber(filename))
        with open(outputFilename, 'w') as f:
            f.write("\nENDMDL\n".join(allCoordinates))

    # because of the way it's split, the last element is empty
    if not atomIds is None:
        coords = getAtomCoord(allCoordinates[:-1], resname, atomIds)
    else:
        coords = getPDBCOM(allCoordinates[:-1], resname)

    outputFilename = getOutputFilename(constants.extractedTrajectoryFolder, filename,
                                       constants.baseExtractedTrajectoryName)
    writeToFile(coords, outputFilename % pathFolder)

def writeFilenamesExtractedCoordinates(pathFolder, resname, atomIds, writeLigandTrajectory, constants):
    if not os.path.exists(constants.extractedTrajectoryFolder % pathFolder):
        os.makedirs(constants.extractedTrajectoryFolder % pathFolder)

    originalPDBfiles = glob.glob(pathFolder+'/*traj*.pdb')
    for filename in originalPDBfiles:
        writeFilenameExtractedCoordinates(filename, resname, atomIds, pathFolder, writeLigandTrajectory, constants)

def parseResname(atomIds, resname):
    if not atomIds is None and len(atomIds) > 0:
        differentResnames = set([atomId.split(":")[-1] for atomId in atomIds])

        if len(differentResnames) > 1:
            sys.exit("Error! Different resnames provided in atomIds!")
        elif len(differentResnames) == 1:
            extractedResname = differentResnames.pop()

    if (atomIds is None or len(atomIds) == 0) and resname == "":
        sys.exit("Either resname or atomId should be provided")
    elif resname == "":
        resname = extractedResname #the atom Id last element is the resname
    elif not atomIds is None and len(atomIds) > 0:
        if extractedResname != resname:
            sys.exit("Residue name in resname and atomId do not match!")
    return resname

def buildFullTrajectory(steps, trajectory, totalSteps, inputTrajectory):
    completeTrajectory = []
    counter = 0
    if len(trajectory) > 0:
        sthWrongInTraj = False
        print inputTrajectory
        for i in range(len(trajectory) - 1):
            try:
                repeated = steps[i+1] - steps[i]
            except IndexError:
                print "sth wrong in trajectory %s. This is likely to disagreement between report and trajecotry files. Please, fix it manually"%inputTrajectory
                sthWrongInTraj = True
                break

            for j in range(repeated):
                snapshot = trajectory[i].split()
                snapshot[0] = str(counter)
                snapshot = ' '.join(snapshot)
                completeTrajectory.append(snapshot)
                counter += 1

        if sthWrongInTraj:
            return completeTrajectory

        if totalSteps == 0:
            iterations = range(1)
        else:
            iterations = range(totalSteps + 1 - counter)

        for i in iterations:
            snapshot = trajectory[-1].split()
            snapshot[0] = str(counter)
            snapshot = ' '.join(snapshot)
            completeTrajectory.append(snapshot)
            counter += 1

    return completeTrajectory

def repeatExtractedSnapshotsInTrajectory(inputTrajectory, constants):
    extractedTrajFolder, trajFilename = os.path.split(inputTrajectory)
    trajectoryNumber = re.sub('\.dat$', '', trajFilename)
    trajectoryNumber = re.sub(constants.baseExtractedTrajectoryName, '', trajectoryNumber)

    origDataFolder = re.sub(constants.extractedTrajectoryFolder%"", "", extractedTrajFolder)
    try:
        reportFile = glob.glob(os.path.join(origDataFolder, constants.reportName + trajectoryNumber))[0]
    except:
        print "folder", origDataFolder
        sys.exit("Couldn't find file that matches: %s"%os.path.join(origDataFolder, constants.reportName + trajectoryNumber))

    with open(inputTrajectory) as f:
        trajectory = f.read().splitlines()

    acceptedSteps = np.loadtxt(reportFile, dtype='int', comments='#', usecols=(1,))

    fullTrajectory = buildFullTrajectory(acceptedSteps, trajectory, totalSteps, inputTrajectory)

    if len(fullTrajectory) > 0:
        outputFilename = os.path.join(constants.outputTrajectoryFolder%origDataFolder, constants.baseExtractedTrajectoryName + trajectoryNumber + '.dat')
        outputFile = open(outputFilename, 'w')
        for snapshot in fullTrajectory:
            outputFile.write("%s\n" % snapshot)
        outputFile.close()

def repeatExtractedSnapshotsInFolder(folder, constants, totalSteps):
    inputTrajectoryFolder = constants.extractedTrajectoryFolder % folder
    outputTrajectoryFolder = constants.outputTrajectoryFolder % folder

    if not os.path.exists(outputTrajectoryFolder):
        os.makedirs(outputTrajectoryFolder)

    inputTrajectories = glob.glob(os.path.join(inputTrajectoryFolder, constants.baseExtractedTrajectoryName + '*'))
    for inputTrajectory in inputTrajectories:
        repeatExtractedSnapshotsInTrajectory(inputTrajectory, constants)

def makeGatheredTrajsFolder(constants):
    if not os.path.exists(constants.gatherTrajsFolder):
        os.makedirs(constants.gatherTrajsFolder)

def gatherTrajs(constants, folder, setNum):
    trajectoriesFilenames = os.path.join(constants.outputTrajectoryFolder%folder, constants.baseExtractedTrajectoryName + "*")
    trajectories = glob.glob(trajectoriesFilenames)
    for inputTrajectory in trajectories:
        trajectoryNumber = extractFilenumber(os.path.split(inputTrajectory)[1])
        if folder != ".": #if not sequential
            setNum = folder
        shutil.copyfile(inputTrajectory, constants.gatherTrajsFilename%(setNum, trajectoryNumber))

def main(folder=".", atomIds=[""], resname="", totalSteps=0, enforceSequential=0, writeLigandTrajectory=True, setNum=0):
    constants = Constants()

    resname = parseResname(atomIds, resname)

    folderWithTrajs = folder

    makeGatheredTrajsFolder(constants)

    #change atomId for list

    if enforceSequential:
        folders = ["."]
    else:
        allFolders = os.listdir(folderWithTrajs)
        folders = [epoch for epoch in allFolders if epoch.isdigit()]
        if len(folders) == 0:
            folders = ["."]

    for folder in folders:
        pathFolder = os.path.join(folderWithTrajs, folder)
        print "Extracting coords from folder %s" % folder
        writeFilenamesExtractedCoordinates(pathFolder, resname, atomIds, writeLigandTrajectory, constants)
        print "Repeating snapshots from folder %s" % folder
        repeatExtractedSnapshotsInFolder(pathFolder, constants, totalSteps)
        print "Gathering trajs in %s" % constants.gatherTrajsFolder
        gatherTrajs(constants, folder, setNum)


if __name__ == "__main__":
    folder, atomIds, resname, enforceSequential, writeLigandTrajectory, totalSteps, setNum = parseArguments()
    main(folder, atomIds, resname, totalSteps, enforceSequential, writeLigandTrajectory, setNum)
