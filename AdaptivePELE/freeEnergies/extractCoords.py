# coding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import os
import argparse
import glob
import re
import socket
import shutil
import sys
import numpy as np
import multiprocess as mp
from AdaptivePELE.atomset import atomset
from AdaptivePELE.freeEnergies import utils
# reload(sys)
# sys.setdefaultencoding('utf-8')


class Constants:
    def __init__(self):
        self.extractedTrajectoryFolder = "%s/extractedCoordinates"
        self.baseExtractedTrajectoryName = "coord_"
        self.reportName = '*report_'
        self.outputTrajectoryFolder = "%s/repeatedExtractedCoordinates"
        self.ligandTrajectoryFolder = "ligand_trajs"
        self.ligandTrajectoryBasename = os.path.join(self.ligandTrajectoryFolder, "traj_ligand_%s.pdb")
        self.gatherTrajsFolder = "allTrajs"
        self.gatherTrajsFilename = os.path.join(self.gatherTrajsFolder, "traj_%s_%s.dat")
        self.gatherNonRepeatedFolder = os.path.join(self.gatherTrajsFolder, "extractedCoordinates")
        self.gatherNonRepeatedTrajsFilename = os.path.join(self.gatherNonRepeatedFolder, "traj_%s_%s.dat")


def parseArguments():
    desc = "Program that extracts residue coordinates for a posterior MSM analysis.\
            It either extracts the resname COM coordinates or those of an atomId, depending on the input.\
            It then fills the rejected steps, which is not done by PELE.\
            Finally, trajectories are gathered together in the same allTrajs folder.\
            It automatically detects whether it is an adaptive or a sequential PELE run by looking for folders\
            with numeric names."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-f", "--folderWithTrajs", default=".",
                        help="Folder with trajectories (or epochs)")
    # parser.add_argument("-atomId", action=AtomIdAction, help="serial:atomName:resname, e.g. 2048:C1:AIN")
    parser.add_argument("-atomIds", nargs='*', help="serial:atomName:resname, e.g. 2048:C1:AIN. May contain more than one atomId")
    parser.add_argument("-resname", default="", help="Ligand resname")
    parser.add_argument("-CA", "--proteinCA", action="store_true", help="Extract protein alpha carbons coordinates")
    parser.add_argument("-s", "--enforceSequential", action="store_true", help="Force the consideration as sequential run (non-adaptive)")
    parser.add_argument("--setNum", type=int, default=0, help="Sets the number to appear in gathered trajectory in order to avoid clashes between different sequential runs. Ignored in adaptive runs.")
    parser.add_argument("-w", "--writeLigandTrajectory", action="store_true", help="It writes a traj_ligand_XXX.pdb file with the ligand coordinates. The user must delete the original trajectory (if wanted)")
    parser.add_argument("-t", "--totalSteps", type=int, default=0, help="Total number of steps in traj. Equivalent to epoch length in adaptive runs")
    parser.add_argument("-nR", "--noRepeat", action="store_true", help="Flag to avoid repeating the rejected steps")
    parser.add_argument("-n", "--numProcessors", type=int, default=None, help="Number of cpus to use")
    # parser.add_argument("-f", nargs='+', help="Files to get coordinates")
    args = parser.parse_args()

    return args.folderWithTrajs, args.atomIds, args.resname, args.proteinCA, args.enforceSequential, args.writeLigandTrajectory, args.totalSteps, args.setNum, args.noRepeat, args.numProcessors


def getCpuCount():
    machine = socket.getfqdn()
    cores = None
    if "bsccv" in machine:
        # life cluster
        cores = os.getenv("SLURM_NTASKS", None)
    elif "mn.bsc" in machine:
        # nord3
        cores = os.getenv("LSB_DJOB_NUMPROC", None)
    elif "bsc.mn" in machine:
        # MNIV
        cores = os.getenv("SLURM_NPROC", None)
    return cores or mp.cpu_count()


def loadAllResnameAtomsInPdb(filename, lig_resname, writeCA):
    prunedFileContent = []
    with open(filename) as f:
        prunedSnapshot = []
        for line in f:
            if utils.is_model(line):
                prunedFileContent.append("".join(prunedSnapshot))
                prunedSnapshot = []
            elif line[17:20] == lig_resname or utils.isAlphaCarbon(line, writeCA):
                prunedSnapshot.append(line)
        if prunedSnapshot:
            prunedFileContent.append("".join(prunedSnapshot))
    return prunedFileContent


def extractFilenumber(filename):
    last = filename.rfind('.')
    first = filename.rfind('_')
    number = re.sub("[^0-9]", "", filename[first+1:last])
    return number


def getOutputFilename(directory, filename, baseOutputFilename):
    filenumber = extractFilenumber(filename)
    return os.path.join(directory, baseOutputFilename+filenumber+".dat")


def getLigandAlphaCarbonsCoords(allCoordinates, lig_resname):
    trajCoords = []
    for coordinates in allCoordinates:
        PDB = atomset.PDB()
        PDB.initialise(coordinates, resname=lig_resname)
        snapshotCoords = [coord for at in PDB.atomList for coord in PDB.atoms[at].getAtomCoords()]
        PDBCA = atomset.PDB()
        PDBCA.initialise(coordinates, type="PROTEIN")
        snapshotCoords.extend([coord for at in PDBCA.atomList for coord in PDBCA.atoms[at].getAtomCoords()])
        trajCoords.append(snapshotCoords)
    return trajCoords


def getPDBCOM(allCoordinates, lig_resname):
    COMs = []
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=lig_resname, heavyAtoms=True)
        COMs.append(pdb.extractCOM())
    return COMs


def getAtomCoord(allCoordinates, lig_resname, atom_Ids):
    coords = []
    # If ever need to speed this up, build a Trajectory class that inherits from PDB
    # and loads the atom according to the position in the snapshot, rather than looking
    # for the atom
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=lig_resname, heavyAtoms=True)
        snapshotcoords = []
        for atomId in atom_Ids:
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


def writeFilenameExtractedCoordinates(filename, lig_resname, atom_Ids, pathFolder, writeLigandTrajectory, constants, writeCA):
    allCoordinates = loadAllResnameAtomsInPdb(filename, lig_resname, writeCA)
    if writeLigandTrajectory:
        outputFilename = os.path.join(pathFolder, constants.ligandTrajectoryBasename % extractFilenumber(filename))
        with open(outputFilename, 'w') as f:
            f.write("\nENDMDL\n".join(allCoordinates))
    if writeCA:
        coords = getLigandAlphaCarbonsCoords(allCoordinates[:-1], lig_resname)
    else:
        if atom_Ids is None or len(atom_Ids) == 0:
            coords = getPDBCOM(allCoordinates, lig_resname)
        else:
            coords = getAtomCoord(allCoordinates, lig_resname, atom_Ids)

    outputFilename = getOutputFilename(constants.extractedTrajectoryFolder, filename,
                                       constants.baseExtractedTrajectoryName)
    writeToFile(coords, outputFilename % pathFolder)


def writeFilenamesExtractedCoordinates(pathFolder, lig_resname, atom_Ids, writeLigandTrajectory, constants, writeCA, pool=None):
    if not os.path.exists(constants.extractedTrajectoryFolder % pathFolder):
        os.makedirs(constants.extractedTrajectoryFolder % pathFolder)

    originalPDBfiles = glob.glob(pathFolder+'/*traj*.pdb')
    workers = []
    for filename in originalPDBfiles:
        if pool is None:
            # serial version
            writeFilenameExtractedCoordinates(filename, lig_resname, atom_Ids, pathFolder, writeLigandTrajectory, constants, writeCA)
        else:
            # multiprocessor version
            workers.append(pool.apply_async(writeFilenameExtractedCoordinates, args=(filename, lig_resname, atom_Ids, pathFolder, writeLigandTrajectory, constants, writeCA)))
    for w in workers:
        w.get()


def parseResname(atom_Ids, lig_resname):
    if atom_Ids is not None and len(atom_Ids) > 0:
        differentResnames = {atomId.split(":")[-1] for atomId in atom_Ids}
        if len(differentResnames) > 1:
            sys.exit("Error! Different resnames provided in atomIds!")
        elif len(differentResnames) == 1:
            extractedResname = differentResnames.pop()

    if (atom_Ids is None or len(atom_Ids) == 0) and lig_resname == "":
        sys.exit("Either resname or atomId should be provided")
    elif lig_resname == "":
        lig_resname = extractedResname  # the atom Id last element is the resname
    elif atom_Ids is not None and len(atom_Ids) > 0:
        if extractedResname != lig_resname:
            sys.exit("Residue name in resname and atomId do not match!")
    return lig_resname


def buildFullTrajectory(steps, trajectory, numtotalSteps, inputTrajectory):
    completeTrajectory = []
    counter = 0
    if len(trajectory) > 0:
        sthWrongInTraj = False
        for i in range(len(trajectory) - 1):
            try:
                repeated = steps[i+1] - steps[i]
            except IndexError:
                print("sth wrong in trajectory %s. This is likely to disagreement between report and trajectory files. Please, fix it manually" % inputTrajectory)
                sthWrongInTraj = True
                break

            for _ in range(repeated):
                snapshot = trajectory[i].split()
                snapshot[0] = str(counter)
                snapshot = ' '.join(snapshot)
                completeTrajectory.append(snapshot)
                counter += 1

        if sthWrongInTraj:
            return completeTrajectory

        if numtotalSteps == 0:
            iterations = list(range(1))
        else:
            iterations = list(range(numtotalSteps + 1 - counter))

        for i in iterations:
            snapshot = trajectory[-1].split()
            snapshot[0] = str(counter)
            snapshot = ' '.join(snapshot)
            completeTrajectory.append(snapshot)
            counter += 1

    return completeTrajectory


def repeatExtractedSnapshotsInTrajectory(inputTrajectory, constants, numtotalSteps):
    extractedTrajFolder, trajFilename = os.path.split(inputTrajectory)
    trajectoryNumber = re.sub(r'\.dat$', '', trajFilename)
    trajectoryNumber = re.sub(constants.baseExtractedTrajectoryName, '', trajectoryNumber)

    origDataFolder = re.sub(constants.extractedTrajectoryFolder % "", "", extractedTrajFolder)
    try:
        reportFile = glob.glob(os.path.join(origDataFolder, constants.reportName + trajectoryNumber))[0]
    except IndexError:
        print("folder", origDataFolder)
        sys.exit("Couldn't find file that matches: %s" % os.path.join(origDataFolder, constants.reportName + trajectoryNumber))

    with open(inputTrajectory) as f:
        trajectory = f.read().splitlines()

    acceptedSteps = np.loadtxt(reportFile, dtype='int', comments='#', usecols=(1,))

    fullTrajectory = buildFullTrajectory(acceptedSteps, trajectory, numtotalSteps, inputTrajectory)

    if len(fullTrajectory) > 0:
        outputFilename = os.path.join(constants.outputTrajectoryFolder % origDataFolder, constants.baseExtractedTrajectoryName + trajectoryNumber + '.dat')
        outputFile = open(outputFilename, 'w')
        for snapshot in fullTrajectory:
            outputFile.write("%s\n" % snapshot)
        outputFile.close()


def repeatExtractedSnapshotsInFolder(folder_name, constants, numtotalSteps, pool=None):
    inputTrajectoryFolder = constants.extractedTrajectoryFolder % folder_name
    outputTrajectoryFolder = constants.outputTrajectoryFolder % folder_name

    if not os.path.exists(outputTrajectoryFolder):
        os.makedirs(outputTrajectoryFolder)

    inputTrajectories = glob.glob(os.path.join(inputTrajectoryFolder, constants.baseExtractedTrajectoryName + '*'))
    workers = []
    for inputTrajectory in inputTrajectories:
        if pool is None:
            # serial version
            repeatExtractedSnapshotsInTrajectory(inputTrajectory, constants, numtotalSteps)
        else:
            # multiprocessor version
            workers.append(pool.apply_async(repeatExtractedSnapshotsInTrajectory, args=(inputTrajectory, constants, numtotalSteps)))
    for w in workers:
        w.get()


def makeGatheredTrajsFolder(constants):
    if not os.path.exists(constants.gatherTrajsFolder):
        os.makedirs(constants.gatherTrajsFolder)
    if not os.path.exists(constants.gatherNonRepeatedFolder):
        os.makedirs(constants.gatherNonRepeatedFolder)


def copyTrajectories(traj_names, destFolderTempletized, folderName):
    for inputTrajectory in traj_names:
        trajectoryNumber = extractFilenumber(os.path.split(inputTrajectory)[1])
        if folderName != ".":  # if not sequential
            setNumber = folderName
        shutil.copyfile(inputTrajectory, destFolderTempletized % (setNumber, trajectoryNumber))


def gatherTrajs(constants, folder_name, setNumber, non_Repeat):
    if non_Repeat:
        trajectoriesFilenames = os.path.join(constants.extractedTrajectoryFolder % folder_name, constants.baseExtractedTrajectoryName + "*")
    else:
        trajectoriesFilenames = os.path.join(constants.outputTrajectoryFolder % folder_name, constants.baseExtractedTrajectoryName + "*")
    trajectories = glob.glob(trajectoriesFilenames)
    copyTrajectories(trajectories, constants.gatherTrajsFilename, folder_name)
    nonRepeatedTrajs = glob.glob(os.path.join(constants.extractedTrajectoryFolder % folder_name, constants.baseExtractedTrajectoryName + "*"))
    copyTrajectories(nonRepeatedTrajs, constants.gatherNonRepeatedTrajsFilename, folder_name)


def main(folder_name=".", atom_Ids="", lig_resname="", numtotalSteps=0, enforceSequential_run=0, writeLigandTrajectory=True, setNumber=0, protein_CA=0, non_Repeat=False, nProcessors=None):

    constants = Constants()

    lig_resname = parseResname(atom_Ids, lig_resname)

    folderWithTrajs = folder_name

    makeGatheredTrajsFolder(constants)

    # change atomId for list

    if enforceSequential_run:
        folders = ["."]
    else:
        allFolders = os.listdir(folderWithTrajs)
        folders = [epoch for epoch in allFolders if epoch.isdigit()]
        if len(folders) == 0:
            folders = ["."]

    if nProcessors is None:
        nProcessors = getCpuCount()
    nProcessors = max(1, nProcessors)
    print("Running extractCoords with %d cores" % (nProcessors))
    pool = mp.Pool()

    for folder_it in folders:
        pathFolder = os.path.join(folderWithTrajs, folder_it)
        print("Extracting coords from folder %s" % folder_it)
        ligand_trajs_folder = os.path.join(pathFolder, constants.ligandTrajectoryFolder)
        if writeLigandTrajectory and not os.path.exists(ligand_trajs_folder):
            os.makedirs(ligand_trajs_folder)
        writeFilenamesExtractedCoordinates(pathFolder, lig_resname, atom_Ids, writeLigandTrajectory, constants, protein_CA, pool=pool)
        if not non_Repeat:
            print("Repeating snapshots from folder %s" % folder_it)
            repeatExtractedSnapshotsInFolder(pathFolder, constants, numtotalSteps, pool=None)
        print("Gathering trajs in %s" % constants.gatherTrajsFolder)
        gatherTrajs(constants, folder_it, setNumber, non_Repeat)


if __name__ == "__main__":
    folder, atomIds, resname, proteinCA, enforceSequential, writeLigandTraj, totalSteps, setNum, nonRepeat, n_processors = parseArguments()
    main(folder, atomIds, resname, totalSteps, enforceSequential, writeLigandTraj, setNum, proteinCA, nonRepeat, n_processors)
