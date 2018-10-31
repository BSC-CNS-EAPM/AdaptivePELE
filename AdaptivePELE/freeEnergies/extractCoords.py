# coding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import os
import argparse
import glob
import re
import shutil
import sys
import mdtraj as md
import numpy as np
from AdaptivePELE.atomset import atomset
from AdaptivePELE.freeEnergies import utils
from AdaptivePELE.utilities import utilities
PARALELLIZATION = True
try:
    import multiprocessing as mp
except ImportError:
    PARALELLIZATION = False
PRODY = True
try:
    import prody as pd
except ImportError:
    PRODY = False


MDTRAJ_FORMATS = set(['.xtc', '.dcd', '.dtr', '.trr', 'mdcrd', 'nc'])


class Constants(object):
    def __init__(self):
        self.extractedTrajectoryFolder = "%s/extractedCoordinates"
        self.baseExtractedTrajectoryName = "coord_"
        self.reportName = '*report_'
        self.baseGatheredFilename = "traj_*.dat"
        self.outputTrajectoryFolder = "%s/repeatedExtractedCoordinates"
        self.ligandTrajectoryFolder = "ligand_trajs"
        self.ligandTrajectoryBasename = os.path.join(self.ligandTrajectoryFolder, "traj_ligand_%s.pdb")
        self.gatherTrajsFolder = "allTrajs"
        self.gatherTrajsFilename = os.path.join(self.gatherTrajsFolder, "traj_%s_%s.dat")
        self.gatherNonRepeatedFolder = os.path.join(self.gatherTrajsFolder, "extractedCoordinates")
        self.gatherNonRepeatedTrajsFilename = os.path.join(self.gatherNonRepeatedFolder, "traj_%s_%s.dat")


class TopologyCompat(object):
    def __init__(self, pdb_file):
        self.path = os.path.split(os.path.abspath(pdb_file))[0]
        self.topologyFiles = pdb_file

    def getTopologyFile(self, epoch, trajectory_number):
        return self.topologyFiles


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
    parser.add_argument("--top", type=str, default=None, help="Topology file for non-pdb trajectories or path to Adaptive topology object")
    parser.add_argument("--sidechains", action="store_true", help="Flag to extract sidechain coordinates")
    parser.add_argument("-sf", "--sidechains_folder", default=".", type=str, help="Folder with the structures to obtain the sidechains to extract")
    parser.add_argument("--serial", action="store_true", help="Flag to deactivate parallelization")
    # parser.add_argument("-f", nargs='+', help="Files to get coordinates")
    args = parser.parse_args()

    return args.folderWithTrajs, args.atomIds, args.resname, args.proteinCA, args.enforceSequential, args.writeLigandTrajectory, args.totalSteps, args.setNum, args.noRepeat, args.numProcessors, args.top, args.sidechains, args.sidechains_folder, args.serial


def loadAllResnameAtomsInPdb(filename, lig_resname, writeCA, sidechains):
    prunedFileContent = []
    sidechains_bool = bool(sidechains)
    with open(filename) as f:
        prunedSnapshot = []
        for line in f:
            if utils.is_model(line):
                prunedFileContent.append("".join(prunedSnapshot))
                prunedSnapshot = []
            elif line[17:20] == lig_resname or utils.isAlphaCarbon(line, writeCA) or utils.isSidechain(line, sidechains_bool, sidechains):
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


def getLigandAlphaCarbonsCoords(allCoordinates, lig_resname, sidechains=False):
    trajCoords = []
    for coordinates in allCoordinates:
        PDB = atomset.PDB()
        PDB.initialise(coordinates, resname=lig_resname)
        snapshotCoords = [coord for at in PDB.atomList for coord in PDB.atoms[at].getAtomCoords()]
        PDBCA = atomset.PDB()
        if not sidechains:
            PDBCA.initialise(coordinates, type="PROTEIN")
        else:
            PDBCA.initialise(coordinates, type="PROTEIN", heavyAtoms=True)
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


def extractIndexesTopology(topology, lig_resname, atoms, writeCA, sidechains):
    selection = []
    if atoms:
        atoms_set = set(atoms)
    template = "%s:%s:%s"
    iline = 0
    bool_sidechains = bool(sidechains)
    with open(topology) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if atoms:
                serial_num = line[6:11].strip()
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                if template % (serial_num, atom_name, residue_name) in atoms_set:
                    selection.append(iline)

            elif (line[17:20] == lig_resname or utils.isAlphaCarbon(line, writeCA) or utils.isSidechain(line, bool_sidechains, sidechains)) and line[76:80].strip().upper() != "H":
                selection.append(iline)
            iline += 1
    return selection


def extractCoordinatesXTCFile(file_name, ligand, atom_Ids, writeCA, topology, selected_indices, sidechains):
    trajectory = md.load(file_name, top=topology)
    if not writeCA and (atom_Ids is None or len(atom_Ids) == 0) and not sidechains:
        # getCOM case
        # convert nm to A
        coordinates = 10*md.compute_center_of_mass(trajectory.atom_slice(selected_indices))
    else:
        coordinates = 10*trajectory.xyz[:, selected_indices, :].reshape((trajectory.n_frames, -1))
    return coordinates


def writeFilenameExtractedCoordinates(filename, lig_resname, atom_Ids, pathFolder, writeLigandTrajectory, constants, writeCA, sidechains, topology=None, indexes=None):
    """
        Process the coordinates of a trajectory

    """
    ext = utilities.getFileExtension(filename)
    if ext == ".pdb":
        allCoordinates = loadAllResnameAtomsInPdb(filename, lig_resname, writeCA, sidechains)
        if writeLigandTrajectory:
            outputFilename = os.path.join(pathFolder, constants.ligandTrajectoryBasename % extractFilenumber(filename))
            with open(outputFilename, 'w') as f:
                f.write("\nENDMDL\n".join(allCoordinates))
        if writeCA:
            coords = getLigandAlphaCarbonsCoords(allCoordinates, lig_resname)
        elif sidechains:
            coords = getLigandAlphaCarbonsCoords(allCoordinates, lig_resname, sidechains=sidechains)
        else:
            if atom_Ids is None or len(atom_Ids) == 0:
                coords = getPDBCOM(allCoordinates, lig_resname)
            else:
                coords = getAtomCoord(allCoordinates, lig_resname, atom_Ids)
    elif ext in MDTRAJ_FORMATS:
        coords = extractCoordinatesXTCFile(filename, lig_resname, atom_Ids, writeCA, topology, indexes, sidechains)
    else:
        raise ValueError("Unrecongnized file extension for %s" % filename)

    outputFilename = getOutputFilename(constants.extractedTrajectoryFolder, filename,
                                       constants.baseExtractedTrajectoryName)
    writeToFile(coords, outputFilename % pathFolder)


def writeFilenamesExtractedCoordinates(pathFolder, lig_resname, atom_Ids, writeLigandTrajectory, constants, writeCA, sidechains, pool=None, topology=None):
    if not os.path.exists(constants.extractedTrajectoryFolder % pathFolder):
        os.makedirs(constants.extractedTrajectoryFolder % pathFolder)

    originalPDBfiles = glob.glob(os.path.join(pathFolder, '*traj*.*'))
    ext = os.path.splitext(originalPDBfiles[0])[1]
    if ext in MDTRAJ_FORMATS:
        if topology is None:
            raise ValueError("Necessary topology not provided!")
        # get topology for the first trajectory
        top_file = topology.getTopologyFile(0, 1)
        indexes = extractIndexesTopology(top_file, lig_resname, atom_Ids, writeCA, sidechains)
    else:
        indexes = None
    workers = []
    for filename in originalPDBfiles:
        if topology is not None:
            epoch, traj_num = get_epoch_traj_num(filename)
            topology_file = topology.getTopologyFile(epoch, traj_num)
        else:
            topology_file = None
        if pool is None:
            # serial version
            writeFilenameExtractedCoordinates(filename, lig_resname, atom_Ids, pathFolder, writeLigandTrajectory, constants, writeCA, sidechains, topology=topology_file, indexes=indexes)
        else:
            # multiprocessor version
            workers.append(pool.apply_async(writeFilenameExtractedCoordinates, args=(filename, lig_resname, atom_Ids, pathFolder, writeLigandTrajectory, constants, writeCA, sidechains, topology_file, indexes)))
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
        for i in range(len(steps) - 1):
            repeated = steps[i+1, 0] - steps[i, 0]
            for _ in range(repeated):
                try:
                    snapshot = trajectory[steps[i, 1]].split()
                except IndexError:
                    print("sth wrong in trajectory %s. This is likely to disagreement between report and trajectory files. Please, fix it manually" % inputTrajectory)
                    sthWrongInTraj = True
                    break
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
        sys.exit("Couldn't find file that matches: %s" % os.path.join(origDataFolder, constants.reportName + trajectoryNumber))

    with open(inputTrajectory) as f:
        trajectory = f.read().splitlines()

    acceptedSteps = np.loadtxt(reportFile, dtype='int', comments='#', usecols=(1, 2))
    if len(acceptedSteps.shape) < 2:
        acceptedSteps = acceptedSteps[np.newaxis, :]

    fullTrajectory = buildFullTrajectory(acceptedSteps, trajectory, numtotalSteps, inputTrajectory)

    if len(fullTrajectory) > 0:
        outputFilename = os.path.join(constants.outputTrajectoryFolder % origDataFolder, constants.baseExtractedTrajectoryName + trajectoryNumber + '.dat')
        with open(outputFilename, "w") as outputFile:
            for snapshot in fullTrajectory:
                outputFile.write("%s\n" % snapshot)


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


def copyTrajectories(traj_names, destFolderTempletized, folderName, setNumber=0, epochNum=None):
    for inputTrajectory in traj_names:
        trajectoryNumber = extractFilenumber(os.path.split(inputTrajectory)[1])
        if folderName != ".":  # if not sequential
            setNumber = folderName
        if epochNum is not None:
            setNumber = epochNum
        shutil.copyfile(inputTrajectory, destFolderTempletized % (setNumber, trajectoryNumber))


def gatherTrajs(constants, folder_name, setNumber, non_Repeat, epochNum=None):
    if not non_Repeat:
        trajectoriesFilenames = os.path.join(constants.outputTrajectoryFolder % folder_name, constants.baseExtractedTrajectoryName + "*")
        trajectories = glob.glob(trajectoriesFilenames)
        copyTrajectories(trajectories, constants.gatherTrajsFilename, folder_name, setNumber)
    nonRepeatedTrajs = glob.glob(os.path.join(constants.extractedTrajectoryFolder % folder_name, constants.baseExtractedTrajectoryName + "*"))
    copyTrajectories(nonRepeatedTrajs, constants.gatherNonRepeatedTrajsFilename, folder_name, setNumber, epochNum=epochNum)


def extractSidechainIndexes_prody(traj, ligand_resname, topology=None):
    if not PRODY:
        raise utilities.UnsatisfiedDependencyException("Prody module not found, will not be able to extract sidechain coordinates")
    atoms = pd.parsePDB(traj)
    sidechains = atoms.select("protein within 5 of resname {}".format(ligand_resname))
    return [atom.getIndex() for atom in sidechains]


def extractSidechainIndexes_mdtraj(traj, lig_resname, topology=None):
    atoms = md.load(traj, top=topology)
    ligand_indices = atoms.top.select("resname '{lig}'".format(lig=lig_resname))
    water_indices = set(atoms.top.select("not protein or not resname '{lig}'".format(lig=lig_resname)))
    # the distance is specified in nm
    sidechains = md.compute_neighbors(atoms, 0.5, ligand_indices)
    sidechains_trajs = []
    for _, sidechain in enumerate(sidechains):
        sidechains_trajs.extend(list(set(sidechain.tolist())-water_indices))
    return sidechains_trajs


def extractSidechainIndexes(trajs, lig_resname, topology=None, pool=None):
    sidechains_trajs = []
    workers = []
    for traj in trajs:
        ext = utilities.getFileExtension(traj)

        if ext == ".pdb":
            if PRODY:
                if pool is None:
                    sidechains_trajs.extend(extractSidechainIndexes_prody(traj, lig_resname))
                else:
                    workers.append(pool.apply_async(extractSidechainIndexes_prody, args=(traj, lig_resname)))
            else:
                if pool is None:
                    sidechains_trajs.extend(extractSidechainIndexes_mdtraj(traj, lig_resname))
                else:
                    workers.append(pool.apply_async(extractSidechainIndexes_mdtraj, args=(traj, lig_resname)))
        elif ext in MDTRAJ_FORMATS:
            epoch, traj_num = get_epoch_traj_num(traj)
            if pool is None:
                sidechains_trajs.extend(extractSidechainIndexes_mdtraj(traj, lig_resname, topology=topology.getTopologyFile(epoch, traj_num)))
            else:
                workers.append(pool.apply_async(extractSidechainIndexes_mdtraj(traj, lig_resname, topology)))
        else:
            raise ValueError("Unrecongnized file extension for %s" % traj)
    for w in workers:
        sidechains_trajs.extend(w.get())
    return list(set(sidechains_trajs))


def get_epoch_traj_num(filename):
    # assumes trajectories come from an Adaptive simulation
    path, traj_name = os.path.split(filename)
    try:
        epoch = int(os.path.split(path)[-1])
    except ValueError:
        # if for some reason epoch number can't be inferred, assume first
        # epoch
        epoch = 0
    try:
        traj_num = utilities.getTrajNum(traj_name)
    except ValueError:
        # if for some reason trajectory number can't be inferred, assume
        # first trajectory
        traj_num = 1
    return epoch, traj_num


def getTopologyObject(topology_file):
    ext = utilities.getFileExtension(topology_file)
    if ext == ".pdb":
        return TopologyCompat(topology_file)
    elif ext == ".pkl":
        return utilities.readClusteringObject(topology_file)
    else:
        raise ValueError("The topology parameter needs to be the path to a pickled Topology object or a pdb!")


def main(folder_name=".", atom_Ids="", lig_resname="", numtotalSteps=0, enforceSequential_run=0, writeLigandTrajectory=True, setNumber=0, protein_CA=0, non_Repeat=False, nProcessors=None, parallelize=True, topology=None, sidechains=False, sidechain_folder="."):

    constants = Constants()

    if topology is not None:
        topology = getTopologyObject(topology)

    lig_resname = parseResname(atom_Ids, lig_resname)

    folderWithTrajs = folder_name

    makeGatheredTrajsFolder(constants)

    if enforceSequential_run:
        folders = ["."]
    else:
        allFolders = os.listdir(folderWithTrajs)
        folders = [epoch for epoch in allFolders if epoch.isdigit()]
        if len(folders) == 0:
            folders = ["."]

    # if multiprocess is not available, turn off parallelization
    parallelize &= PARALELLIZATION

    if parallelize:
        if nProcessors is None:
            nProcessors = utilities.getCpuCount()
        nProcessors = max(1, nProcessors)

        print("Running extractCoords with %d cores" % (nProcessors))
        pool = mp.Pool(nProcessors)
    else:
        pool = None

    sidechains = extractSidechainIndexes(glob.glob(sidechain_folder), lig_resname, topology=topology, pool=pool) if sidechains else []

    for folder_it in folders:
        pathFolder = os.path.join(folderWithTrajs, folder_it)
        print("Extracting coords from folder %s" % folder_it)
        ligand_trajs_folder = os.path.join(pathFolder, constants.ligandTrajectoryFolder)
        if writeLigandTrajectory and not os.path.exists(ligand_trajs_folder):
            os.makedirs(ligand_trajs_folder)
        writeFilenamesExtractedCoordinates(pathFolder, lig_resname, atom_Ids, writeLigandTrajectory, constants, protein_CA, sidechains, pool=pool, topology=topology)
        if not non_Repeat:
            print("Repeating snapshots from folder %s" % folder_it)
            repeatExtractedSnapshotsInFolder(pathFolder, constants, numtotalSteps, pool=None)
        print("Gathering trajs in %s" % constants.gatherTrajsFolder)
        gatherTrajs(constants, folder_it, setNumber, non_Repeat)


if __name__ == "__main__":
    folder, atomIds, resname, proteinCA, enforceSequential, writeLigandTraj, totalSteps, setNum, nonRepeat, n_processors, top, side_chains, sideChain_folder, serial = parseArguments()
    main(folder, atomIds, resname, totalSteps, enforceSequential, writeLigandTraj, setNum, proteinCA, nonRepeat, n_processors, topology=top, sidechains=side_chains, sidechain_folder=sideChain_folder, parallelize=(not serial))
