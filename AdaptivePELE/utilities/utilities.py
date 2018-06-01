from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
from six import reraise as raise_
import os
import sys
import shutil
import numpy as np
import string
import json
from scipy import linalg
import mdtraj as md
try:
    import cPickle as pickle
except ImportError:
    import pickle
from AdaptivePELE.atomset import RMSDCalculator, atomset
from AdaptivePELE.freeEnergies import utils


class UnsatisfiedDependencyException(Exception):
    __module__ = Exception.__module__


def cleanup(tmpFolder):
    """
        Remove folder if exists

        :param tmpFolder: Folder to remove
        :type tmpFolder: str
    """
    if os.path.exists(tmpFolder):
        shutil.rmtree(tmpFolder)


def makeFolder(outputDir):
    """
        Makes folder

        :param outputDir: Folder filename
        :type outputDir: str
    """
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)


def getSnapshots(trajectoryFile, verbose=False, topology=None):
    """
        Gets the snapshots

        :param trajectoryFile: Trajectory filename
        :type trajectoryFile: str
        :param verbose: Add verbose to snapshots
        :type verbose: bool
        :param topology: Topology file object
        :type topology: str

        :returns: iterable -- Snapshots with information
    """
    ext = os.path.splitext(trajectoryFile)[1]
    if ext == ".pdb":
        with open(trajectoryFile, "r") as inputFile:
            inputFileContent = inputFile.read()

        snapshots = inputFileContent.split("ENDMDL")
        if len(snapshots) > 1:
            snapshots = snapshots[:-1]
        if not verbose:
            return snapshots

        remarkInfo = "REMARK 000 File created using PELE++\nREMARK source            : %s\nREMARK original model nr : %d\nREMARK First snapshot is 1, not 0 (as opposed to report)\n%s"
        snapshotsWithInfo = [remarkInfo % (trajectoryFile, i+1, snapshot) for i, snapshot in enumerate(snapshots)]
    elif ext == ".xtc":
        if topology is None:
            raise ValueError("Topology needed for loading xtc files")
        with md.formats.XTCTrajectoryFile(trajectoryFile) as f:
            snapshotsWithInfo, _, _, _ = f.read()
    else:
        raise ValueError("Unrecongnized file extension for %s" % trajectoryFile)
    return snapshotsWithInfo


def getTrajNum(trajFilename):
    """
        Gets the trajectory number

        :param trajFilename: Trajectory filename
        :type trajFilename: str

        :returns: int -- Trajectory number
    """
    return int(trajFilename.split("_")[-1][:-4])


def calculateContactMapEigen(contactMap):
    """
        Calculates eigenvectors and values of an extended contact map

        :param contactMap: Contact map
        :type contactMap: np.array

        :returns: (numpy.ndarray, numpy.ndarray) -- eigenvalues, eigenvectors
    """
    nLig, nCA = contactMap.shape
    extendedCM = np.zeros((nLig+nCA, nLig+nCA))
    extendedCM[nCA:, :nCA] = contactMap
    extendedCM[:nCA, nCA:] = contactMap.T
    assert (extendedCM == extendedCM.T).all(), "Extended ContactMap not symmetric"
    eiv, eic = np.linalg.eigh(extendedCM)
    return eiv, eic


def assertSymmetriesDict(symmetries, PDB):
    """
        Asserts the symmetry list in a PDB

        :param symmetries: List of symmetry groups
        :type symmetries: list of str
        :param PDB: PDB object to check symmetry list against
        :type PDB: :py:class:`.PDB`

        :raise AssertionError: If an atom is not found in the structure
    """
    for group in symmetries:
        for key in group:
            assert key in PDB.atoms, "Symmetry atom %s not found in initial structure" % key
    if symmetries:
        print("Symmetry dictionary correctly defined!")


def getRMSD(traj, nativePDB, resname, symmetries, topology=None):
    """
        Computes the RMSD of a trajectory, given a native and symmetries

        :param traj: Trajecotry filename
        :type traj: str
        :param nativePDB:  Native PDB object
        :type native PDB: :py:class:`.PDB`
        :param resname: Resname to compute its RMSD
        :type resname: str
        :param symmetries: Symmetries dictionary list with independent symmetry groups
        :type symmetries: list of dict
        :param topology: Topology file for non-pdb trajectories
        :type topology: str

        :return: np.array -- Array with the rmsd values of the trajectory
    """

    snapshots = getSnapshots(traj, topology=topology)
    topology_contents = getTopologyFile(topology)
    rmsds = np.zeros(len(snapshots))
    RMSDCalc = RMSDCalculator.RMSDCalculator(symmetries)
    for i, snapshot in enumerate(snapshots):
        snapshotPDB = atomset.PDB()
        snapshotPDB.initialise(snapshot, resname=resname, topology=topology_contents)

        rmsds[i] = RMSDCalc.computeRMSD(nativePDB, snapshotPDB)

    return rmsds


def readClusteringObject(clusteringObjectPath):
    """
        Reads and returns a clustering object

        :param clusteringObjectPath: Clustering object path
        :type clusteringObjectPath: str

        :raise EOFError: If the object is empty

        :returns: :py:class:`.Clustering` -- clusteringObject
    """
    with open(clusteringObjectPath, 'rb') as f:
        try:
            return pickle.load(f)
        except EOFError:
            t, v, tb = sys.exc_info()
            raise_(t, v, tb)


def ensure_connectivity_msm(msm):
    if msm.nstates_full == msm.nstates:
        return msm.stationary_distribution
    else:
        counts = msm.count_matrix_full
        counts += 1/counts.shape[0]
        trans = utils.buildRevTransitionMatrix(counts)
        _, eic = getSortedEigen(trans)
        return getStationaryDistr(eic[:, 0])


def getStationaryDistr(lowestEigenvector):
    absStationary = np.abs(lowestEigenvector)
    return absStationary / absStationary.sum()


def getSortedEigen(T):
    eigenvals, eigenvectors = linalg.eig(T, left=True, right=False)
    sortedIndices = np.argsort(eigenvals)[::-1]

    # reigenvals, reigenvectors = linalg.eig(T.T)
    # rsortedIndices = np.argsort(reigenvals)[::-1]
    return eigenvals[sortedIndices], eigenvectors[:, sortedIndices]


def get_epoch_folders(path):
    """
        List the folders belonging to an adaptive simulation and containing
        trajectories and reports

        :param path: Path where to check for the folders
        :type path: str

        :returns: list -- List of folders belonging to the simulation, sorted

    """
    allFolders = os.listdir(path)
    folders = [epoch for epoch in allFolders if epoch.isdigit()]
    folders.sort(key=int)
    return folders


def gen_atom_name(index):
    # 6760 = 26*26*10
    # The maximum number of unique ids generated by this method is 175760 =
    # 26*26*26*10, however, remember that generally the maximum number of unique
    # atoms that can be used with vmd is around 32000
    ind1 = index//6760
    ind2 = (index % 6760)
    ind3 = ind2 % 260
    return chr(65+ind1)+chr(65+ind2//260)+chr(65+ind3//10)+str(ind3 % 10)


def write_PDB_clusters(pmf_xyzg, title="clusters.pdb", use_beta=False, elements=None):
    templateLine = "HETATM%s %s CLT L 502    %s%s%s  0.75%s          %s  \n"
    if elements is None:
        elements = ["H" for i in range(len(pmf_xyzg))]

    content = ""
    names = []
    for i, line in enumerate(pmf_xyzg):
        number = str(i).rjust(5)
        # number3 = str(i).ljust(3)
        number3 = gen_atom_name(i).ljust(4)
        names.append(number3)
        x = ("%.3f" % line[0]).rjust(8)
        y = ("%.3f" % line[1]).rjust(8)
        z = ("%.3f" % line[2]).rjust(8)
        element = elements[i].rjust(2)
        if use_beta:
            g = ("%.2f" % line[-1]).rjust(6)
        else:
            g = ("%.2f" % 0).rjust(6)
        content += templateLine % (number, number3, x, y, z, g, element)

    with open(title, 'w') as f:
        f.write(content)
    return names


def distanceCOM(coords1, coords2):
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    return np.linalg.norm(coords1-coords2)


def sign(x, tol=1e-7):
    """
        Return the sign of a number

        :param x: Array of number to evaluate the sign
        :type x: numpy.array
        :param tol: Tolerance to define the zero
        :type tol: float

        :returns: int -- Sign of the number
    """
    x[abs(x) < tol] = 0
    return np.sign(x)


def getAtomNames(values):
    """
        Assign to each value an atomname, O for negatives, H for 0 and N for
        positives. This function is created for assign atomnames for custom pdbs

        :param values: Collection of numbers to assing an atom name for
        :type values: iterable

        :returns: list -- List of atom names
    """
    names = ["O", "H", "N"]
    values += 1
    return [names[int(value)] for value in values]


def getReportAndTrajectoryWildcard(JSONdict):
    """
        Extract the trajectory and report filename from the pele control file

        :param JSONdict: Dictionary containing a parsed PELE control file
        :type JSONdict: dict

        :returns: str, str -- Report and trajectory wildcards
    """
    reportWildcard = os.path.split(JSONdict["commands"][0]["PELE_Output"]['reportPath'])[1]
    trajWildcard = os.path.split(JSONdict["commands"][0]["PELE_Output"]['trajectoryPath'])[1]
    trajWildcard = '_%d'.join(os.path.splitext(trajWildcard))
    reportWildcard = '_%d'.join(os.path.splitext(reportWildcard))
    return reportWildcard, trajWildcard


def getPELEControlFileDict(templetizedControlFile):
    """
        Parse a PELE control file into a python dictionary
    """

    with open(templetizedControlFile) as fc:
        peleControlFile = fc.read()

    templateNames = {ele[1]: '"$%s"' % ele[1] for ele in string.Template.pattern.findall(peleControlFile)}
    templateNames.pop("OUTPUT_PATH", None)
    return json.loads(string.Template(peleControlFile).safe_substitute(templateNames)), templateNames


def getMetricsFromReportsInEpoch(reportName, outputFolder, nTrajs):
    """
        Extract the metrics in report file from an epoch to a numpy array
    """
    metrics = []
    for i in range(1, nTrajs):
        report = np.loadtxt(os.path.join(outputFolder, reportName % i))
        if len(report.shape) < 2:
            metrics.append(report.tolist()+[i, 0])
        else:
            traj_line = np.array([i] * report.shape[0])
            snapshot_line = np.array(range(report.shape[0]))
            metrics.extend(np.hstack((report, traj_line[:, np.newaxis], snapshot_line[:, np.newaxis])))
    return np.array(metrics)


def getSASAcolumnFromControlFile(JSONdict):
    for i, metricBlock in enumerate(JSONdict["commands"][0]["PeleTasks"][0]['metrics']):
        if 'sasa' in metricBlock['type'].lower():
            # Report files have 4 fixed columnts, task, steps, accepted steps
            # and energy
            return i+4
    raise ValueError("No SASA metric found in control file!!! Please add it in order to use the moving box feature")


def getTopologyFile(structure):
    """
        Extract the topology information to write structures from xtc format

        :param structure: Pdb file with the topology information
        :type structure: str

        :return: list of str -- The lines of the topology file
    """
    top = []
    with open(structure) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            else:
                top.append("".join([line[:30], "%s%s%s", line[54:]]))
    return top


def write_mdtraj_object_PDB(conformation, output, topology):
    """
        Write a snapshot from a xtc trajectory to pdb

        :param conformation: Mdtraj trajectory object to write
        :type structure: str or Trajectory
        :param output: Output where to write the object
        :type output: str
        :param topology: Topoloy-like object
        :type topology: list
    """
    PDB = atomset.PDB()
    PDB.initialise(conformation, topology=topology)
    PDB.writePDB(output)


def get_mdtraj_object_PDBstring(conformation, topology):
    """
        Get the pdb string of a snapshot from a xtc trajectory to pdb

        :param conformation: Mdtraj trajectory object to write
        :type structure: str or Trajectory
        :param topology: Topoloy-like object
        :type topology: list

        :returns: str -- The pdb representation of a snapshot from a xtc
    """
    PDB = atomset.PDB()
    PDB.initialise(conformation, topology=topology)
    return PDB.get_pdb_string()


def write_xtc_to_pdb(filename, output_file, topology):
    """
        Get the pdb string of a snapshot from a xtc trajectory to pdb

        :param filename: Path to the xtc trajectory
        :type filename: str
        :param output_file: Output where to write the object
        :type output_file: str
        :param topology: Topology file object
        :type topology: str

        :returns: str -- The pdb representation of a snapshot from a xtc
    """
    topology_contents = getTopologyFile(topology)
    xtc_object = md.load(filename, top=topology)
    pdb = atomset.PDB()
    pdb.initialise(xtc_object, topology=topology_contents)
    pdb.writePDB(output_file)


def convert_trajectory_to_pdb(trajectory, topology, output, output_folder):
    """
        Write a trajectory from a non-pdb trajectory to pdb format

        :param trajectory: Trajectory to convert
        :type trajectory: str
        :param topology: Topology file
        :type topology: str
        :param output: Filename of the ouput file
        :type output: str
        :param output_folder: Folder where to store the output trajectory
        :type output_folder: str
    """
    output = os.path.join(output_folder, output)
    topology_contents = getTopologyFile(topology)
    traj = md.load(trajectory, top=topology)
    with open(output, "w") as fw:
        for i in range(traj.n_frames):
            conf = traj.slice(i, copy=False)
            PDB = atomset.PDB()
            PDB.initialise(conf, topology=topology_contents)
            fw.write("MODEL %d\n" % (i+1))
            fw.write(PDB.pdb)
            fw.write("ENDMDL\n")
        fw.write("END\n")


def writeObject(filename, object_to_write):
    with open(filename, "wb") as f:
        pickle.dump(object_to_write, f)
