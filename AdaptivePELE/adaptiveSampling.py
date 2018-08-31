from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import shutil
import os
import json
import time
import glob
import argparse
import atexit
import ast
import signal
import numpy as np
from builtins import range
from AdaptivePELE.constants import blockNames, constants
from AdaptivePELE.atomset import atomset
from AdaptivePELE.utilities import utilities
from AdaptivePELE.utilities.synchronization import ProcessesManager
from AdaptivePELE.validator import controlFileValidator
from AdaptivePELE.spawning import spawning, spawningTypes
from AdaptivePELE.simulation import simulationrunner, simulationTypes
from AdaptivePELE.clustering import clustering, clusteringTypes


def parseArgs():
    parser = argparse.ArgumentParser(description="Perform several iterations of"
                                     " simulations using adaptive sampling to "
                                     "distribute the processors in order "
                                     "to optimize sampling")
    parser.add_argument('controlFile', type=str)
    arg = parser.parse_args()
    return arg


class InitialStructuresError(Exception):
    __module__ = Exception.__module__


def cleanPreviousSimulation(output_path):
    """
        Clean the uneeded data from a previous simulation

        :param output_path: Path where the data is stored
        :type output_path: str
    """
    equilibration_folders = glob.glob(os.path.join(output_path, "equilibration*"))
    for folder in equilibration_folders:
        shutil.rmtree(folder)
    epochs = utilities.get_epoch_folders(output_path)
    for epoch in epochs:
        shutil.rmtree(epoch)


def createMappingForFirstEpoch(initialStructures, topologies, processors):
    """
        Create the topology mapping for the first iteration

        :param initialStructures: List of the initial structures for the first iteration
        :type initialStructures: list
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`
        :param processors: Number of trajectories
        :type processors: int

    """
    topologyMapping = list(range(1, len(initialStructures)))+[0]
    topologyMapping = topologyMapping*int(np.ceil(processors/len(initialStructures)))
    topologies.topologyMap[0] = topologyMapping[:processors]


def writeTopologyFiles(topologies, destination):
    """
        Write the topology files to the desired destination

        :param topologies: List of topology files
        :type topologies: list
        :param destination: Path where to copy the toplogy files
        :type destination: str
    """
    destination = os.path.join(destination, "topology_%d.pdb")
    for i, topology in enumerate(topologies):
        shutil.copy(topology, destination % i)


def checkMetricExitConditionMultipleTrajsinRestart(firstRun, outputFolder, simulationRunner):
    """
        Check the previous simulation data when restarting a simulation with a
        multiple trajectory metric exit condition

        :param firstRun: First epoch to be run in restart simulation
        :type firstRun: int
        :param outputFolder: Folder of the previous simulation data
        :type outputFolder: str
        :param simulationRunner: Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
    """
    if simulationRunner.hasExitCondition():
        if simulationRunner.parameters.exitCondition.type != blockNames.ExitConditionType.metricMultipleTrajs:
            return
        for i in range(firstRun):
            simulationRunner.parameters.exitCondition.checkExitCondition(outputFolder % i)


def filterClustersAccordingToBox(simulationRunnerParams, clusteringObject):
    """
        Filter the clusters to select only the ones whose representative
        structures will fit into the selected box

        :param simulationRunnerParams: :py:class:`.SimulationParameters` Simulation parameters object
        :type simulationRunnerParams: :py:class:`.SimulationParameters`
        :param clusteringObject: Clustering object
        :type clusteringObject: :py:class:`.Clustering`

        :returns list, list: -- list of the filtered clusters, list of bools flagging wether the cluster is selected or not

    """
    box_center = ast.literal_eval(simulationRunnerParams.boxCenter)
    box_radius = simulationRunnerParams.boxRadius
    clustersFiltered = []
    clustersSelected = []
    for cluster in clusteringObject.clusters.clusters:
        if utilities.distanceCOM(box_center, cluster.pdb.getCOM()) < (box_radius-1):
            clustersFiltered.append(cluster)
            clustersSelected.append(True)
        else:
            clustersSelected.append(False)
    return clustersFiltered, clustersSelected


def mergeFilteredClustersAccordingToBox(degeneracy, clustersFiltering):
    """
        Merge the (possibly) partial degeneracy to obtain a complete list,
        the partial list comes from the clusters excluded by the moving box
        :param degeneracy: Degeneracy of the clusters
        :param degeneracy: list
        :param clustersFiltering: List of bools indicating whether a cluster was filtered or not
        :param clustersFiltering: list

        :returns list: -- complete list of cluster degeneracies
    """
    assert len(degeneracy) == sum(clustersFiltering)
    newDegeneracy = []
    for filtered in clustersFiltering:
        if filtered:
            newDegeneracy.append(degeneracy.pop(0))
        else:
            newDegeneracy.append(0)
    return np.array(newDegeneracy)


def expandInitialStructuresWildcard(initialStructuresWildcard):
    """
        Returns the initial structures after expanding the initial structures wildcard

        :param initialStructureWildcard: Wildcard that matches the initial structures
        :type initialStructureWildcard: str

        :return: list of str -- The expanded initial structures
    """
    totalInitialStructures = []
    for initialStructureWildcard in initialStructuresWildcard:
        expandedStructures = glob.glob(initialStructureWildcard)
        totalInitialStructures.extend(map(os.path.abspath, expandedStructures))
    return totalInitialStructures


def checkSymmetryDict(clusteringBlock, initialStructures, resname):
    """
        Check if the symmetries dictionary is valid for the ligand

        :param clusteringBlock: JSON block with the clustering-related parameters
        :type clusteringBlock: json
        :param initialStructures: List with initial structures
        :type initialStructures: list
        :param resname: Residue name of the ligand in the system pdb
        :type resname: str

        :raise AssertionError: If atoms are not found in the structure
     """
    symmetries = clusteringBlock[blockNames.ClusteringTypes.params].get(blockNames.ClusteringTypes.symmetries, {})
    for structure in initialStructures:
        PDB = atomset.PDB()
        PDB.initialise(structure, resname=resname)
        utilities.assertSymmetriesDict(symmetries, PDB)


def fixReportsSymmetry(outputPath, resname, nativeStructure, symmetries, topologies):
    """
        Adds a new column in the report file with the RMSD that takes into account symmetries.
        New reports are stored in the fixedReport_i where i is the number of the report

        :param outputPath: Path where trajectories are found
        :type outputPath: str
        :param resname: Residue name of the ligand in the pdb
        :type resname: str
        :param nativeStructure: Path to the native structure pdb
        :type nativeStructure: str
        :param symmetries: Dictionary containg the symmetries of the ligand
        :type symmetries: dict
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`

        :raise IndexError: If original report file not found in output folder
    """
    outputFilename = "fixedReport_%d"  # move to constants?
    trajName = "*traj*"  # move to constants?
    reportName = "*report_%d"  # move to constants?
    epoch = int(os.path.basename(outputPath))
    trajs = glob.glob(os.path.join(outputPath, trajName))
    nativePDB = atomset.PDB()
    nativePDB.initialise(str(nativeStructure), resname=resname)
    for traj in trajs:
        trajNum = utilities.getTrajNum(traj)
        rmsd = list(utilities.getRMSD(traj, nativePDB, resname, symmetries, topology=topologies.getTopology(epoch, trajNum)))
        try:
            reportFilename = glob.glob(os.path.join(outputPath, reportName) % trajNum)[0]
        except IndexError:
            raise IndexError("File %s not found in folder %s" % (reportName % trajNum, outputPath))
        rmsd.insert(0, "\tCorrected RMSD")
        with open(reportFilename, "r") as f:
            report = f.readlines()
        outfile = open(os.path.join(outputPath, outputFilename % trajNum), "w")
        for line, value in zip(report, rmsd):
            outfile.write(line.rstrip("\n")+str(value)+"\n")
        outfile.close()


def copyInitialStructures(initialStructures, tmpInitialStructuresTemplate, iteration):
    """
        Copies the initial structures from a certain iteration

        :param initialStructures: Name of the initial structures to copy
        :type initialStructures: list of str
        :param tmpInitialStructuresTemplate: Template with the name of the initial strutctures
        :type tmpInitialStructuresTemplate: str
        :param iteration: Epoch number
        :type iteration: int
    """

    for i, name in enumerate(initialStructures):
        shutil.copyfile(name, tmpInitialStructuresTemplate % (iteration, i))


def generateTrajectorySelectionString(epoch, epochOutputPathTempletized):
    """
        Generates the template for the name of the trajectories in a given epoch

        :param epoch: Epoch number
        :type epoch: int
        :param epochOutputPathTempletized: Templetized path where the trajectories of any epoch are stored
        :type epochOutputPathTempletized: str

        :returns: str -- Template for the name of the trajectories in a given
            epoch
    """
    return "[\"" + os.path.join(epochOutputPathTempletized % epoch, constants.trajectoryBasename) + "\"]"


def findFirstRun(outputPath, clusteringOutputObject, simulationRunner):
    """
        Find the last epoch that was properly simulated and clusterized and
        and return the first epoch to run in case of restart

        :param outputPath: Simulation output path
        :type outputPath: str
        :param clusteringOutputObject: Templetized name of the clustering object
        :type clusteringOutputObject: str
        :param simulationRunner: Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner

        :return: int -- Current epoch
    """

    folderWithSimulationData = outputPath
    allFolders = os.listdir(folderWithSimulationData)
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort(reverse=True)

    objectsFound = []
    for epoch in epochFolders:
        if simulationRunner.checkSimulationInterrupted(epoch, outputPath):
            # this should only happen in MD simulations, where checkpoints are
            # periodically written in case the adaptive run dies at
            # mid-simulation, be able to use the already computed trajectories
            return epoch
        if os.path.exists(clusteringOutputObject % epoch):
            objectsFound.append(epoch)
        if objectsFound and epoch < (objectsFound[0]-5):
            break
    while objectsFound:
        epoch = objectsFound.pop(0)
        if checkIntegrityClusteringObject(clusteringOutputObject % epoch):
            return epoch + 1
    return None


def checkIntegrityClusteringObject(objectPath):
    """
        Test whether the found clustering object to reload is a valid object

        :param objectPath: Clustering object path
        :type objectPath: str

        :returns: bool -- True if the found clustering object is valid
    """
    try:
        utilities.readClusteringObject(objectPath)
        return True
    except EOFError:
        return False


def __unicodeToStr(data):
    # convert dict
    if isinstance(data, dict):
        return {__unicodeToStr(key): __unicodeToStr(value) for key, value in data.items()}
    # convert list
    if isinstance(data, list):
        return [__unicodeToStr(val) for val in data]
    # convert unicode to str
    if isinstance(data, unicode):
        return data.encode('utf-8')

    return data


def loadParams(jsonParams):
    """
        Read the control file in JSON format and extract the blocks of simulation,
        general parameters, clustering and spawning

        :param jsonParams: Control file in JSON format from where the parameters will be read
        :type jsonParams: json str
    """
    jsonFile = open(jsonParams, 'r').read()
    # parsedJSON = json.loads(jsonFile, object_hook=__unicodeToStr)
    parsedJSON = json.loads(jsonFile)

    return parsedJSON[blockNames.ControlFileParams.generalParams], parsedJSON[blockNames.ControlFileParams.spawningBlockname],\
        parsedJSON[blockNames.ControlFileParams.simulationBlockname], parsedJSON[blockNames.ControlFileParams.clusteringBlockname]


def saveInitialControlFile(jsonParams, originalControlFile):
    """
        Save the adaptive control file jsonParams in originalControlFile

        :param jsonParams: Input control file in JSON format
        :type jsonParams: str
        :param originalControlFile: Path where to save the control file
        :type originalControlFile: str
    """
    with open(originalControlFile, 'w') as f:
        with open(jsonParams, 'r') as fr:
            jsonFile = fr.read()
        f.write(jsonFile)


def needToRecluster(oldClusteringMethod, newClusteringMethod):
    """
        Check if the parameters have changed in a restart and we need to redo
        the clustering. In particular: type of clustering, theshold calculator
        or distance

        :param oldClusteringMethod: Clustering in a previous simulation before the restart
        :type oldClusteringMethod: :py:class:`.Clustering`
        :param newClusteringMethod: Clustering in the restarted simulation
        :type newClusteringMethod: :py:class:`.Clustering`

        :returns: bool -- If clustering needs to be redone
    """

    # Check 1: change of type
    if oldClusteringMethod.type != newClusteringMethod.type:
        return True

    # Check 2: Change of thresholdCalculator and thresholdDistance
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.rmsd or\
       oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMap:
        return oldClusteringMethod.thresholdCalculator != newClusteringMethod.thresholdCalculator or\
                abs(oldClusteringMethod.contactThresholdDistance - newClusteringMethod.contactThresholdDistance) > 1e-7

    # Check 3: Change of similarity Evaluator in contactMap clustering
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMap:
        return oldClusteringMethod.similarityEvaluator.typeEvaluator != newClusteringMethod.similarityEvaluator.typeEvaluator


def clusterEpochTrajs(clusteringMethod, epoch, epochOutputPathTempletized, topologies):
    """
        Cluster the trajecotories of a given epoch

        :param clusteringMethod: Clustering object
        :type clusteringMethod: :py:class:`.Clustering`
        :param epoch: Number of the epoch to cluster
        :type epoch: int
        :param epochOutputPathTempletized: Path where to find the trajectories
        :type epochOutputPathTempletized: str
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`
"""

    snapshotsJSONSelectionString = generateTrajectorySelectionString(epoch, epochOutputPathTempletized)
    paths = ast.literal_eval(snapshotsJSONSelectionString)
    if len(glob.glob(paths[-1])) == 0:
        sys.exit("No trajectories to cluster! Matching path:%s" % paths[-1])
    clusteringMethod.cluster(paths, topology=topologies, epoch=epoch)


def clusterPreviousEpochs(clusteringMethod, finalEpoch, epochOutputPathTempletized, simulationRunner, topologies):
    """
        Cluster all previous epochs using the clusteringMethod object

        :param clusteringMethod: Clustering object
        :type clusteringMethod: :py:class:`.Clustering`
        :param finalEpoch: Last epoch to cluster (not included)
        :type finalEpoch: int
        :param epochOutputPathTempletized: Path where to find the trajectories
        :type epochOutputPathTempletized: str
        :param simulationRunner: Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
        :param topologies: Topology object containing the set of topologies needed for the simulatioies        :type topologies: :py:class:`.Topology`
"""
    for i in range(finalEpoch):
        simulationRunner.readMappingFromDisk(epochOutputPathTempletized % i)
        topologies.readMappingFromDisk(epochOutputPathTempletized % i, i)
        clusterEpochTrajs(clusteringMethod, i, epochOutputPathTempletized, topologies)


def getWorkingClusteringObjectAndReclusterIfNecessary(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner, topologies):
    """
        It reads the previous clustering method, and, if there are changes,
        it reclusters the previous trajectories. Returns the clustering object to use

        :param firstRun: New epoch to run
        :type firstRun: int
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: :py:class:`.OutputPathConstants`
        :param clusteringBlock: Contains the new clustering block
        :type clusteringBlock: json
        :param spawningParams: Spawning params, to know what reportFile and column to read
        :type spawningParams: :py:class:`.SpawningParams`
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`

        :returns: :py:class:`.Clustering` -- The clustering method to use in the
            adaptive sampling simulation
    """

    lastClusteringEpoch = firstRun - 1
    clusteringObjectPath = outputPathConstants.clusteringOutputObject % (lastClusteringEpoch)
    oldClusteringMethod = utilities.readClusteringObject(clusteringObjectPath)

    clusteringBuilder = clustering.ClusteringBuilder()
    clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                         spawningParams.reportFilename,
                                                         spawningParams.reportCol)

    if needToRecluster(oldClusteringMethod, clusteringMethod):
        print("Reclustering!")
        startTime = time.time()
        clusterPreviousEpochs(clusteringMethod, firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner, topologies)
        endTime = time.time()
        print("Reclustering took %s sec" % (endTime - startTime))
    else:
        clusteringMethod = oldClusteringMethod
        clusteringMethod.setCol(spawningParams.reportCol)
        for ij in range(firstRun):
            topologies.readMappingFromDisk(outputPathConstants.epochOutputPathTempletized % ij, ij)

    return clusteringMethod


def buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock,
                                                         spawningParams, spawningCalculator, simulationRunner, topologies, processManager):
    """
        It reads the previous clustering method, and if there are changes (clustering method or related to thresholds),
        reclusters the previous trajectories. Returns the clustering object to use,
        and the initial structure filenames as strings

        :param firstRun: New epoch to run
        :type firstRun: int
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: str
        :param clusteringBlock: Contains the new clustering block
        :type clusteringBlock: json
        :param spawningParams: Spawning params
        :type spawningParams: :py:class:`.SpawningParams`
        :param spawningCalculator: Spawning calculator object
        :type spawningCalculator: :py:class:`.SpawningCalculator`
        :param simulationRunner: :py:class:`.SimulationRunner` Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`
        :param processManager: Object to synchronize the possibly multiple processes
        :type processManager: :py:class:`.ProcessesManager`

        :returns: :py:class:`.Clustering`, str -- The clustering method to use in the adaptive sampling simulation and the initial structures filenames
    """
    processManager.setStatus(processManager.RUNNING)
    if processManager.isMaster():
        clusteringMethod = getWorkingClusteringObjectAndReclusterIfNecessary(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner, topologies)

        degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.getWorkingProcessors(), firstRun)
        spawningCalculator.log()
        _, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, firstRun, topologies=topologies)
        # for compatibility with old data
        procMapping = [element if element is not None else (0, 0, 0) for element in procMapping]
        topologies.mapEpochTopologies(firstRun, procMapping)
        simulationRunner.updateMappingProcessors(procMapping)
    else:
        clusteringMethod = None
    processManager.setStatus(processManager.WAITING)
    processManager.synchronize()
    initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(simulationRunner.getWorkingProcessors(), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

    return clusteringMethod, initialStructuresAsString


def buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, controlFile, outputPathConstants, clusteringBlock, spawningParams, initialStructures, simulationRunner, processManager):
    """
        Build the clustering object and copies initial structures from control file.
        Returns the clustering object to use and the initial structures filenames as string

        :param debug: In debug, it will not remove the simulations
        :type debug: bool
        :param controlFile: Adaptive sampling control file
        :type controlFile: str
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: :py:class:`.OutputPathConstants`
        :param clusteringBlock: Contains the new clustering block
        :type clusteringBlock: json
        :param spawningParams: Spawning params
        :type spawningParams: :py:class:`.SpawningParams`
        :param initialStructures: Control file initial structures
        :type initialStructures: list
        :param simulationRunner: :py:class:`.SimulationRunner` Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
        :param processManager: Object to synchronize the possibly multiple processes
        :type processManager: :py:class:`.ProcessesManager`

        :returns: :py:class:`.Clustering`, str -- The clustering method to use in the adaptive sampling simulation and the initial structures filenames
    """
    firstRun = 0
    processManager.setStatus(processManager.RUNNING)
    if processManager.isMaster():
        saveInitialControlFile(controlFile, outputPathConstants.originalControlFile)

        copyInitialStructures(initialStructures, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    processManager.setStatus(processManager.WAITING)
    processManager.synchronize()
    initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(len(initialStructures), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

    if processManager.isMaster():
        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                             spawningParams.reportFilename,
                                                             spawningParams.reportCol)
    else:
        clusteringMethod = None
    return clusteringMethod, initialStructuresAsString


def main(jsonParams, clusteringHook=None):
    """
        Main body of the adaptive sampling program.

        :param jsonParams: A string with the name of the control file to use
        :type jsonParams: str
    """

    controlFileValidator.validate(jsonParams)
    generalParams, spawningBlock, simulationrunnerBlock, clusteringBlock = loadParams(jsonParams)

    spawningAlgorithmBuilder = spawning.SpawningAlgorithmBuilder()
    spawningCalculator = spawningAlgorithmBuilder.build(spawningBlock)

    runnerbuilder = simulationrunner.RunnerBuilder()
    simulationRunner = runnerbuilder.build(simulationrunnerBlock)

    restart = generalParams.get(blockNames.GeneralParams.restart, True)
    debug = generalParams.get(blockNames.GeneralParams.debug, False)
    outputPath = generalParams[blockNames.GeneralParams.outputPath]
    initialStructuresWildcard = generalParams[blockNames.GeneralParams.initialStructures]
    writeAll = generalParams.get(blockNames.GeneralParams.writeAllClustering, False)
    nativeStructure = generalParams.get(blockNames.GeneralParams.nativeStructure, '')
    resname = clusteringBlock[blockNames.ClusteringTypes.params].get(blockNames.ClusteringTypes.ligandResname)

    print("================================")
    print("            PARAMS              ")
    print("================================")
    print("Restarting simulations", restart)
    print("Debug:", debug)
    print("Iterations: %d, Mpi processors: %d, Pele steps: %d" % (simulationRunner.parameters.iterations, simulationRunner.parameters.processors, simulationRunner.parameters.peleSteps))
    print("SpawningType:", spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type])
    print("SimulationType:", simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY[simulationRunner.type])
    if simulationRunner.hasExitCondition():
        print("Exit condition:", simulationTypes.EXITCONDITION_TYPE_TO_STRING_DICTIONARY[simulationRunner.parameters.exitCondition.type])
    print("Clustering method:", clusteringBlock[blockNames.ClusteringTypes.type])
    print("Output path: ", outputPath)
    print("Initial Structures: ", initialStructuresWildcard)
    print("================================\n\n")

    initialStructures = expandInitialStructuresWildcard(initialStructuresWildcard)
    if not initialStructures:
        raise InitialStructuresError("No initial structures found!!!")

    if len(initialStructures) > simulationRunner.getWorkingProcessors():
        raise InitialStructuresError("Error: More initial structures than Working Processors found!!!")

    if resname is not None:
        checkSymmetryDict(clusteringBlock, initialStructures, resname)

    outputPathConstants = constants.OutputPathConstants(outputPath)

    if not debug:
        atexit.register(utilities.cleanup, outputPathConstants.tmpFolder)

    simulationRunner.unifyReportNames(spawningCalculator.parameters.reportFilename)
    utilities.makeFolder(outputPath)
    processManager = ProcessesManager(outputPath)
    firstRun = findFirstRun(outputPath, outputPathConstants.clusteringOutputObject, simulationRunner)
    processManager.setStatus(processManager.RUNNING)
    if processManager.isMaster():
        utilities.makeFolder(outputPathConstants.tmpFolder)
        utilities.makeFolder(outputPathConstants.topologies)
        saveInitialControlFile(jsonParams, outputPathConstants.originalControlFile)
    processManager.setStatus(processManager.WAITING)
    processManager.synchronize()

    topologies = utilities.Topology(outputPathConstants.topologies)
    if restart and firstRun is not None:
        topology_files = glob.glob(os.path.join(outputPathConstants.topologies, "topology*.pdb"))
        topologies.setTopologies(topology_files)
        if firstRun == 0:
            createMappingForFirstEpoch(initialStructures, topologies, simulationRunner.getWorkingProcessors())
            clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, jsonParams, outputPathConstants, clusteringBlock, spawningCalculator.parameters, initialStructures, simulationRunner, processManager)
        else:
            clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock, spawningCalculator.parameters, spawningCalculator, simulationRunner, topologies, processManager)
        processManager.setStatus(processManager.RUNNING)
        if processManager.isMaster():
            checkMetricExitConditionMultipleTrajsinRestart(firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner)
        processManager.setStatus(processManager.WAITING)
        processManager.synchronize()

    if firstRun is None or not restart:
        processManager.setStatus(processManager.RUNNING)
        topologies.setTopologies(initialStructures)
        if processManager.isMaster():
            if not debug:
                cleanPreviousSimulation(outputPath)
            writeTopologyFiles(initialStructures, outputPathConstants.topologies)
        processManager.setStatus(processManager.WAITING)
        processManager.synchronize()

        if simulationRunner.parameters.runEquilibration:
            processManager.setStatus(processManager.RUNNING)
            if resname is None:
                raise utilities.RequiredParameterMissingException("Resname not specified in clustering block!!!")
            initialStructures = simulationRunner.equilibrate(initialStructures, outputPathConstants, spawningCalculator.parameters.reportFilename, outputPath, resname, topologies)
            processManager.setStatus(processManager.WAITING)
            processManager.synchronize()
            topologies.setTopologies(initialStructures)
            writeTopologyFiles(initialStructures, outputPathConstants.topologies)
        createMappingForFirstEpoch(initialStructures, topologies, simulationRunner.getWorkingProcessors())

        clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, jsonParams, outputPathConstants, clusteringBlock, spawningCalculator.parameters, initialStructures, simulationRunner, processManager)

    processManager.setStatus(processManager.RUNNING)
    if simulationRunner.parameters.modeMovingBox is not None and simulationRunner.parameters.boxCenter is None:
        simulationRunner.parameters.boxCenter = simulationRunner.selectInitialBoxCenter(initialStructuresAsString, resname)
    if firstRun is None:
        firstRun = 0  # if restart false, but there were previous simulations
    for i in range(firstRun, simulationRunner.parameters.iterations):
        if not processManager.allRunning():
            print("Some process died, killing all replicas!")
            sys.exit(1)
        processManager.setStatus(processManager.RUNNING)
        print("Iteration", i)
        if processManager.isMaster():
            outputDir = outputPathConstants.epochOutputPathTempletized % i
            utilities.makeFolder(outputDir)

            simulationRunner.writeMappingToDisk(outputPathConstants.epochOutputPathTempletized % i)
            topologies.writeMappingToDisk(outputPathConstants.epochOutputPathTempletized % i, i)

        processManager.setStatus(processManager.WAITING)
        processManager.synchronize()
        print("Production run...")
        processManager.setStatus(processManager.RUNNING)
        if not debug:
            simulationRunner.runSimulation(i, outputPathConstants, initialStructuresAsString, topologies, spawningCalculator.parameters.reportFilename)
        processManager.setStatus(processManager.WAITING)
        processManager.synchronize()

        print("Clustering...")
        if processManager.isMaster():
            startTime = time.time()
            clusterEpochTrajs(clusteringMethod, i, outputPathConstants.epochOutputPathTempletized, topologies)
            endTime = time.time()
            print("Clustering ligand: %s sec" % (endTime - startTime))

            if clusteringHook is not None:
                clusteringHook(clusteringMethod, outputPathConstants, simulationRunner, i+1)

        if simulationRunner.parameters.modeMovingBox is not None:
            simulationRunner.getNextIterationBox(outputPathConstants.epochOutputPathTempletized % i, resname, topologies, i)
            if processManager.isMaster():
                clustersList, clustersFiltered = filterClustersAccordingToBox(simulationRunner.parameters, clusteringMethod)
        else:
            if processManager.isMaster():
                clustersList = clusteringMethod.clusters.clusters

        if processManager.isMaster():
            degeneracyOfRepresentatives = spawningCalculator.calculate(clustersList, simulationRunner.getWorkingProcessors(), i)
            spawningCalculator.log()

            if degeneracyOfRepresentatives is not None:
                if simulationRunner.parameters.modeMovingBox is not None:
                    degeneracyOfRepresentatives = mergeFilteredClustersAccordingToBox(degeneracyOfRepresentatives, clustersFiltered)
                print("Degeneracy", degeneracyOfRepresentatives)
                assert len(degeneracyOfRepresentatives) == len(clusteringMethod.clusters.clusters)
            else:
                # When using null or independent spawning the calculate method returns None
                assert spawningCalculator.type in spawningTypes.SPAWNING_NO_DEGENERACY_TYPES, "calculate returned None with spawning type %s" % spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type]

            clusteringMethod.writeOutput(outputPathConstants.clusteringOutputDir % i,
                                         degeneracyOfRepresentatives,
                                         outputPathConstants.clusteringOutputObject % i, writeAll)
            simulationRunner.cleanCheckpointFiles(outputPathConstants.epochOutputPathTempletized % i)

            if i > 0:
                # Remove old clustering object, since we already have a newer one
                try:
                    os.remove(outputPathConstants.clusteringOutputObject % (i-1))
                except OSError:
                    # In case of restart
                    pass

        # Prepare for next pele iteration
        if i != simulationRunner.parameters.iterations-1:
            if spawningCalculator.shouldWriteStructures():
                # Differentiate between null spawning and the rest of spawning
                # methods
                if processManager.isMaster():
                    _, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants,
                                                                                       degeneracyOfRepresentatives,
                                                                                       clusteringMethod, i+1,
                                                                                       topologies=topologies)
                    simulationRunner.updateMappingProcessors(procMapping)
                    topologies.mapEpochTopologies(i+1, procMapping)
                initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(simulationRunner.getWorkingProcessors(),
                                                                                              outputPathConstants.tmpInitialStructuresTemplate,
                                                                                              i+1)

        if processManager.isMaster():
            if clusteringMethod.symmetries and nativeStructure:
                fixReportsSymmetry(outputPathConstants.epochOutputPathTempletized % i, resname,
                                   nativeStructure, clusteringMethod.symmetries, topologies)

            # check exit condition, if defined
            if simulationRunner.hasExitCondition():
                if simulationRunner.checkExitCondition(clusteringMethod, outputPathConstants.epochOutputPathTempletized % i):
                    print("Simulation exit condition met at iteration %d, stopping" % i)
                    processManager.setStatus(processManager.WAITING)
                    processManager.synchronize()
                    # send a signal to all possible adaptivePELE copies to stop
                    for pid in processManager.lockInfo:
                        if pid != processManager.pid:
                            os.kill(pid, signal.SIGTERM)
                    break
                else:
                    print("Simulation exit condition not met at iteration %d, continuing..." % i)
        processManager.setStatus(processManager.WAITING)
        processManager.synchronize()

if __name__ == '__main__':
    args = parseArgs()
    main(args.controlFile)
