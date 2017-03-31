import numpy as np
import sys
import shutil
import os
import json
import time
import glob
from AdaptivePELE.constants import blockNames, constants
from AdaptivePELE.atomset import atomset
from AdaptivePELE.utilities import utilities
from AdaptivePELE.validator import controlFileValidator
from AdaptivePELE.spawning import spawning, spawningTypes
from AdaptivePELE.simulation import simulationrunner, simulationTypes
from AdaptivePELE.clustering import clustering, clusteringTypes
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description="Perform several iterations of"
                                     " simulations using adaptive sampling to "
                                     "distribute the processors in order "
                                     "to optimize sampling")
    parser.add_argument('controlFile', type=str)
    args = parser.parse_args()
    return args

def expandInitialStructuresWildcard(initialStructuresWildcard):
    """
        Returns the initial structures after expanding the initial structures wildcard
    """
    totalInitialStructures = []
    for initialStructureWildcard in initialStructuresWildcard:
        expandedStructures = glob.glob(initialStructureWildcard)
        totalInitialStructures.extend(expandedStructures)
    return totalInitialStructures

def checkSymmetryDict(clusteringBlock, initialStructures, resname):
    """ Check if the symmetries dictionary is valid for the system

        clusteringBlock [In] Block with the clustering-related parameters
        initialStructures [In] List of the initial structures
        resname [In] Residue name of the ligand in the system pdb
     """
    symmetries = clusteringBlock[blockNames.ClusteringTypes.params].get(blockNames.ClusteringTypes.symmetries, {})
    for structure in initialStructures:
        PDB = atomset.PDB()
        PDB.initialise(str(structure), resname=resname)
        utilities.assertSymmetriesDict(symmetries, PDB)


def fixReportsSymmetry(outputPath, resname, nativeStructure, symmetries):
    """
        Correct the RMSD for symmetries in the reports. New reports are stored
        in the fixedReport_i where i is the number of the report

        outputPath [In] Path where to find the trajectories
        resname [In] Residue name of the ligand in the pdb
        nativeStructure [In] Path to the native structure pdb
        symmetries [In] Dictionary containg the symmetries of the ligand
    """
    outputFilename = "fixedReport_%d"
    trajName = "*traj*.pdb"
    reportName = "*report_%d"
    trajs = glob.glob(os.path.join(outputPath, trajName))
    nativePDB = atomset.PDB()
    nativePDB.initialise(str(nativeStructure), resname=resname)
    for traj in trajs:
        trajNum = utilities.getTrajNum(traj)
        rmsd = list(utilities.getRMSD(traj, nativePDB, resname, symmetries))
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
       Copy the initial structures from a certain iteration

       initialStructures [In] Name of the initial structures to copy
       tmpInitialStructuresTemplate [In] Template with the name of the initial
       strutctures
       iteration [In] Number of epoch
    """

    for i, name in enumerate(initialStructures):
        shutil.copyfile(name, tmpInitialStructuresTemplate % (iteration, i))


def createMultipleComplexesFilenames(numberOfSnapshots, tmpInitialStructuresTemplate, iteration):
    """
       Create the string to substitute the complexes in the pele control file

       numberOfSnapshots [In] Number of complexes to write
       tmpInitialStructuresTemplate [In] Template with the name of the initial
       strutctures
       iteration [In] Number of epoch
    """
    jsonString = "\n"
    for i in range(numberOfSnapshots-1):
        jsonString += constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, i)) + ",\n"
    jsonString += constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, numberOfSnapshots-1))
    return jsonString


def generateSnapshotSelectionStringLastRound(currentEpoch, epochOutputPathTempletized):
    """ Generate the template for the name of the trajectories in the current
        epoch
        currentEpoch [In] Epoch number
        epochOutputPathTempletized [In] Template for the path where the
        trajectories of any epoch are stored
    """
    return " \"" + os.path.join(epochOutputPathTempletized % currentEpoch, constants.trajectoryBasename) + "\""


def findFirstRun(outputPath, clusteringOutputObject):
    """
        Find the last epoch that was properly simulated and clusterized and
        and return the first epoch to run in case of restart

        outputPath [In] Simulation output path
        clusteringOutputObject [In] Templetized name of the clustering object
    """

    folderWithSimulationData = outputPath
    allFolders = os.listdir(folderWithSimulationData)
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort(reverse=True)

    objectsFound = []
    for epoch in epochFolders:
        if os.path.exists(clusteringOutputObject % epoch):
            objectsFound.append(epoch)
        if objectsFound and epoch < (objectsFound[0]-5):
            break
    while objectsFound:
        epoch = objectsFound.pop(0)
        if checkIntegrityClusteringObject(clusteringOutputObject % epoch):
            return epoch + 1
    return 0


def checkIntegrityClusteringObject(objectPath):
    """
        Test wheter the found clustering object to reload is a valid object

        :param objectPath: Path to the clustering object
        :param objectPath: str
        :returns: bool -- True if the clustering object found is valid
    """
    try:
        utilities.readClusteringObject(objectPath)
        return True
    except EOFError:
        return False


def loadParams(jsonParams):
    """
        Read the control file in JSON format and extract the blocks of simulation,
        general parameters, clustering and spawning

        jsonParams [In] Control file in JSON format from where the parameters
        will be read
    """
    jsonFile = open(jsonParams, 'r').read()
    parsedJSON = json.loads(jsonFile)

    return parsedJSON[blockNames.ControlFileParams.generalParams], parsedJSON[blockNames.ControlFileParams.spawningBlockname],\
        parsedJSON[blockNames.ControlFileParams.simulationBlockname], parsedJSON[blockNames.ControlFileParams.clusteringBlockname]


def saveInitialControlFile(jsonParams, originalControlFile):
    """ Save the control file jsonParams in originalControlFile

        jsonParams [In] Input control file in JSON format
        originalControlFile [In] Path where to save the control file"""
    file = open(originalControlFile, 'w')
    jsonFile = open(jsonParams, 'r').read()
    file.write(jsonFile)


def needToRecluster(oldClusteringMethod, newClusteringMethod):
    """ Check if the parameters have changed in a restart and we need to redo
        the clustering

        oldClusteringMethod [In] clusteringMethod in the previous simulation
        newClusteringMethod [In] clusteringMethod specified in the restart simulation
    """
    # Check 1: change of type
    if oldClusteringMethod.type != newClusteringMethod.type:
        return True

    # Check 2: Change of thresholdCalculator and thresholdDistance
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contacts or\
       oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMapAccumulative:
        return oldClusteringMethod.thresholdCalculator != newClusteringMethod.thresholdCalculator or\
                abs(oldClusteringMethod.contactThresholdDistance - newClusteringMethod.contactThresholdDistance) > 1e-7

    # Check 3: Change of nClusters in contactMapAgglomerative
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMapAgglomerative:
        return oldClusteringMethod.nclusters != newClusteringMethod.nclusters


def clusterEpochTrajs(clusteringMethod, epoch, epochOutputPathTempletized, processorsToClusterMapping):
    """ Cluster the trajecotories of a given epoch into the clusteringMethod
        object
        clusteringMethod [In] Object that holds all the clustering-related
        information
        epoch [In] Number of the epoch to cluster
        epochOutputPathTempletized [In] Path where to find the trajectories"""

    snapshotsJSONSelectionString = generateSnapshotSelectionStringLastRound(epoch, epochOutputPathTempletized)
    snapshotsJSONSelectionString = "[" + snapshotsJSONSelectionString + "]"
    paths = eval(snapshotsJSONSelectionString)
    if len(glob.glob(paths[-1])) == 0:
        sys.exit("No more trajectories to cluster")
    clusteringMethod.cluster(paths, processorsToClusterMapping)


def clusterPreviousEpochs(clusteringMethod, finalEpoch, epochOutputPathTempletized, simulationRunner):
    """ Cluster the epoch until finalEpoch into the clusteringMethod

        clusteringMethod [In] Object that holds all the clustering-related
        information
        finalEpoch [In] Number of the final epoch to cluster (Not included)
        epochOutputPathTempletized [In] Path where to find the trajectories"""
    for i in range(finalEpoch):
        simulationRunner.readMappingFromDisk(epochOutputPathTempletized % i)
        clusterEpochTrajs(clusteringMethod, i, epochOutputPathTempletized, simulationRunner.processorsToClusterMapping)


def getWorkingClusteringObject(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner):
    """
        It reads the previous clustering method, and, if there are changes,
        it reclusters the previous trajectories. Returns the clustering object to use

        firstRun [In] New epoch to run
        outputPathConstants [In] Contains outputPath-related constants
        clusteringBlock [In] Contains the new clustering block
        spawningParams [In] Contains spawning params, to know what reportFile and column to read

        returns clusteringMethod, the clustering method to use in the adaptive sampling simulation
    """

    lastClusteringEpoch = firstRun - 1
    clusteringObjectPath = outputPathConstants.clusteringOutputObject % (lastClusteringEpoch)
    oldClusteringMethod = utilities.readClusteringObject(clusteringObjectPath)

    clusteringBuilder = clustering.ClusteringBuilder()
    clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                         spawningParams.reportFilename,
                                                         spawningParams.reportCol)

    if needToRecluster(oldClusteringMethod, clusteringMethod):
        print "Reclustering!"
        startTime = time.time()
        clusterPreviousEpochs(clusteringMethod, firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner)
        endTime = time.time()
        print "Reclustering took %s sec" % (endTime - startTime)
    else:
        clusteringMethod = oldClusteringMethod
        clusteringMethod.setCol(spawningParams.reportCol)

    return clusteringMethod

def buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock,
        spawningParams, spawningCalculator, simulationRunner):
    """
        It reads the previous clustering method, and, if there are changes,
        it reclusters the previous trajectories. Returns the clustering object to use,
        and the initial structure filenames as strings

        firstRun [In] New epoch to run
        outputPathConstants [In] Contains outputPath-related constants
        clusteringBlock [In] Contains the new clustering block
        spawningParams [In] Spawning params
        spawningCalculator [In] Spawning calculator object
        simulationRunner [In] Simulation runner object

        returns clusteringMethod, the clustering method to use in the adaptive sampling simulation
        returns initialStructuresAsString, the initial structures filenames in a string
    """
    clusteringMethod = getWorkingClusteringObject(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner)

    if not hasattr(clusteringMethod, "conformationNetwork"):
        clusteringMethod.epoch = firstRun

    degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, firstRun)
    spawningCalculator.log()
    print "Degeneracy", degeneracyOfRepresentatives
    seedingPoints = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, firstRun)
    initialStructuresAsString = createMultipleComplexesFilenames(seedingPoints, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    simulationRunner.updateMappingProcessors(degeneracyOfRepresentatives, clusteringMethod)

    return clusteringMethod, initialStructuresAsString

def buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, outputPath, controlFile, outputPathConstants, clusteringBlock, spawningParams, initialStructures):
    """
        Build the clustering object and copies initial structures from control file.
        Returns the clustering object to use and the initial structures filenames as string

        debug [In] Whether to delete or not the previous simulation
        outputPath [In] Simulation output path
        controlFile [In] Adaptive sampling control file
        outputPathConstants [In] Contains outputPath-related constants
        clusteringBlock [In] Contains the new clustering block
        spawningParams [In] Spawning params
        initialStructures [In] Control file initial structures

        returns clusteringMethod, the clustering method to use in the adaptive sampling simulation
        returns initialStructuresAsString, the initial structures filenames in a string
    """
    if not debug:
        shutil.rmtree(outputPath)
    utilities.makeFolder(outputPath)
    saveInitialControlFile(controlFile, outputPathConstants.originalControlFile)

    firstRun = 0
    copyInitialStructures(initialStructures, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    initialStructuresAsString = createMultipleComplexesFilenames(len(initialStructures), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

    clusteringBuilder = clustering.ClusteringBuilder()
    clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                         spawningParams.reportFilename,
                                                         spawningParams.reportCol)
    initialClusters = clusteringMethod.clusterInitialStructures(initialStructures)
    return clusteringMethod, initialStructuresAsString, initialClusters


def preparePeleControlFile(i, outputPathConstants, simulationRunner, peleControlFileDictionary):
    """
        Substitute the parameters in the pele control file specified with the
        provided in the control file

        i [In] Epoch number
        outputPathConstants [In] Object that has as attributes constant related
        to the outputPath that will be used to create the working control file
        simulationRunner [In] Simulation runner object
        peleControlFileDictionary [In] Dictonary containing the values of the
        parameters to substitute in the control file
    """
    outputDir = outputPathConstants.epochOutputPathTempletized % i
    utilities.makeFolder(outputDir)
    peleControlFileDictionary["OUTPUT_PATH"] = outputDir
    peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + i*simulationRunner.parameters.processors
    simulationRunner.makeWorkingControlFile(outputPathConstants.tmpControlFilename % i, peleControlFileDictionary)


def main(jsonParams):
    """
        Main body of the adaptive sampling program.

        :params jsonParams: A string with the name of the control file to use
        :type jsonParams: str
    """

    controlFileValidator.validate(jsonParams)
    generalParams, spawningBlock, simulationrunnerBlock, clusteringBlock = loadParams(jsonParams)

    spawningAlgorithmBuilder = spawning.SpawningAlgorithmBuilder()
    spawningCalculator, spawningParams = spawningAlgorithmBuilder.build(spawningBlock)

    runnerbuilder = simulationrunner.RunnerBuilder()
    simulationRunner = runnerbuilder.build(simulationrunnerBlock)

    clusteringType = clusteringBlock[blockNames.ClusteringTypes.type]
    restart = generalParams[blockNames.GeneralParams.restart]
    debug = generalParams[blockNames.GeneralParams.debug]
    outputPath = generalParams[blockNames.GeneralParams.outputPath]
    initialStructuresWildcard = generalParams[blockNames.GeneralParams.initialStructures]
    writeAll = generalParams[blockNames.GeneralParams.writeAllClustering]
    nativeStructure = generalParams.get(blockNames.GeneralParams.nativeStructure, '')
    resname = str(clusteringBlock[blockNames.ClusteringTypes.params][blockNames.ClusteringTypes.ligandResname])


    print "================================"
    print "            PARAMS              "
    print "================================"
    print "Restarting simulations", generalParams[blockNames.GeneralParams.restart]
    print "Debug:", generalParams[blockNames.GeneralParams.debug]

    print "Iterations: %d, Mpi processors: %d, Pele steps: %d" % (simulationRunner.parameters.iterations, simulationRunner.parameters.processors, simulationRunner.parameters.peleSteps)

    print "SpawningType:", spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type]

    print "SimulationType:", simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY[simulationRunner.type]
    if simulationRunner.hasExitCondition():
        print "Exit condition:", simulationTypes.EXITCONDITION_TYPE_TO_STRING_DICTIONARY[simulationRunner.parameters.exitCondition.type]
    print "Clustering method:", clusteringType

    print "Output path: ", outputPath
    print "Initial Structures: ", initialStructuresWildcard
    print "================================\n\n"

    print "wildcard", initialStructuresWildcard
    initialStructures = expandInitialStructuresWildcard(initialStructuresWildcard)
    checkSymmetryDict(clusteringBlock, initialStructures, resname)

    outputPathConstants = constants.OutputPathConstants(outputPath)

    # if not debug: atexit.register(utilities.cleanup, outputPathConstants.tmpFolder)

    utilities.makeFolder(outputPath)
    utilities.makeFolder(outputPathConstants.tmpFolder)
    saveInitialControlFile(jsonParams, outputPathConstants.originalControlFile)

    firstRun = findFirstRun(outputPath, outputPathConstants.clusteringOutputObject)

    if restart and firstRun != 0:
        clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock, spawningParams, spawningCalculator, simulationRunner)
    else:
        firstRun = 0  # if restart false, but there were previous simulations
        clusteringMethod, initialStructuresAsString, initialClusters = buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, outputPath, jsonParams, outputPathConstants, clusteringBlock, spawningParams, initialStructures)
        simulationRunner.makeInitialMapping(initialClusters)


    peleControlFileDictionary = {"COMPLEXES": initialStructuresAsString, "PELE_STEPS": simulationRunner.parameters.peleSteps}

    for i in range(firstRun, simulationRunner.parameters.iterations):
        print "Iteration", i

        print "Preparing control file..."
        preparePeleControlFile(i, outputPathConstants, simulationRunner, peleControlFileDictionary)

        print "Production run..."
        if not debug:
            startTime = time.time()
            simulationRunner.runSimulation(outputPathConstants.tmpControlFilename % i)
            endTime = time.time()
            print "PELE %s sec" % (endTime - startTime)

        simulationRunner.writeMappingToDisk(outputPathConstants.epochOutputPathTempletized % i)

        print "Clustering..."
        startTime = time.time()
        clusterEpochTrajs(clusteringMethod, i, outputPathConstants.epochOutputPathTempletized, simulationRunner.processorsToClusterMapping)
        endTime = time.time()
        print "Clustering ligand: %s sec" % (endTime - startTime)


        degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, i)
        spawningCalculator.log()
        print "Degeneracy", degeneracyOfRepresentatives
        simulationRunner.updateMappingProcessors(degeneracyOfRepresentatives, clusteringMethod)


        clusteringMethod.writeOutput(outputPathConstants.clusteringOutputDir % i,
                                     degeneracyOfRepresentatives,
                                     outputPathConstants.clusteringOutputObject % i, writeAll)

        if i > 0:
            # Remove old clustering object, since we already have a newer one
            try:
                os.remove(outputPathConstants.clusteringOutputObject % (i-1))
            except OSError:
                # In case of restart
                pass

        # Prepare for next pele iteration
        if i != simulationRunner.parameters.iterations-1:
            numberOfSeedingPoints = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, i+1)
            initialStructuresAsString = createMultipleComplexesFilenames(numberOfSeedingPoints, outputPathConstants.tmpInitialStructuresTemplate, i+1)
            peleControlFileDictionary["COMPLEXES"] = initialStructuresAsString

        if clusteringMethod.symmetries and nativeStructure:
            fixReportsSymmetry(outputPathConstants.epochOutputPathTempletized % i, resname,
                               nativeStructure, clusteringMethod.symmetries)

        # check exit condition, if defined
        if simulationRunner.hasExitCondition() and simulationRunner.checkExitCondition(clusteringMethod):
            print "Simulation exit condition met at iteration %d" % i
            break
    # utilities.cleanup
    # utilities.cleanup(outputPathConstants.tmpFolder)

if __name__ == '__main__':
    args = parseArgs()
    main(args.controlFile)
