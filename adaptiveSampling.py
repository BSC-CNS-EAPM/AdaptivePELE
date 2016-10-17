import numpy as np
import sys
import shutil
import os
import json
import time
import glob
import pickle
from constants import blockNames, constants
from atomset import atomset
from utilities import utilities
from validator import controlFileValidator
from spawning import spawning, spawningTypes
from simulation import simulationrunner, simulationTypes
from clustering import clustering, clusteringTypes


def checkSymmetryDict(clusteringBlock, initialStructures, resname):
    utilities.generateReciprocalAtoms(clusteringBlock[blockNames.ClusteringTypes.params].get(blockNames.ClusteringTypes.symmetries, {}))
    symmetries = clusteringBlock[blockNames.ClusteringTypes.params].get(blockNames.ClusteringTypes.symmetries, {})
    for structure in initialStructures:
        PDB = atomset.PDB()
        PDB.initialise(structure, resname=resname)
        utilities.assertSymmetriesDict(symmetries, PDB)


def fixReportsSymmetry(outputPath, resname, nativeStructure, symmetries):
    outputFilename = "fixedReport_%d"
    trajName = "*traj*.pdb"
    reportName = "*report_%d"
    trajs = glob.glob(os.path.join(outputPath, trajName))
    nativePDB = atomset.PDB()
    nativePDB.initialise(nativeStructure, resname=resname)
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
    for i, name in enumerate(initialStructures):
        shutil.copyfile(name, tmpInitialStructuresTemplate % (iteration, i))


def createMultipleComplexesFilenames(numberOfSnapshots, tmpInitialStructuresTemplate, iteration):
    jsonString = "\n"
    for i in range(numberOfSnapshots-1):
        jsonString += constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, i)) + ",\n"
    jsonString += constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, numberOfSnapshots-1))
    return jsonString


def generateSnapshotSelectionStringLastRound(currentEpoch, outputPathTempletized):
    return " \"" + os.path.join(outputPathTempletized % currentEpoch, constants.trajectoryBasename) + "\""


def writeInitialStructuresAccordingToDegeneracy(tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clustering, iteration):

    counts = 0
    for i, cluster in enumerate(clustering.clusters.clusters):
        for j in range(int(degeneracyOfRepresentatives[i])):
            outputFilename = tmpInitialStructuresTemplate % (iteration, counts)
            print 'Writing to ', outputFilename, 'cluster', i
            cluster.writePDB(outputFilename)

            counts += 1

    print "counts & cluster centers", counts, np.where(np.array(degeneracyOfRepresentatives) > 0)[0].size
    return counts


def findFirstRun(outputPath, clusteringOutputObject):
    """ Assumes that the outputPath is XXX/%d """

    folderWithSimulationData = outputPath
    allFolders = os.listdir(folderWithSimulationData)
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort(reverse=True)

    for epoch in epochFolders:
        if os.path.exists(clusteringOutputObject % epoch):
            return epoch + 1
    return 0


def loadParams(jsonParams):
    """
        TODO: change for variables in a block names file,
        and work it out a bit more
    """
    jsonFile = open(jsonParams, 'r').read()
    parsedJSON = json.loads(jsonFile)

    return parsedJSON[blockNames.ControlFileParams.generalParams], parsedJSON[blockNames.ControlFileParams.spawningBlockname],\
        parsedJSON[blockNames.ControlFileParams.simulationBlockname], parsedJSON[blockNames.ControlFileParams.clusteringBlockname]


def saveInitialControlFile(jsonParams, originalControlFile):
    file = open(originalControlFile, 'w')
    jsonFile = open(jsonParams, 'r').read()
    file.write(jsonFile)

def needToRecluster(oldClusteringMethod, newClusteringMethod):
    #Check 1: change of type
    if oldClusteringMethod.type != newClusteringMethod.type:
        return True

    #Check 2: Change of thresholdCalculator and thresholdDistance
    if  oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contacts or\
        oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMapAccumulative:
            return oldClusteringMethod.thresholdCalculator != newClusteringMethod.thresholdCalculator or\
                abs(oldClusteringMethod.contactThresholdDistance - newClusteringMethod.contactThresholdDistance) > 1e-7

    #Check 3: Change of nClusters in contactMapAgglomerative
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMapAgglomerative:
        return oldClusteringMethod.nclusters != newClusteringMethod.nclusters

def clusterEpochTrajs(clusteringMethod, epoch, outputPathTempletized):
    snapshotsJSONSelectionString = generateSnapshotSelectionStringLastRound(epoch, outputPathTempletized)
    snapshotsJSONSelectionString = "[" + snapshotsJSONSelectionString + "]"
    paths = eval(snapshotsJSONSelectionString)
    if len(glob.glob(paths[-1])) == 0:
        sys.exit("No more trajectories to cluster")
    clusteringMethod.cluster(paths)

def clusterPreviousEpochs(clusteringMethod, finalEpoch, outputPathTempletized):
    for i in range(finalEpoch):
        clusterEpochTrajs(clusteringMethod, i, outputPathTempletized)


def main(jsonParams=None):
    if jsonParams is None:
        jsonParams = sys.argv[1]

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
    initialStructures = generalParams[blockNames.GeneralParams.initialStructures]
    writeAll = generalParams[blockNames.GeneralParams.writeAllClustering]
    nativeStructure = generalParams.get(blockNames.GeneralParams.nativeStructure, '')
    resname = clusteringBlock[blockNames.ClusteringTypes.params][blockNames.ClusteringTypes.ligandResname]


    print "================================"
    print "            PARAMS              "
    print "================================"
    print "Restarting simulations", generalParams[blockNames.GeneralParams.restart]
    print "Debug:", generalParams[blockNames.GeneralParams.debug]

    print "Iterations: %d, Mpi processors: %d, Pele steps: %d" % (simulationRunner.parameters.iterations, simulationRunner.parameters.processors, simulationRunner.parameters.peleSteps)

    print "SpawningType:", spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type]

    print "SimulationType:", simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY[simulationRunner.type]
    print "Clustering method:", clusteringType

    print "Output path: ", outputPath
    print "Initial Structures: ", initialStructures
    print "================================\n\n"

    checkSymmetryDict(clusteringBlock, initialStructures, resname)

    outputPathConstants = constants.OutputPathConstants(outputPath)

    # if not debug: atexit.register(utilities.cleanup, outputPathConstants.tmpFolder)

    utilities.makeFolder(outputPath)
    utilities.makeFolder(outputPathConstants.tmpFolder)
    saveInitialControlFile(jsonParams, outputPathConstants.originalControlFile)

    # print variable epsilon information
    if spawningParams.epsilon is not None:
        epsilon_file = open("epsilon_values.txt", "w")
        epsilon_file.write("Iteration\tEpsilon\n")
        epsilon_file.close()

    if restart:
        firstRun = findFirstRun(outputPath, outputPathConstants.clusteringOutputObject)

        if firstRun != 0:
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
                clusterPreviousEpochs(clusteringMethod, firstRun, outputPathConstants.outputPathTempletized)
                endTime = time.time()
                print "Reclustering took %s sec" % (endTime - startTime)
            else:
                clusteringMethod = oldClusteringMethod
                clusteringMethod.setCol(spawningParams.reportCol)

            degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, firstRun)
            spawningCalculator.log()
            print "Degeneracy", degeneracyOfRepresentatives

            seedingPoints = writeInitialStructuresAccordingToDegeneracy(outputPathConstants.tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clusteringMethod, firstRun)

            initialStructuresAsString = createMultipleComplexesFilenames(seedingPoints, outputPathConstants.tmpInitialStructuresTemplate, firstRun)

    if not restart or firstRun == 0:
        # if restart and firstRun = 0, it must check into the initial structures
        # Choose initial structures
        if not debug:
            shutil.rmtree(outputPath)
        utilities.makeFolder(outputPath)
        firstRun = 0
        saveInitialControlFile(jsonParams, outputPathConstants.originalControlFile)
        copyInitialStructures(initialStructures, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
        initialStructuresAsString = createMultipleComplexesFilenames(len(initialStructures), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                             spawningParams.reportFilename,
                                                             spawningParams.reportCol)

    peleControlFileDictionary = {"COMPLEXES": initialStructuresAsString, "PELE_STEPS": simulationRunner.parameters.peleSteps}

    # Make control file
    outputDir = outputPathConstants.outputPathTempletized % firstRun
    utilities.makeFolder(outputDir)
    peleControlFileDictionary["OUTPUT_PATH"] = outputDir
    peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + firstRun*simulationRunner.parameters.processors
    # simulationRunner.parameters needn't be passed to makeWorkingControlFile
    simulationRunner.makeWorkingControlFile(simulationRunner.parameters.templetizedControlFile, outputPathConstants.tmpControlFilename % firstRun, peleControlFileDictionary)

    for i in range(firstRun, simulationRunner.parameters.iterations):
        print "Iteration", i
        if spawningParams.epsilon is not None:
            epsilon_file = open("epsilon_values.txt", "a")
            epsilon_file.write("%d\t%f\n" % (i, spawningParams.epsilon))
            epsilon_file.close()

        print "Production run..."
        if not debug:
            startTime = time.time()
            simulationRunner.runSimulation(outputPathConstants.tmpControlFilename % i)
            endTime = time.time()
            print "PELE %s sec" % (endTime - startTime)

        print "Clustering..."
        startTime = time.time()
        clusterEpochTrajs(clusteringMethod, i, outputPathConstants.outputPathTempletized)
        endTime = time.time()
        print "Clustering ligand: %s sec" % (endTime - startTime)

        degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, i)
        spawningCalculator.log()
        print "Degeneracy", degeneracyOfRepresentatives

        numberOfSeedingPoints = writeInitialStructuresAccordingToDegeneracy(outputPathConstants.tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clusteringMethod, i+1)

        clusteringMethod.writeOutput(outputPathConstants.clusteringOutputDir % i,
                                     degeneracyOfRepresentatives,
                                     outputPathConstants.clusteringOutputObject % i, writeAll)
        if i > 0:
            # Remove old clustering object, since we already have a newer one
            try:
                os.remove(outputPathConstants.clusteringOutputObject%(i-1))
            except OSError:
                # In case of restart
                pass

        # Prepare for next pele iteration
        if i != simulationRunner.parameters.iterations-1:
            initialStructuresAsString = createMultipleComplexesFilenames(numberOfSeedingPoints, outputPathConstants.tmpInitialStructuresTemplate, i+1)
            peleControlFileDictionary["COMPLEXES"] = initialStructuresAsString

            outputDir = outputPathConstants.outputPathTempletized % (i+1)
            utilities.makeFolder(outputDir)  # PELE does not do it automatically
            peleControlFileDictionary["OUTPUT_PATH"] = outputDir
            peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + (i+1)*simulationRunner.parameters.processors
            # simulationRunner.parameters needn't be passed to makeWorkingControlFile
            simulationRunner.makeWorkingControlFile(simulationRunner.parameters.templetizedControlFile, outputPathConstants.tmpControlFilename % (i+1), peleControlFileDictionary)

        if clusteringMethod.symmetries and nativeStructure:
            fixReportsSymmetry(outputPathConstants.outputPathTempletized % i, resname,
                               nativeStructure, clusteringMethod.symmetries)
    # utilities.cleanup
    # utilities.cleanup(outputPathConstants.tmpFolder)

if __name__ == '__main__':
    main()
