import numpy as np
import sys
import string
import shutil
import os
import subprocess
import json
#from pyproct.clustering.clustering import Clustering
import math
import atexit
import argparse
import time
import glob
import pickle
import blockNames
import atomset
import utilities
import multiprocessing
from functools import partial

from spawning import spawning, spawningTypes
from simulation import simulationrunner, simulationTypes
from clustering import clustering

def copyInitialStructures(initialStructures, tmpInitialStructuresTemplate, iteration):
    for i, name in enumerate(initialStructures):
        shutil.copyfile(name, tmpInitialStructuresTemplate%(iteration,i))

def createMultipleComplexesFilenames(numberOfSnapshots, inputFileTemplate, tmpInitialStructuresTemplate, iteration):
    jsonString = "\n"
    for i in range(numberOfSnapshots-1):
        jsonString += inputFileTemplate%(tmpInitialStructuresTemplate%(iteration,i)) + ",\n"
    jsonString += inputFileTemplate%(tmpInitialStructuresTemplate%(iteration,numberOfSnapshots-1))
    return jsonString

def generateSnapshotSelectionStringLastRound(currentEpoch, outputPathTempletized, trajectoryWildcard):
    return " \"" + os.path.join(outputPathTempletized%currentEpoch, trajectoryWildcard) + "\""

def countEpochAcceptedSteps(folder):
    acceptedSteps = os.popen("tail -qn 1 " + os.path.join(folder,"*report_*") + " | awk '{a+=$3}END{print a}'").read()
    return acceptedSteps

def makeClusterRepresentativesInitialStructures(tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clustering, iteration):

    counts = 0
    for i, cluster in enumerate(clustering.clusters.clusters):
        for j in range(int(degeneracyOfRepresentatives[i])):
            outputFilename = tmpInitialStructuresTemplate%(iteration,counts)
            print 'Writing to ', outputFilename, 'cluster', i
            cluster.writePDB(outputFilename)

            counts += 1

    print "counts & cluster centers", counts, np.where(np.array(degeneracyOfRepresentatives) > 0)[0].size
    return counts

def findFirstRun(outputPath, CLUSTERING_OUTPUT_OBJECT):
    """ Assumes that the outputPath is XXX/%d """

    folderWithSimulationData = outputPath
    allFolders = os.listdir(folderWithSimulationData)
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort(reverse=True)

    for epoch in epochFolders:
        if os.path.exists(CLUSTERING_OUTPUT_OBJECT%epoch):
            return epoch
    if epoch <= 0:
        return 0

def loadParams(jsonParams):
    """
        TODO: change for variables in a block names file, and work it out a bit more
    """
    jsonFile = open(jsonParams, 'r').read()
    parsedJSON = json.loads(jsonFile)

    return parsedJSON[blockNames.ControlFileParams.generalParams], parsedJSON[blockNames.ControlFileParams.spawningBlockname],\
            parsedJSON[blockNames.ControlFileParams.simulationBlockname], parsedJSON[blockNames.ControlFileParams.clusteringBlockname]

def saveInitialControlFile(jsonParams, originalControlFile):
    file = open(originalControlFile, 'w')
    jsonFile = open(jsonParams, 'r').read()
    file.write(jsonFile)

def main(jsonParams=None):
    if jsonParams is None:
        jsonParams = sys.argv[1]

    generalParams, spawningBlock, simulationrunnerBlock, clusteringBlock = loadParams(jsonParams)

    spawningAlgorithmBuilder = spawning.SpawningAlgorithmBuilder()
    startingConformationsCalculator, spawningParams = spawningAlgorithmBuilder.build(spawningBlock)

    runnerbuilder = simulationrunner.RunnerBuilder()
    simulationRunner = runnerbuilder.build(simulationrunnerBlock)

    clusteringType = clusteringBlock[blockNames.ClusteringTypes.type]

    restart = generalParams[blockNames.GeneralParams.restart]
    debug = generalParams[blockNames.GeneralParams.debug]
    outputPath = generalParams[blockNames.GeneralParams.outputPath]
    initialStructures = generalParams[blockNames.GeneralParams.initialStructures]

    print "================================"
    print "            PARAMS              "
    print "================================"
    print "Restarting simulations", generalParams[blockNames.GeneralParams.restart]
    print "Debug:", generalParams[blockNames.GeneralParams.debug]

    print "Iterations: %d, Mpi processors: %d, Pele steps: %d"%(simulationRunner.parameters.iterations, simulationRunner.parameters.processors, simulationRunner.parameters.peleSteps)

    print "SpawningType:", spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[startingConformationsCalculator.type]

    print "SimulationType:", simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY[simulationRunner.type]
    print "Clustering method:", clusteringType

    print "Output path: ", outputPath
    print "Initial Structures: ", initialStructures
    print "================================\n\n"


    #PRIVATE CONSTANTS
    ORIGINAL_CONTROLFILE = os.path.join(outputPath, "originalControlFile.conf")

    outputPathTempletized = os.path.join(outputPath, "%d")

    tmpFolder = "tmp_" +  outputPath.replace("/", "_")

    inputFileTemplate = "{ \"files\" : [ { \"path\" : \"%s\" } ] }"
    tmpInitialStructuresTemplate = tmpFolder+"/initial_%d_%d.pdb"
    tmpControlFilename = tmpFolder+"/controlFile%d.conf"
    tmpControlFilename_pyproct = tmpFolder+"/pyproct_controlFile%d.conf"
    tmpControlFilename_pyproct_thisEpoch = tmpFolder+"/pyproct_thisEpoch_controlFile%d.conf"
    trajectoryBasename = "*traj*"
    ligandTrajectoryBasename = "ligandSnapshots"

    trajectoryWildcard = trajectoryBasename + "_%d.pdb"
    ligandTrajectory = ligandTrajectoryBasename + "_%s.pdb" #e.g. ligandSnapshots_%s.pdb



    PYPROCT = "pyproct.main"
    PYPROCT_OUTPUT_DIR = outputPathTempletized+"/pyproct"
    PYPROCT_OUTPUT_DIR_SINGLE_EPOCH = outputPathTempletized+"/pyproct_thisEpoch"
    PYPROCT_REPRESENTATIVE_OUTPUT = PYPROCT_OUTPUT_DIR + "/results/representatives.pdb"
    PYPROCT_RESULTS_OUTPUT = PYPROCT_OUTPUT_DIR + "/results/results.json"
    PYPROCT_COMPRESSION_OUTPUT = PYPROCT_OUTPUT_DIR_SINGLE_EPOCH + "/results/compressed.pdb"

    PYPROCT_COMPRESSION_BLOCK = ",\"compression\":{\n\
            \"keep_remarks\": \"ALL\",\n\
            \"add_source_details\": true,\n\
            \"final_number_of_frames\": %d,\n\
            \"file\": \"compressed.pdb\",\n\
            \"type\":\"KMEDOIDS\"\n\
        }"#note leading comma

    CLUSTERING_OUTPUT_DIR = os.path.join(outputPathTempletized, "clustering")
    CLUSTERING_OUTPUT_OBJECT = os.path.join(CLUSTERING_OUTPUT_DIR, "object.pkl")
    #END PRIVATE CONSTANTS


    #if not debug: atexit.register(utilities.cleanup, tmpFolder)

    utilities.makeFolder(outputPath)
    utilities.makeFolder(tmpFolder)
    saveInitialControlFile(jsonParams, ORIGINAL_CONTROLFILE)
    
    #print variable epsilon information
    if not spawningParams.epsilon is None:
        epsilon_file = open("epsilon_values.txt","w")
        epsilon_file.write("Iteration\tEpsilon\n")
        epsilon_file.close()

    if restart:
        firstRun = findFirstRun(outputPath, CLUSTERING_OUTPUT_OBJECT)

        if firstRun != 0:
            with open(CLUSTERING_OUTPUT_OBJECT%(firstRun), 'rb') as f:
                clusteringMethod = pickle.load(f)

            degeneracyOfRepresentatives = startingConformationsCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, firstRun)
            startingConformationsCalculator.log()
            print "Degeneracy", degeneracyOfRepresentatives

            seedingPoints = makeClusterRepresentativesInitialStructures(tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clusteringMethod, firstRun)

            initialStructuresAsString = createMultipleComplexesFilenames(seedingPoints, inputFileTemplate, tmpInitialStructuresTemplate, firstRun)
            
    if not restart or firstRun == 0: # if restart and firstRun = 0, it must check into the initial structures
        #Choose initial structures
        if not debug: shutil.rmtree(outputPath)
        utilities.makeFolder(outputPath)
        firstRun = 0
        saveInitialControlFile(jsonParams, ORIGINAL_CONTROLFILE)
        copyInitialStructures(initialStructures, tmpInitialStructuresTemplate, firstRun)
        initialStructuresAsString = createMultipleComplexesFilenames(len(initialStructures), inputFileTemplate, tmpInitialStructuresTemplate, firstRun)


    peleControlFileDictionary = {"COMPLEXES":initialStructuresAsString, "PELE_STEPS":simulationRunner.parameters.peleSteps}

    #Make control file
    outputDir = outputPathTempletized%firstRun
    utilities.makeFolder(outputDir)
    peleControlFileDictionary["OUTPUT_PATH"] = outputDir
    peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + firstRun*simulationRunner.parameters.processors
    #simulationRunner.parameters needn't be passed to makeWorkingControlFile
    simulationRunner.makeWorkingControlFile(simulationRunner.parameters.templetizedControlFile, tmpControlFilename%firstRun, peleControlFileDictionary) 

    for i in range(firstRun, simulationRunner.parameters.iterations):
        print "Iteration", i
        if not spawningParams.epsilon is None:
            epsilon_file = open("epsilon_values.txt","a")
            epsilon_file.write("%d\t%f\n"%(i,spawningParams.epsilon))
            epsilon_file.close()
        
        print "Production run..."
        if not debug:
            startTime = time.time() 
            simulationRunner.runSimulation(tmpControlFilename%i)
            endTime = time.time() 
            print "PELE %s sec" % (endTime - startTime)


        print "Clustering..."
        startTime = time.time() 
        #snapshotsJSONSelectionString = generateSnapshotSelectionString(i, outputPathTempletized, trajectoryBasename)
        snapshotsJSONSelectionString = generateSnapshotSelectionStringLastRound(i, outputPathTempletized, trajectoryBasename)
        snapshotsJSONSelectionString = "[" + snapshotsJSONSelectionString + "]"
        paths = eval(snapshotsJSONSelectionString)
        if len(glob.glob(paths[-1])) == 0: sys.exit("No more trajectories to cluster")
        if i == 0:
            clusteringBuilder = clustering.ClusteringBuilder()
            clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                                 spawningParams.reportFilename,
                                                                 spawningParams.reportCol)
        #else:
        #    #CAN'T THIS BE REMOVED????
        #    with open(CLUSTERING_OUTPUT_OBJECT%(i-1), 'rb') as input:
        #        clusteringMethod = pickle.load(input)
        clusteringMethod.cluster(paths)
        endTime = time.time() 
        print "Clustering ligand: %s sec" % (endTime - startTime)


        degeneracyOfRepresentatives = startingConformationsCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, i)
        startingConformationsCalculator.log()
        print "Degeneracy", degeneracyOfRepresentatives

        numberOfSeedingPoints = makeClusterRepresentativesInitialStructures(tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clusteringMethod, i+1)

        clusteringMethod.writeOutput(CLUSTERING_OUTPUT_DIR%i, degeneracyOfRepresentatives, CLUSTERING_OUTPUT_OBJECT%i)
        if i > 0:
            # Remove old clustering object, since we already have a newer one
            os.remove(CLUSTERING_OUTPUT_OBJECT%(i-1))

        #Prepare for next pele iteration
        if i != simulationRunner.parameters.iterations-1:
            initialStructuresAsString = createMultipleComplexesFilenames(numberOfSeedingPoints, inputFileTemplate, tmpInitialStructuresTemplate, i+1)
            peleControlFileDictionary["COMPLEXES"] = initialStructuresAsString

            outputDir = outputPathTempletized%(i+1)
            utilities.makeFolder(outputDir) #PELE does not do it automatically
            peleControlFileDictionary["OUTPUT_PATH"] = outputDir
            peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + (i+1)*simulationRunner.parameters.processors
            #simulationRunner.parameters needn't be passed to makeWorkingControlFile
            simulationRunner.makeWorkingControlFile(simulationRunner.parameters.templetizedControlFile, tmpControlFilename%(i+1), peleControlFileDictionary) 

    #utilities.cleanup
    #utilities.cleanup(tmpFolder)

if __name__ == '__main__':
    main()
