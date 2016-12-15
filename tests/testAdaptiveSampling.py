import adaptiveSampling
import unittest
import pickle
from atomset import atomset
import shutil
import os


class TestadaptiveSampling(unittest.TestCase):

    def checkClusteringObjects(self, goldenPath, outputPath):
        goldenPathObject = os.path.join(goldenPath, "%d/clustering/object.pkl")
        outputPathObject = os.path.join(outputPath, "%d/clustering/object.pkl")
        for i in range(2, 3):
            with open(goldenPathObject % i, 'rb') as f:
                goldenCluster = pickle.load(f)
            with open(outputPathObject % i, 'rb') as f2:
                outputCluster = pickle.load(f2)
            goldenCluster.clusters.printClusters()
            outputCluster.clusters.printClusters()
            self.assertEqual(outputCluster, goldenCluster)

    def checkStartingConformations(self, goldenPath, outputPath):
        goldenPathInitial = os.path.join(goldenPath, "%d/initial_%d.pdb")
        outputPathInitial = os.path.join(outputPath, "initial_%d_%d.pdb")

        j = 0
        for j in range(3):
            for ij in range(1, 5):
                goldenInitial = atomset.PDB()
                goldenInitial.initialise(goldenPathInitial % (j, ij))

                outputInitial = atomset.PDB()
                if j == 0:
                    # initial structures will always be "initial_0_0.pdb"
                    outputInitial.initialise(outputPathInitial % (j, 0))
                else:
                    outputInitial.initialise(outputPathInitial % (j, ij % 4))

                self.assertEqual(goldenInitial, outputInitial)

    def checkTrajectories(self, goldenPath, outputPath):
        goldenPathTrajectory = os.path.join(goldenPath, "%d/trajectory_%d.pdb")
        outputPathTrajectory = os.path.join(outputPath, "%d/trajectory_%d.pdb")

        for epoch in range(3):
            for i in range(1, 5):
                goldenTrajFile = open(goldenPathTrajectory % (epoch, i), 'r')
                goldenTraj = goldenTrajFile.read()

                goldenTrajFile.close()
                outputTrajFile = open(outputPathTrajectory % (epoch, i), 'r')
                outputTraj = outputTrajFile.read()
                outputTrajFile.close()

                self.assertEqual(goldenTraj, outputTraj)

    def integrationTest(self, controlFile, goldenPath, outputPath):
        # Function to test --> integration test
        adaptiveSampling.main(controlFile)

        # Assertions
        tmpFolder = "tmp_" + outputPath.replace("/", "_")

        self.checkClusteringObjects(goldenPath, outputPath)
        self.checkStartingConformations(goldenPath, tmpFolder)
        self.checkTrajectories(goldenPath, outputPath)

        # cleanup
        shutil.rmtree(outputPath)
        shutil.rmtree(tmpFolder)

    def testIntegration1(self):
        """
            Simulations are not run, trajectories and reports are copied
        """
        controlFile = "tests/data/3ptb_data/integrationTest1.conf"
        goldenPath = "tests/data/3ptb_data/originTest1"
        outputPath = "tests/data/3ptb_data/Test1"

        self.integrationTest(controlFile, goldenPath, outputPath)

    def testIntegration2(self):
        """
            Simulations are not run, trajectories and reports are copied
        """
        controlFile = "tests/data/3ptb_data/integrationTest2.conf"
        goldenPath = "tests/data/3ptb_data/srcTest2Epsilon"
        outputPath = "tests/data/3ptb_data/Test2"

        self.integrationTest(controlFile, goldenPath, outputPath)

    def testIntegration3(self):
        """
            Simulations are actually run
        """
        controlFile = "tests/data/3ptb_data/integrationTest3.conf"
        goldenPath = "tests/data/3ptb_data/originTest3"
        outputPath = "tests/data/3ptb_data/Test3"

        try:
            self.integrationTest(controlFile, goldenPath, outputPath)
        except SystemExit:
            # Catch error for not having PELE installed
            print ("Warning! There was a sysExit in test3, this is usually "
                   "caused by not having PELE installed, so it can be ignored "
                   "if the test are not running on MareNostrum or life")
            shutil.rmtree(outputPath)
            shutil.rmtree("tmp_tests_data_3ptb_data_Test3/")
