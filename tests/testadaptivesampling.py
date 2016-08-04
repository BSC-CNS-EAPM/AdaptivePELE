import adaptiveSampling
import unittest
import glob
import pickle
import atomset
import shutil
import os

class TestadaptiveSampling(unittest.TestCase):
    def integrationTest(self, controlFile, goldenPath, outputPath):
        #Function to test --> integration test
        adaptiveSampling.main(controlFile)

        #Assertions

        #check clustering objects
        goldenPathObject = os.path.join(goldenPath, "%d/clustering/object.pkl")
        outputPathObject = os.path.join(outputPath, "%d/clustering/object.pkl")
        for i in range(3):
            with open(goldenPathObject%i,'rb') as f:
                goldenCluster = pickle.load(f)
            with open(outputPathObject%i,'rb') as f2:
                outputCluster = pickle.load(f2)
            self.assertEqual(outputCluster, goldenCluster)


        #check initial structures
        goldenPathInitial=os.path.join(goldenPath, "%d/initial_%d.pdb")
        tmpFolder = "tmp_" +  outputPath.replace("/", "_")
        outputPathInitial=os.path.join(tmpFolder, "initial_%d_%d.pdb")

        j = 0
        for j in range(3):
            for ij in range(1,5):
                goldenInitial = atomset.PDB()
                goldenInitial.initialise(goldenPathInitial%(j,ij))

                outputInitial = atomset.PDB()
                if j == 0:
                    #initial structures will always be "initial_0_0.pdb"
                    outputInitial.initialise(outputPathInitial%(j,0))
                else:
                    outputInitial.initialise(outputPathInitial%(j,ij%4))

                self.assertEqual(goldenInitial, outputInitial)

        #cleanup
        shutil.rmtree(outputPath)
        shutil.rmtree(tmpFolder)

    def testIntegration1(self):
        controlFile = "tests/data/3ptb_data/integrationTest1.conf"
        goldenPath="tests/data/3ptb_data/originTest1"
        outputPath="tests/data/3ptb_data/Test1"

        self.integrationTest(controlFile, goldenPath, outputPath)

    def testIntegration2(self):
        controlFile = "tests/data/3ptb_data/integrationTest2.conf"
        goldenPath="tests/data/3ptb_data/srcTest2Epsilon"
        outputPath="tests/data/3ptb_data/Test2"

        self.integrationTest(controlFile, goldenPath, outputPath)
