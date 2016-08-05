import unittest
import tests.testSpawning as tSpawning
import tests.testAtomset as tAtomset
import tests.testClustering as tClustering
import tests.testAdaptiveSampling as tAdaptive

def main():
    testSuite = unittest.TestSuite()

    testSuite.addTest(unittest.makeSuite(tSpawning.TestSpawningCalculator))
    testSuite.addTest(unittest.makeSuite(tAtomset.atomsetTest))
    testSuite.addTest(unittest.makeSuite(tClustering.clusteringTest))
    testSuite.addTest(unittest.makeSuite(tAdaptive.TestadaptiveSampling))

    runner = unittest.TextTestRunner()
    runner.run(testSuite)

if __name__ == "__main__":
    main()
