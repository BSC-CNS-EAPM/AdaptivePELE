import unittest
import tests.testSpawning as tSpawning
import tests.testAtomset as tAtomset
import tests.testClustering as tClustering
import tests.testAdaptiveSampling as tAdaptive
import tests.testThresholdcalculator as tThreshold
import tests.testDensityCalculator as tDensity

def main():
    testSuite = unittest.TestSuite()

    # testSuite.addTest(unittest.makeSuite(tSpawning.TestSpawningCalculator))
    testSuite.addTest(unittest.makeSuite(tAtomset.atomsetTest))
    # testSuite.addTest(unittest.makeSuite(tClustering.clusteringTest))
    # testSuite.addTest(unittest.makeSuite(tAdaptive.TestadaptiveSampling))
    # testSuite.addTest(unittest.makeSuite(tThreshold.thresholdCalculatorTest))
    # testSuite.addTest(unittest.makeSuite(tDensity.densityCalculatorTest))

    runner = unittest.TextTestRunner()
    runner.run(testSuite)

if __name__ == "__main__":
    main()
