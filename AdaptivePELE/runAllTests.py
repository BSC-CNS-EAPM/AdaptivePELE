import unittest
from AdaptivePELE.tests import testSpawning as tSpawning
from AdaptivePELE.tests import testAtomset as tAtomset
from AdaptivePELE.tests import testClustering as tClustering
from AdaptivePELE.tests import testAdaptiveSampling as tAdaptive
from AdaptivePELE.tests import testThresholdcalculator as tThreshold
from AdaptivePELE.tests import testDensityCalculator as tDensity


def main():
    testSuite = unittest.TestSuite()

    testSuite.addTest(unittest.makeSuite(tAtomset.atomsetTest))
    testSuite.addTest(unittest.makeSuite(tSpawning.TestSpawningCalculator))
    testSuite.addTest(unittest.makeSuite(tThreshold.thresholdCalculatorTest))
    testSuite.addTest(unittest.makeSuite(tDensity.densityCalculatorTest))
    testSuite.addTest(unittest.makeSuite(tClustering.clusteringTest))
    testSuite.addTest(unittest.makeSuite(tAdaptive.TestadaptiveSampling))

    runner = unittest.TextTestRunner()
    runner.run(testSuite)

if __name__ == "__main__":
    main()
