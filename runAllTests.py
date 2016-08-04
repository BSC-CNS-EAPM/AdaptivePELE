import unittest
import tests.testStartingConformationsCalculator as tStartConf
import tests.testatomset as tatomset
import tests.testClustering as tClustering
import tests.testadaptivesampling as tadaptive

def main():
    testSuite = unittest.TestSuite()

    testSuite.addTest(unittest.makeSuite(tStartConf.TestStartingConformationsCalculator))
    testSuite.addTest(unittest.makeSuite(tatomset.atomsetTest))
    testSuite.addTest(unittest.makeSuite(tClustering.clusteringTest))
    testSuite.addTest(unittest.makeSuite(tadaptive.TestadaptiveSampling))

    runner = unittest.TextTestRunner()
    runner.run(testSuite)

if __name__ == "__main__":
    main()
