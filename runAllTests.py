import unittest
import tests.testStartingConformationsCalculator as tStartConf

def main():
    testSuite = unittest.TestSuite()

    testSuite.addTest(unittest.makeSuite(tStartConf.TestStartingConformationsCalculator))

    runner = unittest.TextTestRunner()
    runner.run(testSuite)

if __name__ == "__main__":
    main()
