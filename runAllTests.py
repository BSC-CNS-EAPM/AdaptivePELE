import unittest
import testStartingConformationsCalculator

def main():
    testSuite = unittest.TestSuite()

    testSuite.addTest(unittest.makeSuite(testStartingConformationsCalculator.TestStartingConformationsCalculator))

    runner = unittest.TextTestRunner()
    runner.run(testSuite)

if __name__ == "__main__":
    main()
