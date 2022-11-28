from __future__ import absolute_import, division, print_function, unicode_literals
import unittest

from AdaptivePELE.simulation import simulationrunner

class TestSimulationRunner(unittest.TestCase):
    def testMultipleControlFileSelection(self):
        params = simulationrunner.SimulationParameters()
        params.templetizedControlFile = ["file1.txt", "file2.txt", "file3.txt"]
        runner = simulationrunner.PeleSimulation(params)
        files_out = [runner.getControlFileForEpoch(i) for i in range(10)]
        golden = ["file1.txt", "file2.txt", "file3.txt"]*3+["file1.txt"]
        self.assertEqual(files_out, golden)

    def testOneControlFileSelection(self):
        params = simulationrunner.SimulationParameters()
        params.templetizedControlFile = ["file1.txt"]
        runner = simulationrunner.PeleSimulation(params)
        files_out = [runner.getControlFileForEpoch(i) for i in range(10)]
        golden = ["file1.txt"]*10
        self.assertEqual(files_out, golden)

    def testMultipleControlFileSelectionEven(self):
        params = simulationrunner.SimulationParameters()
        params.templetizedControlFile = ["file1.txt", "file2.txt"]
        runner = simulationrunner.PeleSimulation(params)
        files_out = [runner.getControlFileForEpoch(i) for i in range(10)]
        golden = ["file1.txt", "file2.txt"]*5
        self.assertEqual(files_out, golden)

    def testMultipleControlFileSelectionLessThan(self):
        params = simulationrunner.SimulationParameters()
        params.templetizedControlFile = [f"file{i}.txt" for i in range(10)]
        runner = simulationrunner.PeleSimulation(params)
        files_out = [runner.getControlFileForEpoch(i) for i in range(5)]
        golden = [f"file{i}.txt" for i in range(5)]
        self.assertEqual(files_out, golden)
