from __future__ import absolute_import, division, print_function, unicode_literals
import unittest
import shutil
import os
import AdaptivePELE.adaptiveSampling as adaptiveSampling


class TestMD(unittest.TestCase):

    def check_succesful_simulation(self, output, epochs):
        for epoch in range(epochs):
            self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "summary.txt")))
        self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "object.pkl")))

    def testOpenMM3ptb(self):
        output_path = "tests/data/openmm_3ptb"
        controlFile = "tests/data/templetized_controlFile_3ptb_md.conf"

        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2)
        # cleanup
        shutil.rmtree(output_path)

    def testOpenMM1ab1(self):
        output_path = "tests/data/openmm_1ab1"
        controlFile = "tests/data/templetized_controlFile_1ab1_md.conf"
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2)
        # cleanup
        shutil.rmtree(output_path)
