import adaptiveSampling
import unittest
import glob
import pickle
import atomset
import shutil

class TestadaptiveSampling(unittest.TestCase):
    def test_integration(self):
        adaptiveSampling.main("tests/data/3ptb_data/integrationTest1.conf")
        golden_path="tests/data/3ptb_data/originTest1/%d/clustering/object.pkl"
        output_path="tests/data/3ptb_data/Test1/%d/clustering/object.pkl"
        for i in range(3):
            with open(golden_path%i,'rb') as f:
                golden_cluster = pickle.load(f)
            with open(output_path%i,'rb') as f2:
                output_cluster = pickle.load(f2)
            self.assertEqual(output_cluster, golden_cluster)
        golden_path_initial="tests/data/3ptb_data/originTest1/%d/initial_%d.pdb"
        output_path_initial="tmp_tests_data_3ptb_data_Test1/initial_%d_%d.pdb"

        j = 0
        for ij in range(1,5):
            golden_initial = atomset.PDB()
            golden_initial.initialise(golden_path_initial%(j,ij))
            output_initial = atomset.PDB()
            output_initial.initialise(output_path_initial%(j,0))
            self.assertEqual(golden_initial, output_initial)
        for j in range(1,3):
            for ij in range(1,5):
                golden_initial = atomset.PDB()
                golden_initial.initialise(golden_path_initial%(j,ij))
                output_initial = atomset.PDB()
                output_initial.initialise(output_path_initial%(j,ij%4))
                self.assertEqual(golden_initial, output_initial)
        shutil.rmtree("tests/data/3ptb_data/Test1/")
        shutil.rmtree("tmp_tests_data_3ptb_data_Test1/")

