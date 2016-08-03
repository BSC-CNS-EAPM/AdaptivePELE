import adaptiveSampling
import unittest
import glob
import pickle

class TestadaptiveSampling(unittest.TestCase):
    def test_integration(self):
        adaptiveSampling.main("tests/data/3ptb_data/integrationTest1.conf")
        golden_path="tests/data/3ptb_data/originTest1/"
        output_path="tests/data/3ptb_data/Test1/"
        for i in range(3):
            with open(golden_path+i+"/object.pkl",'rb') as f:
                golden_cluster = pickle.load(f)
            with open(output_path+i+"/object.pkl",'rb') as f2:
                output_cluster = pickle.load(f2)
            self.assertEqual(output_cluster, golden_cluster)

