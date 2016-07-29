import atomset
import unittest
import numpy as np

class atomsetTest(unittest.TestCase):
    """ For the moment the tests include loading from file and string, resname
    and atomname selection. It uses a toy pdb of only 5 lines located in
    tests/data
    TODO: heavyAtoms and type selection
    """
    def testPDB_from_file(self):
        pdb = atomset.PDB()
        pdb.initialise("tests/data/pdb_test.pdb")
        total_mass = 66.0382
        COM_array = np.array([516.1336264, 1373.048894, 602.7150822])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass
        np.testing.assert_array_almost_equal(COM, pdb.getCOM())
        self.assertAlmostEqual(total_mass, pdb.totalMass)
    def testPDB_from_str(self):
        pdbfile = open("tests/data/pdb_test.pdb", "r")
        pdbstring = pdbfile.read()
        pdbfile.close()
        pdb = atomset.PDB()
        pdb.initialise(pdbstring)
        total_mass = 66.0382
        COM_array = np.array([516.1336264, 1373.048894, 602.7150822])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass
        np.testing.assert_array_almost_equal(COM, pdb.getCOM())
        self.assertAlmostEqual(total_mass, pdb.totalMass)
    def testPDB_from_file_sel_resname(self):
        pdb = atomset.PDB()
        pdb.initialise("tests/data/pdb_test.pdb", resname="CYS")
        total_mass = 24.0214
        COM_array = np.array([198.2005714, 496.7745627, 247.4804735])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass
        self.assertAlmostEqual(total_mass, pdb.totalMass)
        np.testing.assert_array_almost_equal(COM, pdb.getCOM())

    def testPDB_from_str_sel_resname(self):
        pdbfile = open("tests/data/pdb_test.pdb", "r")
        pdbstring = pdbfile.read()
        pdbfile.close()
        pdb = atomset.PDB()
        pdb.initialise(pdbstring, resname="CYS")
        total_mass = 24.0214
        COM_array = np.array([198.2005714, 496.7745627, 247.4804735])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass
        self.assertAlmostEqual(total_mass, pdb.totalMass)
        np.testing.assert_array_almost_equal(COM, pdb.getCOM())
    def testPDB_from_file_sel_atomname(self):
        pdb = atomset.PDB()
        pdb.initialise("tests/data/pdb_test.pdb", atomname="CA")
        total_mass = 24.0214
        COM_array = np.array([195.3420248, 490.6731271, 217.3816593])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass
        self.assertAlmostEqual(total_mass, pdb.totalMass)
        np.testing.assert_array_almost_equal(COM, pdb.getCOM())

    def testPDB_from_str_sel_atomname(self):
        pdbfile = open("tests/data/pdb_test.pdb", "r")
        pdbstring = pdbfile.read()
        pdbfile.close()
        pdb = atomset.PDB()
        pdb.initialise(pdbstring, atomname="CA")
        total_mass = 24.0214
        COM_array = np.array([195.3420248, 490.6731271, 217.3816593])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass
        self.assertAlmostEqual(total_mass, pdb.totalMass)
        np.testing.assert_array_almost_equal(COM, pdb.getCOM())

        
