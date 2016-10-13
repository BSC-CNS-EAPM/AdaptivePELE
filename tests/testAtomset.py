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
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb")

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom1 = atomset.Atom("ATOM      1  N   ASN A   1       7.920  22.268   9.257 1.00 15.18           N1+")
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom3 = atomset.Atom("ATOM      3  CA  ASN A   1       7.870  19.937   8.524 1.00 16.63           C  ")
        atom4 = atomset.Atom("ATOM      4  O   ASN A   1       7.030  20.308   7.700 1.00 16.13           O  ")
        atom5 = atomset.Atom("ATOM      5  CB  CYS A   1       8.108  20.445  11.030 1.00 16.53           C  ")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2, atom3.id: atom3, atom4.id: atom4, atom5.id: atom5}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_from_str(self):
        # preparation
        pdbfile = open("tests/data/pdb_test.pdb", "r")
        pdbstring = pdbfile.read()
        pdbfile.close()
        pdb = atomset.PDB()

        # function to test
        pdb.initialise(pdbstring)

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom1 = atomset.Atom("ATOM      1  N   ASN A   1       7.920  22.268   9.257 1.00 15.18           N1+")
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom3 = atomset.Atom("ATOM      3  CA  ASN A   1       7.870  19.937   8.524 1.00 16.63           C  ")
        atom4 = atomset.Atom("ATOM      4  O   ASN A   1       7.030  20.308   7.700 1.00 16.13           O  ")
        atom5 = atomset.Atom("ATOM      5  CB  CYS A   1       8.108  20.445  11.030 1.00 16.53           C  ")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2, atom3.id: atom3, atom4.id: atom4, atom5.id: atom5}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_resname(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb", resname="CYS")

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom5 = atomset.Atom("ATOM      5  CB  CYS A   1       8.108  20.445  11.030 1.00 16.53           C  ")
        goldenAtomsDict = {atom2.id: atom2, atom5.id: atom5}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_atomname(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb", atomname="CA")

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom3 = atomset.Atom("ATOM      3  CA  ASN A   1       7.870  19.937   8.524 1.00 16.63           C  ")
        goldenAtomsDict = {atom2.id: atom2, atom3.id: atom3}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_type_protein(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test_ligand.pdb",
                       type=atomset.PDB.typeProtein)

        # assertion
        pdbContent = "MODEL        1\n\
ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  \n\
ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  \n\
TER\n\
HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  \n\
TER\n\
HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+\n\
TER\n\
ENDMDL\n\
END   \n"

        atom1 = atomset.Atom("ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  ")
        goldenAtomsDict = {atom1.id: atom1}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_type_hetero(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test_ligand.pdb", type=atomset.PDB.typeHetero)

        # assertion
        pdbContent = "MODEL        1\n\
ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  \n\
ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  \n\
TER\n\
HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  \n\
TER\n\
HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+\n\
TER\n\
ENDMDL\n\
END   \n"

        atom1 = atomset.Atom("HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  ")
        atom2 = atomset.Atom("HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_type_heavyAtoms(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test_ligand.pdb", heavyAtoms=False)

        # assertion
        pdbContent = "MODEL        1\n\
ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  \n\
ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  \n\
TER\n\
HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  \n\
TER\n\
HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+\n\
TER\n\
ENDMDL\n\
END   \n"

        atom1 = atomset.Atom("ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  ")
        atom2 = atomset.Atom("ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  ")
        atom3 = atomset.Atom("HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  ")
        atom4 = atomset.Atom("HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2, atom3.id: atom3, atom4.id: atom4}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_COM(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb")

        # assertion
        total_mass = 66.0382
        COM_array = np.array([516.1336264, 1373.048894, 602.7150822])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass

        np.testing.assert_array_almost_equal(COM, pdb.getCOM())
        self.assertAlmostEqual(total_mass, pdb.totalMass)

    def testPDB_write(self):
        # preparation
        pdb = atomset.PDB()
        pdb.initialise("tests/data/pdb_test.pdb")

        # function to test
        pdb.writePDB("tests/data/pdb_test_write.pdb")

        # assertion
        pdbtestfile = open("tests/data/pdb_test_write.pdb", "r")
        pdbtestsstr = pdbtestfile.read()
        pdbtestfile.close()
        self.assertEqual(pdb.pdb, pdbtestsstr)

    def testPDB_RMSD(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/ain_native_fixed.pdb", resname='AIN')
        pdb_traj = atomset.PDB()
        pdb_traj.initialise("tests/data/ain_trajectory.pdb", resname='AIN')

        # function to test
        RMSD = atomset.computeRMSD(pdb_native, pdb_traj)
        golden_RMSD = 3.928617
        self.assertAlmostEqual(RMSD, golden_RMSD, 5)

    def testPDB_RMSD_symmetries(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/ain_native_fixed.pdb", resname='AIN')
        pdb_traj = atomset.PDB()
        pdb_traj.initialise("tests/data/ain_trajectory.pdb", resname='AIN')

        # function to test
        RMSD = atomset.computeRMSD(pdb_native, pdb_traj, {"1733:O1:AIN":"1735:O2:AIN", "1735:O2:AIN":"1733:O1:AIN"})
        golden_RMSD = 3.860743
        self.assertAlmostEqual(RMSD, golden_RMSD, 5)

    def testPDB_contacts(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/native_ain.pdb")

        # function to test
        contacts = pdb_native.countContacts("AIN", 8)
        golden_contacts = 19
        self.assertEqual(contacts, golden_contacts)

    def testPDB_contactmap(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/pdb_test_contact.pdb")

        # function to test
        contact_map, contacts = pdb_native.createContactMap("AIN", 8)
        golden_contact_map = np.array([[1, 0, 0, 0], [0, 1, 1, 1]])
        golden_contacts = pdb_native.countContacts("AIN", 8)
        np.testing.assert_array_equal(contact_map, golden_contact_map)
        self.assertEqual(golden_contacts, contacts)
