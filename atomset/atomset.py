import numpy as np
import re
import StringIO
from scipy import sparse


class Atom:
    _chargePattern = re.compile("[0-9]|\+|\-")
    _ATOM_WEIGHTS = {"H": 1.00794,
                    "D": 2.01410178,  # deuterium
                    "HE": 4.00,
                    "LI": 6.941,
                    "BE": 9.01,
                    "B": 10.811,
                    "C": 12.0107,
                    "N": 14.0067,
                    "O": 15.9994,
                    "F": 18.998403,
                    "NE": 20.18,
                    "NA": 22.989769,
                    "MG": 24.305,
                    "AL": 26.98,
                    "SI": 28.09,
                    "P": 30.973762,
                    "S": 32.065,
                    "CL": 35.453,
                    "AR": 39.95,
                    "K": 39.0983,
                    "CA": 40.078,
                    "SC": 44.96,
                    "TI": 47.87,
                    "V": 50.94,
                    "CR": 51.9961,
                    "MN": 54.938045,
                    "FE": 55.845,
                    "CO": 58.93,
                    "NI": 58.6934,
                    "CU": 63.546,
                    "ZN": 65.409,
                    "GA": 69.72,
                    "GE": 72.64,
                    "AS": 74.9216,
                    "SE": 78.96,
                    "BR": 79.90,
                    "KR": 83.80,
                    "RB": 85.47,
                    "SR": 87.62,
                    "Y": 88.91,
                    "ZR": 91.22,
                    "NB": 92.91,
                    "W": 95.94,  # Molybdenum?  Not sure why it's not always MO
                    "MO": 95.94,
                    "TC": 98.0,
                    "RU": 101.07,
                    "RH": 102.91,
                    "PD": 106.42,
                    "AG": 107.8682,
                    "CD": 112.411,
                    "IN": 114.82,
                    "SN": 118.71,
                    "SB": 121.76,
                    "TE": 127.60,
                    "I": 126.90447,
                    "XE": 131.29,
                    "CS": 132.91,
                    "BA": 137.33,
                    "PR": 140.91,
                    "EU": 151.96,
                    "GD": 157.25,
                    "TB": 158.93,
                    "IR": 192.22,
                    "PT": 195.084,
                    "AU": 196.96657,
                    "HG": 200.59,
                    "PB": 207.2,
                    "U": 238.03}

    def __init__(self, atomContent):
        """ Create an atom from a pdb line

            :param atomContent: Line of the pdb from which the atom will be created
            :type atomContent: str
        """
        # atomContent = atomContent.split()
        if len(atomContent) > 6 and (atomContent[:4] == 'ATOM' or atomContent[:6] == 'HETATM'):
            self.atomSerial = atomContent[6:11].strip()
            self.name = atomContent[12:16].strip()
            self.resname = atomContent[17:20].strip()
            self.resChain = atomContent[21]
            self.resnum = atomContent[22:26].strip()
            self.x = float(atomContent[30:38])
            self.y = float(atomContent[38:46])
            self.z = float(atomContent[46:54])

            self.type = re.sub(self._chargePattern, "", atomContent[76:80]).strip().upper()
            self.mass = self._ATOM_WEIGHTS[self.type]

            if atomContent.startswith('ATOM'):
                self.protein = True
            else:
                self.protein = False

            self.id = self.atomSerial + ":" + self.name + ":" + self.resname
            # self.id = self.atomSerial

    def isHeavyAtom(self):
        """
            Check if Atom is a heavy atom

            :returns: bool -- True if Atom is heavy atom, false otherwise
        """
        return self.type != 'H'

    def isProtein(self):
        """
            Check if Atom is a protein atom

            :returns: bool -- True if Atom is a protein atom, false otherwise
        """
        return self.protein

    def isHeteroAtom(self):
        """
            Check if Atom is an hetero atom

            :returns: bool -- True if Atom is an hetero atom, false otherwise
        """
        return not self.protein

    def printAtom(self):
        """
            Print Atom information
        """
        print self.atomSerial, self.name, self.resname, self.resChain, self.resnum, self.x, self.y, self.z, self.type, self.mass
        # print self.atomSerial, self.name, self.resname, self.resChain, self.resnum, self.r, self.type, self.mass

    def __eq__(self, atom2):
        return self.id == atom2.id

    def __lt__(self, atom2):
        return self.serial < atom2.serial

    def __str__(self):
        return "%s: %s %s %s [%f, %f, %f] %s %f" % (self.id, self.atomSerial,
                                                    self.resChain, self.resnum,
                                                    self.x, self.y, self.z,
                                                    self.type, self.mass)
        # return "%s: %s %s %s [%f, %f, %f] %s %f"%(self.id, self.atomSerial, self.resChain, self.resnum, self.r[0], self.r[1], self.r[2], self.type, self.mass)

    def getAtomCoords(self):
        """
            Get the coordinates of the atom

            :returns: numpy.Array -- Array with the coordinate of the atom
        """
        return np.array([self.x, self.y, self.z])

    def squaredDistance(self, atom2):
        """
            Calculate the squared distance between two atoms

            :param atom2: Second Atom to whom the distance will be calculated
            :type atom2: Atom
            :returns: float -- The distance between the atoms
        """
        d = (self.x - atom2.x)**2 + (self.y - atom2.y)**2 + (self.z - atom2.z)**2
        return d


class PDB:
    _typeProtein = "PROTEIN"
    _typeHetero = "HETERO"
    _typeAll = "ALL"

    def __init__(self):
        """
            Object that will contain the information of a PDB file. Has to call
            the initialise method to load the file
        """
        self.atoms = {}
        # {atomId: atom, ...}
        # Where atomId := serail:atomName:resName
        self.totalMass = 0
        self.pdb = ""
        self.com = None
        self.centroid = None

        # Necessary for contactMaps
        self.atomList = []

    def __eq__(self, other):
        """ Compare two pdb strings, remark lines should be ignored and only the
        atoms and its information should be compared"""
        pdb1 = [element.strip() for element in self.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
        pdb2 = [element.strip() for element in other.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
        return pdb1 == pdb2

    def initialise(self, PDBstr, heavyAtoms=True, resname="", atomname="", type="ALL"):
        """
            Load the information from a PDB file or a string with the PDB
            contents

            :param PDBstr: may be a path to the PDB file or a string with the contents of the PDB
            :type PDBstr: str
            :param heavyAtoms: wether to consider only heavy atoms (True if onl y heavy atoms have to be considered)
            :type heavyAtoms: bool
            :param resname: Residue name to select from the pdb (will only select the residues with that name)
            :type resname: str
            :param resname: Residue name to select from the pdb (will only select the residues with that name)
            :type resname: str
            :param type: type of atoms to select: may be ALL, PROTEIN or HETERO
            :type type: str
            :raises: ValueError if the pdb contained no atoms
        """
        PDBContent = StringIO.StringIO(readPDB(PDBstr))  # Using StringIO
        # creates a buffer that can handle a pdb file or a string containing
        # the PDB
        self.pdb = PDBContent.read()  # in case one wants to write it

        stringWithPDBContent = self.pdb.split('\n')

        for atomLine in stringWithPDBContent:

            # HUGE optimisation (not to create an atom each time ~ 2e-5 s/atom)
            if resname != "" and not atomLine[17:20].strip() == resname:
                # optimisation
                continue
            if atomname != "" and not atomLine[12:16].strip() == atomname:
                # optimisation
                continue

            atom = Atom(atomLine)
            # Here atom will be not null, empty or not.
            # With "try", we prune empty atoms
            try:
                if (not heavyAtoms or atom.isHeavyAtom()) and\
                   (type == self._typeAll or (type == self._typeProtein and atom.isProtein()) or (type == self._typeHetero and atom.isHeteroAtom())):
                        self.atoms.update({atom.id: atom})
                        self.atomList.append(atom.id)
            except:
                pass
        if self.atoms == {}:
            raise ValueError('The input pdb file/string was empty, no atoms loaded!')

    def computeTotalMass(self):
        """
            Calculate the total mass of the PDB
        """
        self.totalMass = 0
        for atomId, atom in self.atoms.items():
            self.totalMass += atom.mass

    def printAtoms(self):
        """
            Print Atom information for all the atoms in the PDB
        """
        for atom in self.atoms.values():
            print atom  # atom.printAtom()

    def getNumberOfAtoms(self):
        """
            Get the number of Atoms in the PDB

            :returns: int -- Number of atoms in the PDB
        """
        return len(self.atoms)

    def getAtom(self, atomId):
        """
            Get an Atom in the PDB by its id

            :param atomId: Id of the Atom (in the format "atomserial:atomName:resname")
            :type atomId: str
            :returns: int -- Number of atoms in the PDB
            :raises: KeyError if the id is not in the PDB
        """
        return self.atoms[atomId]

    def extractCOM(self):
        """
            Calculate the PDB's center of mass

            :returns: list -- List with the center of mass coordinates
        """
        if not self.totalMass:
            self.computeTotalMass()
        COM = [0., 0., 0.]
        for atomId, atom in self.atoms.items():
            COM[0] += atom.mass * atom.x
            COM[1] += atom.mass * atom.y
            COM[2] += atom.mass * atom.z

        COM[0] /= self.totalMass
        COM[1] /= self.totalMass
        COM[2] /= self.totalMass
        self.com = COM
        return COM

    def getCOM(self):
        """
            Get the PDB's center of mass

            :returns: list -- List with the center of mass coordinates
        """
        if self.com is None:
            return self.extractCOM()
        else:
            return self.com

    def extractCentroid(self):
        """
            Calculate the PDB centroid

            :returns: List -- List with the centroid coordinates
        """
        centroid = [0., 0., 0.]
        for atomId, atom in self.atoms.items():
            centroid[0] += atom.x
            centroid[1] += atom.y
            centroid[2] += atom.z

        n = float(len(self.atoms))
        centroid[0] /= n
        centroid[1] /= n
        centroid[2] /= n
        self.centroid = centroid
        return centroid

    def getCentroid(self):
        """
            Get the PDB's centroid

            :returns: list -- List with the centroid coordinates
        """
        if self.centroid is None:
            return self.extractCentroid()
        else:
            return self.centroid

    def writePDB(self, path):
        """
            Write the pdb contents of the file from wich the PDB object was
            created

            :param path: Path of the file where to write the pdb
            :type path: str
        """
        file = open(path, 'w')
        file.write(self.pdb)
        file.close()

    def countContacts(self, ligandResname, contactThresholdDistance):
        """
            Count the number of alpha carbons that are in contact with the
            protein (i.e. less than contactThresholdDistance Amstrogms away)

            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: str
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        contactThresholdDistance2 = contactThresholdDistance**2

        ligandPDB = PDB()
        ligandPDB.initialise(self.pdb, resname=ligandResname, heavyAtoms=True)

        alphaCarbonsPDB = PDB()
        alphaCarbonsPDB.initialise(self.pdb, type=self._typeProtein,
                                   atomname="CA")

        # count contacts
        contacts = set([])
        for ligandAtomId, ligandAtom in ligandPDB.atoms.items():
            # can be optimised with cell list
            for proteinAtomId, proteinAtom in alphaCarbonsPDB.atoms.items():
                dist2 = ligandAtom.squaredDistance(proteinAtom)
                if dist2 < contactThresholdDistance2:
                    contacts.update([proteinAtomId])

        return len(contacts)

    def createContactMap(self, ligandResname, contactThresholdDistance):
        """
            Create the contact map of the protein and ligand. The contact map is
            a boolean matrix that has as many rows as the number of ligand heavy
            atoms and as many columns as the number of alpha carbons. The
            value is one if the ligand atom and the alpha carbons are less than
            contactThresholdDistance Amstrongs away

            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: str
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: numpy.Array -- The contact map of the ligand and the protein
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        contactThresholdDistance2 = contactThresholdDistance**2

        ligandPDB = PDB()
        ligandPDB.initialise(self.pdb, resname=ligandResname, heavyAtoms=True)

        alphaCarbonsPDB = PDB()
        alphaCarbonsPDB.initialise(self.pdb, type=self._typeProtein,
                                   atomname="CA")

        # empty contact map, rows are atoms of the ligand, columns are protein
        # alpha carbons
        contactMap = np.zeros((len(ligandPDB.atomList),
                               len(alphaCarbonsPDB.atomList)), dtype=bool)
        contacts = set([])
        for rowind, ligandAtom in enumerate(ligandPDB.atomList):
            for colind, proteinAtom in enumerate(alphaCarbonsPDB.atomList):
                dist2 = ligandPDB.atoms[ligandAtom].squaredDistance(alphaCarbonsPDB.atoms[proteinAtom])
                if dist2 < contactThresholdDistance2:
                    contactMap[rowind, colind] = True
                    contacts.update([proteinAtom])
        return contactMap, len(contacts)

    def contactMapnew(self,ligandResname, contactThresholdDistance):
        contactThresholdDistance2 = contactThresholdDistance**2

        ligandPDB = PDB()
        ligandPDB.initialise(self.pdb, resname=ligandResname, heavyAtoms=True)

        alphaCarbonsPDB = PDB()
        alphaCarbonsPDB.initialise(self.pdb, type=self._typeProtein,
                                   atomname="CA")
        nLigand = len(ligandPDB.atomList)
        nAlpha = len(alphaCarbonsPDB.atomList)

        alphaCarbonCoords = map(lambda x: alphaCarbonsPDB.atoms[x].getAtomCoords(), alphaCarbonsPDB.atomList)
        alphaCarbonCoords = np.array(alphaCarbonCoords)
        ligandPositions = []
        alphaPositions = []
        for rowind, ligandAtomId in enumerate(ligandPDB.atomList):
            ligandAtom = ligandPDB.atoms[ligandAtomId]
            distanceCoords = (alphaCarbonCoords - ligandAtom.getAtomCoords())**2
            distanceCoords = distanceCoords.sum(axis=1)
            alphaPosCurrent = np.where(distanceCoords < contactThresholdDistance2)[0]
            alphaPositions.extend(alphaPosCurrent)
            ligandPositions.extend(np.zeros_like(alphaPosCurrent)+rowind)
        data = np.ones_like(alphaPositions, dtype=bool)
        Coo_matrix = sparse.coo_matrix((data, (ligandPositions, alphaPositions)),
                                       shape=(nLigand, nAlpha), dtype=bool)
        contactMap = Coo_matrix.todense()
        nContacts = len(np.where(contactMap.sum(axis=0))[0])
        return contactMap, Coo_matrix

def computeRMSD2(PDB1, PDB2, symmetries={}):
    """
        Compute the squared RMSD between two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :param symmetries: Dictionary with elements atomId:symmetricalAtomId corresponding with the symmetrical atoms
        :type symmetries: dict
        :returns: float -- The squared RMSD between two PDB
    """
    rmsd = 0
    for atom1Id, atom1 in PDB1.atoms.iteritems():
        d2 = []
        # HANDLE THE CASE WHEN ATOM2 IS NOT FOUND
        atom2 = PDB2.getAtom(atom1Id)
        d2 = atom1.squaredDistance(atom2)
        if symmetries:
            symAtom2Id = symmetries.get(atom1Id)
            if symAtom2Id:
                symAtom2 = PDB2.getAtom(symAtom2Id)
                d2sym = atom1.squaredDistance(symAtom2)
            else:
                d2sym = d2
        else:
            d2sym = d2

        rmsd += min(d2, d2sym)
    n = len(PDB1.atoms.items())
    return rmsd/n


def computeRMSD(PDB1, PDB2, symmetries={}):
    """
        Compute the RMSD between two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :param symmetries: Dictionary with elements atomId:symmetricalAtomId corresponding with the symmetrical atoms
        :type symmetries: dict
        :returns: float -- The squared RMSD between two PDB
    """
    return np.sqrt(computeRMSD2(PDB1, PDB2, symmetries))


def computeCOMDifference(PDB1, PDB2):
    """
        Compute the difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The distance between the centers of mass between two PDB
    """
    return computeCOMDifference(PDB1, PDB2)

def computeCOMSquaredDifference(PDB1, PDB2):
    """
        Compute the squared difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The squared distance between the centers of mass between two PDB
    """
    COM1 = PDB1.getCOM()
    COM2 = PDB2.getCOM()

    dx =  COM1[0] - COM2[0]
    dy =  COM1[1] - COM2[1]
    dz =  COM1[2] - COM2[2]

    return dx*dx + dy*dy + dz*dz

def computeSquaredCentroidDifference(PDB1, PDB2):
    """
        Compute the centroid squared difference between two PDBs

        :param PDB1: First PDB
        :type PDB1: PDB
        :param PDB2: Second PDB
        :type PDB2: PDB
        :returns: float -- The squared centroid distance between two PDB
    """
    centroid1 = PDB1.getCentroid()
    centroid2 = PDB2.getCentroid()

    dx =  centroid1[0] - centroid2[0]
    dy =  centroid1[1] - centroid2[1]
    dz =  centroid1[2] - centroid2[2]

    return dx*dx + dy*dy + dz*dz

def readPDB(pdbfile):
    """
        Helper function, parses a string with PDB content or the path of a pdb file into a string

        :param pdbfile: A string with PDB content or the path of a pdb file
        :type pdbfile: str
        :returns: str -- A string with PDB content
    """
    try:
        return open(pdbfile, "r").read()
    except IOError:
        return pdbfile
