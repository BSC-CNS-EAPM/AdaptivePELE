from __future__ import unicode_literals
import numpy as np
import re
import mdtraj as md
from io import StringIO, open
cimport cython
cimport numpy as np
# try:
#     # Check if the basestring type if available, this will fail in python3
#     basestring
# except NameError:
#     basestring = str


cdef class Atom:
    _chargePattern = re.compile("[0-9]|\+|\-")
    _ATOM_WEIGHTS = {u"H": 1.00794,
                    u"D": 2.01410178,  # deuterium
                    u"HE": 4.00,
                    u"LI": 6.941,
                    u"BE": 9.01,
                    u"B": 10.811,
                    u"C": 12.0107,
                    u"N": 14.0067,
                    u"O": 15.9994,
                    u"F": 18.998403,
                    u"NE": 20.18,
                    u"NA": 22.989769,
                    u"MG": 24.305,
                    u"AL": 26.98,
                    u"SI": 28.09,
                    u"P": 30.973762,
                    u"S": 32.065,
                    u"CL": 35.453,
                    u"AR": 39.95,
                    u"K": 39.0983,
                    u"CA": 40.078,
                    u"SC": 44.96,
                    u"TI": 47.87,
                    u"V": 50.94,
                    u"CR": 51.9961,
                    u"MN": 54.938045,
                    u"FE": 55.845,
                    u"CO": 58.93,
                    u"NI": 58.6934,
                    u"CU": 63.546,
                    u"ZN": 65.409,
                    u"GA": 69.72,
                    u"GE": 72.64,
                    u"AS": 74.9216,
                    u"SE": 78.96,
                    u"BR": 79.90,
                    u"KR": 83.80,
                    u"RB": 85.47,
                    u"SR": 87.62,
                    u"Y": 88.91,
                    u"ZR": 91.22,
                    u"NB": 92.91,
                    u"W": 95.94,  # Molybdenum?  Not sure why it's not always MO
                    u"MO": 95.94,
                    u"TC": 98.0,
                    u"RU": 101.07,
                    u"RH": 102.91,
                    u"PD": 106.42,
                    u"AG": 107.8682,
                    u"CD": 112.411,
                    u"IN": 114.82,
                    u"SN": 118.71,
                    u"SB": 121.76,
                    u"TE": 127.60,
                    u"I": 126.90447,
                    u"XE": 131.29,
                    u"CS": 132.91,
                    u"BA": 137.33,
                    u"PR": 140.91,
                    u"EU": 151.96,
                    u"GD": 157.25,
                    u"TB": 158.93,
                    u"IR": 192.22,
                    u"PT": 195.084,
                    u"AU": 196.96657,
                    u"HG": 200.59,
                    u"PB": 207.2,
                    u"U": 238.03}
    def __init__(self, basestring atomContent=u""):
        """ Create an atom from a pdb line

            :param atomContent: Line of the pdb from which the atom will be created
            :type atomContent: basestring
        """
        # Force string attributes to be unicode strings
        self.atomSerial = u""
        self.name = u""
        self.resname = u""
        self.resChain = u""
        self.resnum = u""
        self.type = u""
        # atomContent = atomContent.split()
        if len(atomContent) > 6 and (atomContent[:4] == u'ATOM' or atomContent[:6] == u'HETATM'):
            self.atomSerial = atomContent[6:11].strip()
            self.name = atomContent[12:16].strip()
            self.resname = atomContent[17:20].strip()
            self.resChain = atomContent[21]
            self.resnum = atomContent[22:26].strip()
            self.x = float(atomContent[30:38])
            self.y = float(atomContent[38:46])
            self.z = float(atomContent[46:54])

            self.type = re.sub(self._chargePattern, u"", atomContent[76:80]).strip().upper()
            self.mass = self._ATOM_WEIGHTS[self.type]

            if atomContent.startswith(u'ATOM'):
                self.protein = True
            else:
                self.protein = False

            self.id = self.atomSerial + u":" + self.name + u":" + self.resname
            # self.id = self.atomSerial

    def __getstate__(self):
        # Copy the object's state from
        state = {u"atomSerial": self.atomSerial, u"name": self.name, u"x": self.x,
                 u"y": self.y, u"z": self.z, u"mass": self.mass, u"type": self.type,
                 u"resname": self.resname, u"resChain": self.resChain,
                 u"resnum": self.resnum, u"protein": self.protein, u"id": self.id}
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.atomSerial = state[u'atomSerial']
        self.name = state[u'name']
        self.resname = state[u'resname']
        self.resnum = state[u'resnum']
        self.resChain = state[u'resChain']
        self.type = state[u'type']
        self.id = state[u'id']
        self.mass = state[u'mass']
        self.x = state[u'x']
        self.y = state[u'y']
        self.z = state[u'z']
        self.protein = state[u'protein']

    def set_properties(self, bint isProtein, int atomSerial, basestring atomName, basestring resName, basestring resNum, float x, float y, float z, basestring element, int resChain):
        self.atomSerial = u"%d" % atomSerial
        self.name = atomName
        self.resname = resName
        self.resChain = u"%d" % resChain
        self.resnum = resNum
        self.x = x
        self.y = y
        self.z = z

        self.type = element
        self.mass = self._ATOM_WEIGHTS[self.type]

        self.protein = isProtein

        self.id = self.atomSerial + u":" + self.name + u":" + self.resname

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

    def __richcmp__(self, Atom atom2, int op):
        if op == 2:
            #equality
            return self.id == atom2.id
        elif op == 3:
            return self.id != atom2.id
        elif op == 1:
            if self.id == atom2.id:
                return True
            else:
                return self.serial < atom2.serial
        elif op == 5:
            if self.id == atom2.id:
                return True
            else:
                return self.serial > atom2.serial
        elif op == 0:
            return self.serial < atom2.serial
        elif op == 4:
            return self.serial > atom2.serial

    def __str__(self):
        return u"%s: %s %s %s [%f, %f, %f] %s %f" % (self.id, self.atomSerial,
                                                    self.resChain, self.resnum,
                                                    self.x, self.y, self.z,
                                                    self.type, self.mass)

    def getAtomCoords(self):
        """
            Get the coordinates of the atom

            :returns: numpy.Array -- Array with the coordinate of the atom
        """
        return np.array([self.x, self.y, self.z])

    def squaredDistance(self, Atom atom2):
        """
            Calculate the squared distance between two atoms

            :param atom2: Second Atom to whom the distance will be calculated
            :type atom2: Atom
            :returns: float -- The distance between the atoms
        """
        return (self.x - atom2.x)**2 + (self.y - atom2.y)**2 + (self.z - atom2.z)**2


cdef class PDB:
    _typeProtein = u"PROTEIN"
    _typeHetero = u"HETERO"
    _typeAll = u"ALL"
    _typeCM = u"CM"

    #Atoms to be used in the contact map
    CMAtoms = {u"ALA": u"empty", u"VAL": u"empty", u"LEU": u"empty", u"ILE": u"empty",
               u"MET": u"empty", u"PRO": u"empty", u"PHE": u"CZ", u"TYR": u"OH",
               u"TRP": u"CH2", u"SER": u"empty", u"THR": u"empty", u"CYS": u"empty",
               u"ASN": u"empty", u"GLN": u"empty", u"LYS": u"NZ", u"HIS": u"CE1",
               u"HIE": u"CE1", u"HID": u"CE1", u"HIP": u"CE1", u"ARG": u"NE",
               u"ASP": u"OD1", u"GLU": u"OE1", u"GLY": u"empty"}
    ATOM_LINE_TEMPLATE = u"%s%s %s %s %s%s%s   %.3f%.3f%.3f%.2f%.2f          %s   "

    def __init__(self):
        """
            Object that will contain the information of a PDB file. Has to call
            the initialise method to load the file
        """
        self.atoms = {}
        # {atomId: atom, ...}
        # Where atomId := serial:atomName:resName
        self.totalMass = 0
        # ensure every string is unicode
        self.pdb = None
        self.com = None
        self.centroid = None

        # Necessary for contactMaps
        self.atomList = []

    def __richcmp__(self, object other, int op):
        """
            Compare two pdb strings, remark lines should be ignored and only the
            atoms and its information should be compared
        """
        cdef list pdb1, pdb2
        if op == 2:
            if self.isfromPDBFile() and other.isfromPDBFile():
                pdb1 = [element.strip() for element in self.pdb.split(u'\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                pdb2 = [element.strip() for element in other.pdb.split(u'\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                return pdb1 == pdb2
            else:
                return self.atoms == other.atoms
        elif op == 3:
            if self.isfromPDBFile() and other.isfromPDBFile():
                pdb1 = [element.strip() for element in self.pdb.split('\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                pdb2 = [element.strip() for element in other.pdb.split('\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                return pdb1 != pdb2
            else:
                return self.atoms != other.atoms
        else:
            print "No boolean operator available for PDB apart from equality"

    def __getstate__(self):
        # Copy the object's state from
        state = {u"atoms": self.atoms, u"atomList": self.atomList,
                 u"com": self.com, u"centroid": self.centroid,
                 u"totalMass": self.totalMass, u"pdb": self.pdb}
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.atoms = state[u'atoms']
        self.atomList = state[u'atomList']
        self.com = state.get(u'com')
        self.centroid = state.get(u'centroid')
        self.totalMass = state[u'totalMass']
        self.pdb = state[u'pdb']

    def isfromPDBFile(self):
        return isinstance(self.pdb, basestring)

    def _initialisePDB(self, basestring PDBstr, bint heavyAtoms=True, basestring resname=u"", basestring atomname=u"", basestring type=u"ALL", basestring chain=u"", int resnum = 0):
        """
            Load the information from a PDB file or a string with the PDB
            contents

            :param PDBstr: may be a path to the PDB file or a string with the contents of the PDB
            :type PDBstr: basestring
            :param heavyAtoms: wether to consider only heavy atoms (True if onl y heavy atoms have to be considered)
            :type heavyAtoms: bool
            :param resname: Residue name to select from the pdb (will only select the residues with that name)
            :type resname: basestring
            :param atomname: Residue name to select from the pdb (will only select the atoms with that name)
            :type atomname: basestring
            :param type: type of atoms to select: may be ALL, PROTEIN or HETERO
            :type type: basestring
            :param chain: Chain name to select from the pdb (will only select the atoms with that name)
            :type chain: basestring
            :param resnum: Residue number to select from the pdb (will only select the atoms with that name)
            :type atomname: int
            :raises: ValueError if the pdb contained no atoms
        """
        cdef object PDBContent
        cdef list stringWithPDBContent
        cdef int atomLineNum
        cdef basestring atomName, resName, atomLine, resnumStr
        cdef Atom atom
        if resnum == 0:
            resnumStr = u""
        else:
            resnumStr = u"%d" % (resnum)
        PDBContent = StringIO(readPDB(PDBstr))  # Using StringIO
        # creates a buffer that can handle a pdb file or a string containing
        # the PDB
        self.pdb = PDBContent.read()  # in case one wants to write it

        stringWithPDBContent = self.pdb.split(u'\n')
        for atomLine in stringWithPDBContent:
            if not atomLine.startswith(u"ATOM") and not atomLine.startswith(u"HETATM"):
                continue
            if type == self._typeCM:
                atomName = atomLine[12:16].strip()
                resName = atomLine[17:20].strip()
                if resName not in self.CMAtoms:
                    continue
                if atomName != u"CA" and atomName != self.CMAtoms[resName]:
                    continue
            else:
                # HUGE optimisation (not to create an atom each time ~ 2e-5 s/atom)
                if resname != u"" and not atomLine[17:20].strip() == resname:
                    continue
                if atomname != u"" and not atomLine[12:16].strip() == atomname:
                    continue
                if chain != u"" and not atomLine[21:22].strip() == chain:
                    continue
                if resnumStr != u"" and not atomLine[22:26].strip() == resnumStr:
                    continue

            atom = Atom(atomLine)
            # Here atom will be not null, empty or not.
            # With "try", we prune empty atoms
            try:
                if (not heavyAtoms or atom.isHeavyAtom()) and\
                   (type == self._typeAll or type == self._typeCM or (type == self._typeProtein and atom.isProtein()) or (type == self._typeHetero and atom.isHeteroAtom())):
                        self.atoms.update({atom.id: atom})
                        self.atomList.append(atom.id)
            except:
                pass
        if self.atoms == {}:
            raise ValueError('The input pdb file/string was empty, no atoms loaded!')

    def _initialiseXTC(self, object frame, bint heavyAtoms=True, basestring resname=u"", basestring atomname=u"", basestring type=u"ALL", basestring chain=u"", int resnum = 0):
        """
            Load the information from a loaded XTC file into a  mdtraj Trajectory

            :param PDBstr: may be a path to the PDB file or a string with the contents of the PDB
            :type PDBstr: basestring
            :param heavyAtoms: wether to consider only heavy atoms (True if onl y heavy atoms have to be considered)
            :type heavyAtoms: bool
            :param resname: Residue name to select from the pdb (will only select the residues with that name)
            :type resname: basestring
            :param atomname: Residue name to select from the pdb (will only select the atoms with that name)
            :type atomname: basestring
            :param type: type of atoms to select: may be ALL, PROTEIN or HETERO
            :type type: basestring
            :param chain: Chain name to select from the pdb (will only select the atoms with that name)
            :type chain: basestring
            :param resnum: Residue number to select from the pdb (will only select the atoms with that name)
            :type atomname: int
            :raises: ValueError if the pdb contained no atoms
        """
        cdef list stringWithPDBContent
        cdef int atomLineNum, atomIndex, atomSerial, resChain
        cdef basestring atomName, resName, atomLine, resnumStr, selection_string, element
        cdef Atom atom
        # cdef np.ndarray[int, ndim=1] selection_indexes
        cdef set selection_indexes
        cdef bint isProtein
        cdef object chain_obj, atomProv
        cdef float x, y, z
        selection_string = self.createSelectionString(heavyAtoms, resname, atomname, type, chain, resnum)
        print(selection_string)
        if resnum == 0:
            resnumStr = u""
        else:
            resnumStr = str(resnum)
        self.pdb = frame  # in case one wants to write it
        selection_indexes = set(self.pdb.topology.select(selection_string))
        for chain_obj in self.pdb.topology.chains:
            resChain = chain_obj.index
            for atomProv in chain_obj.atoms:
                if atomProv.index not in selection_indexes:
                    continue
                isProtein = atomProv.is_backbone or atomProv.is_sidechain
                atomSerial = atomProv.serial
                atomName = atomProv.name
                resName = atomProv.residue.name
                resNum = "%d" % atomProv.residue.resSeq
                x = self.pdb.xyz[0, atomProv.index, 0] * 10
                y = self.pdb.xyz[0, atomProv.index, 1] * 10
                z = self.pdb.xyz[0, atomProv.index, 2] * 10
                element = atomProv.element.symbol.upper()
                atom = Atom()
                atom.set_properties(isProtein, atomSerial, atomName, resName, resNum, x, y, z, element, resChain)
                self.atoms.update({atom.id: atom})
                self.atomList.append(atom.id)
        if self.atoms == {}:
            raise ValueError('The input pdb file/string was empty, no atoms loaded!')

    def createSelectionString(self, bint heavyAtoms=True, basestring resname=u"", basestring atomname=u"", basestring type=u"ALL", basestring chain=u"", int resnum = 0):
        cdef list selection = []
        print(heavyAtoms, resname, atomname, type, chain, resnum)
        if type == u"CM":
            for res in self.CMAtoms:
                if self.CMAtoms[res] != u"empty":
                    selection.append(u"(resname %s and name %s)" % (res, self.CMAtoms[res]))
            return u"protein and (%s)" % u" or ".join(selection)
        if atomname != u"":
            selection.append(u"name %s" % atomname)
        if heavyAtoms:
            selection.append(u"not element H")
        if resname != u"":
            selection.append(u"resname %s" % resname)
        if resnum != 0:
            selection.append(u"residue %d" % resnum)
        if type == u"HETERO":
            selection.append(u"not protein")
        elif type == u"PROTEIN":
            selection.append(u"protein")
        if selection != []:
            return u" and ".join(selection)
        else:
            return u"all"

    def initialise(self, object coordinates, bint heavyAtoms=True, basestring resname=u"", basestring atomname=u"", basestring type=u"ALL", basestring chain=u"", int resnum = 0):
        """
            Wrapper function
        """
        if isinstance(coordinates, basestring):
            self._initialisePDB(coordinates, heavyAtoms, resname, atomname, type, chain, resnum)
        else:
            self._initialiseXTC(coordinates, heavyAtoms, resname, atomname, type, chain, resnum)

    def computeTotalMass(self):
        """
            Calculate the total mass of the PDB
        """
        cdef int atomNum
        self.totalMass = 0.0
        for atomNum in range(len(self.atomList)):
            atom = self.atoms[self.atomList[atomNum]]
            self.totalMass += atom.mass

    def printAtoms(self):
        """
            Print Atom information for all the atoms in the PDB
        """
        cdef Atom atom
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
            :type atomId: basestring
            :returns: int -- Number of atoms in the PDB
            :raises: KeyError if the id is not in the PDB
        """
        return self.atoms[atomId]

    def __len__(self):
        return len(self.atomList)


    def __getitem__(self, atomId):
        return self.atoms[atomId]

    def __setitem__(self, atomId, atom):
        self.atoms[atomId] = atom

    def __delitem__(self, atomId):
        self.atoms.pop(atomId)
        self.atomList.remove(atomId)

    def __iter__(self):
        for atomId in self.atomList:
            yield self.atoms[atomId]

    def extractCOM(self):
        """
            Calculate the PDB's center of mass

            :returns: list -- List with the center of mass coordinates
        """
        if not self.totalMass:
            self.computeTotalMass()
        cdef list COM
        cdef int atomNum
        COM = [0., 0., 0.]
        for atomName in self.atoms:
            atom = self.atoms[atomName]
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
        cdef list centroid
        cdef double n
        cdef int atomNum
        centroid = [0., 0., 0.]
        for atomNum in range(len(self.atomList)):
            atom = self.atoms[self.atomList[atomNum]]
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

    def writePDB(self, basestring path, list topology=[]):
        """
            Write the pdb contents of the file from wich the PDB object was
            created

            :param path: Path of the file where to write the pdb
            :type path: basestring
        """
        cdef object fileHandle, atom
        cdef basestring prevLine = None
        if self.isfromPDBFile():
            with open(path, 'w', encoding="utf-8") as fileHandle:
                fileHandle.write(self.pdb)
        else:
            with open(path, 'w', encoding="utf-8") as fileHandle:
                # This might be problematic?, for the moment write 1
                fileHandle.write("MODEL 1\n")
                for line, atom in zip(topology, self.pdb.topology.atoms):
                    if prevLine is not None and (prevLine[21] != line[21] or (prevLine[22:26] != line[22:26] and ("HOH" == line[17:20] or "HOH" == prevLine[17:20]))):
                        fileHandle.write("TER\n")
                    x, y, z = tuple(self.pdb.xyz[0, atom.index])
                    x = (u"%.3f" % x).rjust(8)
                    y = (u"%.3f" % y).rjust(8)
                    z = (u"%.3f" % z).rjust(8)
                    fileHandle.write(line % (x, y, z))
                    prevLine = line
                fileHandle.write("ENDMDL\n")
                fileHandle.write("END\n")

    def countContacts(self, basestring ligandResname, int contactThresholdDistance, int ligandResnum=0, basestring ligandChain=""):
        """
            Count the number of alpha carbons that are in contact with the
            protein (i.e. less than contactThresholdDistance Amstrogms away)

            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: basestring
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        cdef double contactThresholdDistance2,dist2
        contactThresholdDistance2= contactThresholdDistance**2

        cdef PDB ligandPDB, alphaCarbonsPDB

        ligandPDB = PDB()
        ligandPDB.initialise(self.pdb, resname=ligandResname, resnum=ligandResnum, chain=ligandChain, heavyAtoms=True)

        alphaCarbonsPDB = PDB()
        alphaCarbonsPDB.initialise(self.pdb, type=self._typeProtein,
                                   atomname=u"CA")
        # count contacts
        cdef set contacts = set([])
        cdef int rowind, colind
        cdef basestring proteinAtomId
        cdef Atom ligandAtom, proteinAtom
        for rowind in range(len(ligandPDB.atomList)):
        # can be optimised with cell list
            ligandAtom = ligandPDB.atoms[ligandPDB.atomList[rowind]]
            for colind in range(len(alphaCarbonsPDB.atomList)):
                proteinAtomId = alphaCarbonsPDB.atomList[colind]
                proteinAtom = alphaCarbonsPDB.atoms[proteinAtomId]
                dist2 = ligandAtom.squaredDistance(proteinAtom)
                if dist2 < contactThresholdDistance2:
                    contacts.update([proteinAtomId])

        return len(contacts)


def computeCOMDifference(PDB1, PDB2):
    """
        Compute the difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The distance between the centers of mass between two PDB
    """
    return np.sqrt(computeCOMSquaredDifference(PDB1, PDB2))


def computeCOMSquaredDifference(PDB PDB1, PDB PDB2):
    """
        Compute the squared difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The squared distance between the centers of mass between two PDB
    """
    cdef list COM1, COM2
    COM1 = PDB1.getCOM()
    COM2 = PDB2.getCOM()

    dx = COM1[0] - COM2[0]
    dy = COM1[1] - COM2[1]
    dz = COM1[2] - COM2[2]

    return dx*dx + dy*dy + dz*dz


def computeSquaredCentroidDifference(PDB PDB1, PDB PDB2):
    """
        Compute the centroid squared difference between two PDBs

        :param PDB1: First PDB
        :type PDB1: PDB
        :param PDB2: Second PDB
        :type PDB2: PDB
        :returns: float -- The squared centroid distance between two PDB
    """
    cdef list centroid1, centroid2
    centroid1 = PDB1.getCentroid()
    centroid2 = PDB2.getCentroid()

    dx = centroid1[0] - centroid2[0]
    dy = centroid1[1] - centroid2[1]
    dz = centroid1[2] - centroid2[2]

    return dx*dx + dy*dy + dz*dz

def readPDB(pdbfile):
    """
        Helper function, parses a string with PDB content or the path of a pdb file into a string

        :param pdbfile: A string with PDB content or the path of a pdb file
        :type pdbfile: basestring
        :returns: basestring -- A string with PDB content
    """
    try:
        return open(pdbfile, "rt").read()
    except IOError:
        return pdbfile
