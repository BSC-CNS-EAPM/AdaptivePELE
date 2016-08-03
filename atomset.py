import numpy as np
import re
import StringIO

class Atom:
    ATOM_WEIGHTS = {"H":1.00794,
                    "D":2.01410178, # deuterium
                    "HE":4.00,
                    "LI":6.941,
                    "BE":9.01,
                    "B":10.811,
                    "C":12.0107,
                    "N":14.0067,
                    "O":15.9994,
                    "F":18.998403,
                    "NE":20.18,
                    "NA":22.989769,
                    "MG":24.305,
                    "AL":26.98,
                    "SI":28.09,
                    "P":30.973762,
                    "S":32.065,
                    "CL":35.453,
                    "AR":39.95,
                    "K":39.0983,
                    "CA":40.078,
                    "SC":44.96,
                    "TI":47.87,
                    "V":50.94,
                    "CR":51.9961,
                    "MN":54.938045,
                    "FE":55.845,
                    "CO":58.93,
                    "NI":58.6934,
                    "CU":63.546,
                    "ZN":65.409,
                    "GA":69.72,
                    "GE":72.64,
                    "AS":74.9216,
                    "SE":78.96,
                    "BR":79.90,
                    "KR":83.80,
                    "RB":85.47,
                    "SR":87.62,
                    "Y":88.91,
                    "ZR":91.22,
                    "NB":92.91,
                    "W":95.94, #Molybdenum?  Not sure why it's not always MO
                    "MO":95.94,
                    "TC":98.0,
                    "RU":101.07,
                    "RH":102.91,
                    "PD":106.42,
                    "AG":107.8682,
                    "CD":112.411,
                    "IN":114.82,
                    "SN":118.71,
                    "SB":121.76,
                    "TE":127.60,
                    "I":126.90447,
                    "XE":131.29,
                    "CS":132.91,
                    "BA":137.33,
                    "PR":140.91,
                    "EU":151.96,
                    "GD":157.25,
                    "TB":158.93,
                    "IR":192.22,
                    "PT":195.084,
                    "AU":196.96657,
                    "HG":200.59,
                    "PB":207.2,
                    "U":238.03}

    def __init__(self, atomContent):
        #atomContent = atomContent.split()
        if len(atomContent) > 6 and (atomContent[:4] == 'ATOM' or atomContent[:6] == 'HETATM'):
            self.atomSerial = atomContent[6:11].strip()
            self.name = atomContent[12:16].strip()
            self.resname = atomContent[17:20].strip()
            self.resChain = atomContent[21]
            self.resnum = atomContent[22:26].strip()
            self.x = float(atomContent[30:38])
            self.y = float(atomContent[38:46])
            self.z = float(atomContent[46:54])

            chargePattern = re.compile("[0-9]|\+|\-")
            self.type = re.sub(chargePattern, "", atomContent[76:80]).strip().upper()
            self.mass = self.ATOM_WEIGHTS[self.type]

            if atomContent.startswith('ATOM'):
                self.protein = True
            else:
                self.protein = False

            self.id = self.atomSerial + ":" + self.name + ":" + self.resname
            #self.id = self.atomSerial

    def isHeavyAtom(self):
        return self.type != 'H'

    def isProtein(self):
        return self.protein

    def isHeteroAtom(self):
        return not self.protein

    def printAtom(self):
        print self.atomSerial, self.name, self.resname, self.resChain, self.resnum, self.x, self.y, self.z, self.type, self.mass
        #print self.atomSerial, self.name, self.resname, self.resChain, self.resnum, self.r, self.type, self.mass

    def __eq__(self, atom2):
        return self.id == atom2.id

    def __lt__(self, atom2):
        return self.serial < atom2.serial

    def __str__(self):
        return "%s: %s %s %s [%f, %f, %f] %s %f"%(self.id, self.atomSerial, self.resChain, self.resnum, self.x, self.y, self.z, self.type, self.mass)
        #return "%s: %s %s %s [%f, %f, %f] %s %f"%(self.id, self.atomSerial, self.resChain, self.resnum, self.r[0], self.r[1], self.r[2], self.type, self.mass)

    def squaredDistance(self, atom2):
        """
            Distance between two atoms
        """
        d = (self.x - atom2.x)**2 + (self.y - atom2.y)**2 + (self.z - atom2.z)**2
        return d

class PDB:
    """
        Dictionary with atoms
        {atomId: atom, ...}
        Where atomId := atomName:resName
    """
    typeProtein = "PROTEIN"
    typeHetero = "HETERO"
    typeAll = "ALL"

    def __init__(self):
        self.atoms = {}
        self.totalMass = 0
        self.pdb = ""
        self.com = 0
        
    def __eq__(self, other):
        """ Compare two pdb strings, remark lines should be ignored and only the atoms and its information should be compared"""
        pdb1 = [element.strip() for element in self.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
        pdb2 = [element.strip() for element in other.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
        return pdb1 == pdb2

    def initialise(self, PDBstr, heavyAtoms=True, resname="", atomname="", type=typeAll):
        """
            heavyAtoms=True --> just consider heavy atoms
            PDBstr may be a path to the PDB file or a string with
            the contents of the PDB
        """
        PDBContent = StringIO.StringIO(readPDB(PDBstr)) # Using StringIO
        # creates a buffer that can handle a pdb file or a string containing
        # the PDB
        self.pdb = PDBContent.read() #in case one wants to write it

        stringWithPDBContent = self.pdb.split('\n')

        for atomLine in stringWithPDBContent:

            #HUGE optimisation (not to create an atom each time ~ 2e-5 s/atom)
            if resname != "" and not atomLine[17:20].strip() == resname: continue #optimisation
            if atomname != "" and not atomLine[12:16].strip() == atomname: continue #optimisation

              
            atom = Atom(atomLine)
            
            #Here atom will be not null, empty or not. With "try", we prune empty atoms
            try:
                if (not heavyAtoms or atom.isHeavyAtom()) and\
                    (type==self.typeAll or (type==self.typeProtein and atom.isProtein()) or (type==self.typeHetero and atom.isHeteroAtom()) ):
                    self.atoms.update({atom.id:atom})
            except:
                pass

        if self.atoms == {}:
            raise ValueError('The input pdb file/string was empty, no atoms loaded!')

    def computeTotalMass(self):
        self.totalMass = 0
        for atomId, atom in self.atoms.items():
            self.totalMass += atom.mass

    def printAtoms(self):
        for atom in self.atoms.values():
            print atom #atom.printAtom()

    def getAtom(self, atomId):
        return self.atoms[atomId]

    def extractCOM(self):
        self.computeTotalMass()
        COM = np.array([0.,0.,0.])
        for atomId, atom in self.atoms.items():
            COM += atom.mass * np.array([atom.x, atom.y, atom.z])
        COM /= self.totalMass
        self.com = COM
        return COM

    def getCOM(self):
        return self.com

    def writePDB(self, path):
        file = open(path, 'w')
        file.write(self.pdb)
        file.close()

    def countContacts(self, ligandResname, contactThresholdDistance):
        contactThresholdDistance2 = contactThresholdDistance**2

        ligandPDB = PDB()
        ligandPDB.initialise(self.pdb, resname=ligandResname, heavyAtoms=True)
        
        alphaCarbonsPDB = PDB()
        alphaCarbonsPDB.initialise(self.pdb, type = self.typeProtein, atomname = "CA")

        #count contacts
        contacts = set([])
        for ligandAtomId, ligandAtom in ligandPDB.atoms.items():
            for proteinAtomId, proteinAtom in alphaCarbonsPDB.atoms.items(): #can be optimised with cell list
                dist2 = ligandAtom.squaredDistance(proteinAtom)
                if dist2 < contactThresholdDistance2:
                    contacts.update([proteinAtomId])

        return len(contacts)

def computeRMSD2(PDB1, PDB2, symmetries={}):
    """
        Uses atom.id to match atoms from different pdbs

        Symmetries: Dictionary with elements atomId:symmetricalAtomId
    """
    rmsd = 0
    for atom1Id, atom1 in PDB1.atoms.iteritems():
        d2 = []
        ##HANDLE THE CASE WHEN ATOM2 IS NOT FOUND
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
        Uses atom.id to match atoms from different pdbs

        Symmetries: Dictionary with elements atomId:symmetricalAtomId
    """
    return np.sqrt(computeRMSD2(PDB1, PDB2, symmetries))

def computeCOMDifference(PDB1, PDB2):
    COM1 = PDB1.extractCOM()
    COM2 = PDB2.extractCOM()
    return np.linalg.norm(COM1 - COM2)

def getPDBSnapshots(file, verbose=False):
    inputFile = open(file, "r")
    inputFileContent = inputFile.read()

    snapshots = inputFileContent.split("ENDMDL")[:-1]
    if not verbose: 
        return snapshots

    remarkInfo = "REMARK 000 File created using PELE++\nREMARK source            : %s\nREMARK original model nr : %d\nREMARK First snapshot is 1, not 0 (as opposed to report)\n"
    snapshotsWithInfo  = [remarkInfo%(file, i+1) + snapshot for i,snapshot in enumerate(snapshots)]
    return snapshotsWithInfo

def readPDB(pdbfile):
        ##Finish more robust PDB initialization
        
        try:
           return open(pdbfile, "r").read()
        except IOError:
           return pdbfile
