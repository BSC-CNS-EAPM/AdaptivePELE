from __future__ import print_function
import os.path
import numpy as np
import glob
import AdaptivePELE.constants


class PDBLoadException(Exception):
    __module__ = Exception.__module__


class PDBManager:
    """
    Class that loads and processes a given pdb file
    """
    # Residue Names that are recognized by the AMBER forcefields.

    VALID_RESNAMES = {"ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY",
                      "HID", "HIE", "HIP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR",
                      "TRP", "TYR", "VAL", "CALA", "CARG", "CASN", "CASP", "CCYS", "CCYX", "CGLN", "CGLU",
                      "CGLY", "CHID", "CHIE", "CHIP", "CILE", "CLEU", "CLYS", "CMET", "CPHE", "CPRO", "CSER",
                      "CTHR", "CTRP", "CTYR", "CVAL", "NHE", "NME", "ACE", "NALA", "NARG", "NASN", "NASP",
                      "NCYS", "NCYX", "NGLN", "NGLU", "NGLY", "NHID", "NHIE", "NHIP", "NILE", "NLEU", "NLYS",
                      "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL", "HIS"}

    VALID_NUCLEIC = {"DA", "DA3", "DA5", "DAN", "DC", "DC3", "DC5", "DCN", "DG", "DG3", "DG5", "DGN", "DT",
                     "DT3", "DT5", "DTN", "RA", "RA3", "RA5", "RAN", "RC", "RC3", "RC5", "RCN", "RG", "RG3",
                     "RG5", "RGN", "RU", "RU3", "RU5", " RUN"}

    VALID_ION = {"AG", "AL", "Ag", "BA", "BR", "Be", "CA", "CD", "CE", "CL", "CO", "CR", "CS", "CU", "CU1",
                 "Ce", "Cl-", "Cr", "Dy", "EU", "EU3", "Er", "F", "FE", "FE2", "GD3", "H3O+", "HE+", "HG",
                 "HZ+", "Hf", "IN", "IOD", "K", "K+", "LA", "LI", "LU", "MG", "MN", "NA", "NH4", "NI", "Na+",
                 "Nd", "PB", "PD", "PR", "PT", "Pu", "RB", "Ra", "SM", "SR", "Sm", "Sn", "TB", "TL", "Th",
                 "Tl", "Tm", "U4+", "V2+", "Y", "YB2", "ZN", "Zr"}

    # Modified residues for which there is a template in AdaptivePELE/constants/MDtemplates/
    TEMPLATE_PATH = os.path.join("".join(AdaptivePELE.constants.__path__), "MDtemplates/amber_*.lib")

    VALID_MODIFIED_RES = [name.split("_")[-1][:-4] for name in glob.glob(TEMPLATE_PATH)]

    def __init__(self, PDBtoLoad, resname):
        """

        :param PDBtoLoad: Path to the pdb to load into memory
        :type PDBtoLoad: str
        :param resname: Name of the ligand
        :type resname: str
        """
        self.PDBtoLoad = PDBtoLoad
        self.resname = resname
        self.Protein = Structure(parent=None, ID="protein")
        self.Ligand = Structure(parent=None, ID=self.resname)
        self.Other = Structure(parent=None, ID="other")
        # Dictionary with the information that appears in each of the pdb fields. Used to avoid using magic numbers
        self.POSITIONS = {"DBREF": 0, "IDCODE": 1, "ATOMNAME": 2, "RESNAME": 3, "CHAINID": 4, "RESNUMBER": 5,
                          "COORDX": 6, "COORDY": 7, "COORDZ": 8, "OCUPANCY": 9, "BFACTOR": 10}
        # List with the cysteines bonded with disulphite bonds
        self.bondedCYS = []
        # Set with the modified residues
        self.modified_res = set()
        # Dictionary with the valid Atom names for each of the heavy atoms of each residue
        self.AtomTemplates = self.loadTemplates()
        self.loadPDB()

    def loadTemplates(self):
        """
        Method that loads the Templates with the heavy atoms information
        :return: Dictionary with each residue as key and its heavy atom names has values
        """
        templates_dict = {}
        template = os.path.join("".join(AdaptivePELE.constants.__path__), "MDtemplates/Template_PDB")
        with open(template, "r") as temp:
            for line in temp:
                line = line.split()
                templates_dict.setdefault(line[0], set()).add(line[1])
        return templates_dict

    def loadPDB(self):
        """
        Method that loads the pdb into memory
        """
        with open(self.PDBtoLoad, "r") as inp:
            currentStructure = None
            currentChain = None
            currentResidue = None
            currentChainName = ''
            currentResidueNumber = ''
            for line in inp:
                line = line.rstrip()
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # This line splits the pdb into each one of its fields
                    columns = [line[:6].strip(), line[6:11].strip(), line[12:16].strip(), line[17:20].strip(),
                               line[21].strip(), line[22:27].strip(), line[30:38].strip(), line[38:46].strip(),
                               line[46:54].strip(), line[54:60].strip(), line[60:66].strip()]

                    # Extract the additional letter that is placed sometimes into the residues that have alternative positions
                    if columns[self.POSITIONS["OCUPANCY"]]:
                        if float(columns[self.POSITIONS["OCUPANCY"]]) < 1:
                            if len(columns[self.POSITIONS["RESNAME"]]) > 3:
                                columns[self.POSITIONS["RESNAME"]] = columns[self.POSITIONS["RESNAME"]][1:]


                    if columns[self.POSITIONS["RESNAME"]] in self.VALID_RESNAMES:
                        if currentStructure != self.Protein:
                            currentStructure = self.Protein
                            currentChainName = columns[self.POSITIONS["CHAINID"]]
                            currentChain = Chain(currentStructure, currentChainName)

                    elif columns[self.POSITIONS["RESNAME"]] in self.VALID_MODIFIED_RES:
                        self.modified_res.add(columns[self.POSITIONS["RESNAME"]])
                        if currentStructure != self.Protein:
                            currentStructure = self.Protein
                            currentChainName = columns[self.POSITIONS["CHAINID"]]
                            currentChain = Chain(currentStructure, currentChainName)

                    elif columns[self.POSITIONS["RESNAME"]] == self.resname:
                        if currentStructure != self.Ligand:
                                currentStructure = self.Ligand
                                currentChainName = columns[self.POSITIONS["CHAINID"]]
                                currentChain = Chain(currentStructure, currentChainName)

                    elif columns[self.POSITIONS["RESNAME"]] in self.VALID_ION or columns[self.POSITIONS["RESNAME"]]\
                            in ["WAT", "HOH"] or columns[self.POSITIONS["RESNAME"]] in self.VALID_NUCLEIC:
                        if currentStructure != self.Other:
                            currentStructure = self.Other
                            currentChainName = columns[self.POSITIONS["CHAINID"]]
                            currentChain = Chain(currentStructure, currentChainName)

                    else:
                        raise PDBLoadException("Residue %s of Chain %s is not in templates" % (columns[self.POSITIONS["RESNAME"]],
                                                                              columns[self.POSITIONS["CHAINID"]]))

                    if currentChainName != columns[self.POSITIONS["CHAINID"]]:
                        currentChainName = columns[self.POSITIONS["CHAINID"]]
                        currentChain = Chain(currentStructure, currentChainName)

                    if currentResidueNumber != columns[self.POSITIONS["RESNUMBER"]]:
                        currentResidueNumber = columns[self.POSITIONS["RESNUMBER"]]
                        currentResidue = Residue(currentChain, columns[self.POSITIONS["RESNAME"]], currentResidueNumber)

                    Atom(currentResidue, columns[self.POSITIONS["ATOMNAME"]],[columns[self.POSITIONS["COORDX"]],
                                                                              columns[self.POSITIONS["COORDY"]],
                                                                              columns[self.POSITIONS["COORDZ"]]],
                         columns[self.POSITIONS["OCUPANCY"]],
                         columns[self.POSITIONS["BFACTOR"]])

    def writePDB(self, finalpdb, *args):
        """
        Method that writes a pdb of the selected structures
        :param outputpath: name of the path to save the pdb
        :type: str
        :param outputname: name of the output pdb
        :type: str
        :param args: structure objects to save into the pdb
        """
        atomNumber = 0
        with open(finalpdb, "w") as out_pdb:
            for structure in args:
                if structure.id == "protein":
                    Atomtype = "ATOM"
                else:
                    Atomtype = "HETATM"
                for chain in structure:
                    for residue in chain:
                        for atom in residue:
                            atomNumber += 1
                            out_pdb.write("%-6s%5s %4s %3s %s%4s    %8s%8s%8s\n" % (Atomtype, atomNumber, atom.id, residue.id, chain.id, residue.num, atom.coords[0], atom.coords[1], atom.coords[2]))
                    out_pdb.write("TER\n")

    def writeAll(self, outputpath, outputname):
        """
        Method that writes the whole pdb

        :param outputpath: name of the path to save the pdb
        :type: str
        :param outputname: name of the output pdb
        :type: str
        :return: the path of the new pdb
        """
        finalpdb = os.path.join(outputpath, outputname)
        self.writePDB(finalpdb, self.Protein, self.Ligand, self.Other)
        return finalpdb

    def checkprotonation(self):
        """
        Method that change the HIS names to the correct name according to its actual protonation state
        """
        print("Checking the Histidine protonation State")
        for chain in self.Protein:
            for residue in chain:
                residue.checkHISProtonationState()

    def renumber(self, starting_number=1):
        """
        Renumbers the each one of the chains starting from the starting_number

        :param starting_number: initial number for the structure
        :type starting_number: int
        """
        resnumber = starting_number
        for chain in self.Protein:
            for residue in chain:
                residue.renumber(resnumber)
                resnumber += 1

    def joinChains(self):
        """
        This method unifies all the chains of the protein in one single chain.
        This is done because Tleap doesn't support chain identifiers, and, if different chains are provided,
        Tleap renumbers them in an aribtary way, making impossible to keep track of the residue numbers.
        """
        mainChain = self.Protein[0]
        for chain in self.Protein[1:]:
            mainChain.append(chain)
            self.Protein.remove(chain)

    def getDisulphideBondsforTleapTemplate(self):
        """
        This method converts the bonded cysteines information into a string to be substitued into the Tleap template
        :return: string with the bonded cysteines written with Tleap syntaxis
        """
        tleapString = "bond COMPLX.%s.SG COMPLX.%s.SG\n"
        bonds_to_return = []
        for cys_pair in self.bondedCYS:
            bonds_to_return.append(tleapString % (cys_pair[0].num, cys_pair[1].num))
        return "".join(bonds_to_return)

    def loadDisulphideBonds(self):
        """
        Method that checks the cysteines that are bonded using the euclidian distance between them as discriminator
        """
        cysteines = {}
        cheked_cys = set()
        for chain in self.Protein:
            for residue in chain:
                if residue.id == "CYS" or residue.id == "CYX":
                    for atom in residue:
                        if atom.id == "SG":
                            cysteines[residue] = atom.coords
        for cystiene in cysteines:
            cheked_cys.add(cystiene)
            for cystiene_2 in cysteines:
                if cystiene_2 not in cheked_cys:
                    coords1 = np.array([float(x) for x in cysteines[cystiene]])
                    coords2 = np.array([float(x) for x in cysteines[cystiene_2]])
                    distance = np.linalg.norm(coords1-coords2)
                    # 2.1 is the threshold for bonding
                    if distance <= 2.1:
                        self.bondedCYS.append((cystiene, cystiene_2))
        self.renameBondedCysteines()

    def renameBondedCysteines(self):
        # Method that renames the cysteines
        print("%d disulphide bounds found" % len(self.bondedCYS))
        for cys_pair in self.bondedCYS:
            print("Disulphide bound between CYS number %s and CYS number %s" % (cys_pair[0].num, cys_pair[1].num))
            for cys in cys_pair:
                print("Cysteine number %s renamed to CYX" % cys.num)
                cys.rename("CYX")

    def checkMissingAtoms(self):
        # Method that check that all the heavy atoms are in the templates and also checks that all residues have all the heavy atoms
        for chain in self.Protein:
            for i, residue in enumerate(chain):
                if residue.id in self.AtomTemplates:
                    atomsNames = set(residue.getChildNames())
                    if i == 0:
                        templateAtoms = self.AtomTemplates["N%s" % residue.id]
                    else:
                        templateAtoms = self.AtomTemplates[residue.id]
                    extra_atoms = atomsNames.difference(templateAtoms)
                    missing_atoms = templateAtoms.difference(atomsNames)
                    for atom in extra_atoms:
                        if atom != "OXT":
                            if "H" in atom:
                                print("Warning: Atom %s of Residue %s in chain %s not in Templates.\nRemoving Hydrogen" % (atom, residue.id, chain.id))
                                residue.remove(residue[atom])
                            else:
                                raise PDBLoadException("ERROR: Atom %s of Residue %s in chain %s not in Templates" % (atom, residue.id, chain.id))
                    for atom in missing_atoms:
                        print("Warning: Residue %s of chain %s doesn't have the Atom %s" %(residue.id, chain.id, atom))

    def getModifiedResiduesTleapTemplate(self):
        """
        Method that creates an string to be used in Tleap to load the forcefield parameters for the modified residues found
        :return: string with the templates to load in Tleap syntaxis
        """
        templatespath = os.path.join("".join(AdaptivePELE.constants.__path__), "MDtemplates/amber_%s.lib")
        tleapString = "loadoff %s\n" % templatespath
        stringToReturn = []
        for res in self.modified_res:
            stringToReturn.append(tleapString % res)
        return "".join(stringToReturn)

    def correctAlternativePositions(self):
        # This method selects the positions with higher occupancy when alternative positions are found
        All = (self.Protein, self.Ligand, self.Other)
        for structure in All:
            for chain in structure:
                for residue in chain:
                    residue.alternativepositions()

    def checkgaps(self):
        # Method that checks for possible gaps in the structure by checking the number of the residues
        for chain in self.Protein:
            prev_residue = 0
            for residue in chain:
                if prev_residue and residue.num > prev_residue + 1:
                    print("Warning: Possible gap found in chain %s between residue %s and %s" % (chain.id, prev_residue, residue.num))
                prev_residue = residue.num

    def checkLigand(self):
        for chain in self.Ligand:
            for residue in chain:
                for atom in residue:
                    if atom.id.startswith("CL"):
                        oldname = atom.id
                        atom.id = "Cl%s" % oldname[2:]
                        print("Atom %s of %s rename to %s" % (oldname, self.resname, atom.id))

    def preparePDBforMD(self):
        """
        Method that prepares the pdb to be used in adaptivePELE MD simulation

        """
        # Check if the structure has posible gaps
        self.checkgaps()
        # Select the positions with higher occupancy if alternative positions are found
        self.correctAlternativePositions()
        # Load the information of disulphidebonds
        self.loadDisulphideBonds()
        # check the protonation states of the histidines
        self.checkprotonation()
        # Rename atoms from the ligand to match parmchk atom names
        self.checkLigand()
        # Make a unique chain for the protein to avoid problems with Tleap
        # Because Tleap doesn't support chain ids
        self.joinChains()
        # Renumber the pdb and remove insertion codes
        self.renumber()
        # Check for missing atoms
        self.checkMissingAtoms()


class PDBase:
    """
    Abstract class to hold the different pdb structures
    """
    def __init__(self, parent, ID):
        """
        :param parent: Object that is one step above in the hirearchy.
        The used hirearchy is the following: Structure - Chain - Residue - Atom
        :param ID: Name of the object to instance
        """
        self.parent = parent
        self.id = ID
        # List with the objects that are one step bellow
        self.childs = []
        self.child_dict = {}

    def setHirearchy(self):
        # Method that loads the current object to its father childs list
        self.parent.childs.append(self)
        self.parent.child_dict[self.id] = self

    def rename(self, newname):
        # Method to rename the object
        self.id = newname

    def __iter__(self):
        for child in self.childs:
            yield child

    def __getitem__(self, index):
        try:
            return self.childs[index]
        except TypeError:
            return self.child_dict[index]

    def getChildNames(self):
        # Method that returns a list with the names of all the child objects that are one step bellow
        return [child.id for child in self.childs]

    def append(self, target):
        self.childs = self.childs + target.childs
        for child in target:
            child.parent = self

    def remove(self, target):
        self.childs.remove(target)


class Structure(PDBase):
    pass


class Chain(PDBase):
    def __init__(self, parent, ID):
        PDBase.__init__(self, parent, ID)
        self.setHirearchy()


class Residue(PDBase):
    def __init__(self, parent, ID, number):
        PDBase.__init__(self, parent, ID)
        self.num, self.insertionCode = self.set_number(number)
        self.setHirearchy()

    def set_number(self, number):
        insertion = None
        try:
            num = int(number)
        except ValueError:
            num = int(number[:-1])
            insertion = number[-1]
        return num, insertion


    def checkHISProtonationState(self):
        """
        Method that renames the histidines according to their protonation states
        """
        if self.id == "HIS":
            HD1, HE1 = False, False
            for atom in self.childs:
                if atom.id == "HD1":
                    HD1 = True
                elif atom.id == "HE1":
                    HE1 = True
            if HD1 and HE1:
                self.rename("HIP")
                print("Histidine number %s renamed to HIP" % self.num)
            elif HD1:
                self.rename("HID")
                print("Histidine number %s renamed to HID" % self.num)
            elif HE1:
                self.rename("HIE")
                print("Histidine number %s renamed to HIE" % self.num)

    def renumber(self, new_number):
        self.num = new_number

    def alternativepositions(self):
        atomsToRemove = []
        alternativepositions = {}
        for atom in self:
            if atom.ocupancy:
                if atom.ocupancy < 1:
                    alternativepositions.setdefault(atom.id, []).append(atom)
        for key in alternativepositions:
            bigger = 0
            lastAtom = None
            for value in alternativepositions[key]:
                if lastAtom:
                    if value.ocupancy > bigger:
                        bigger = value.ocupancy
                        atomsToRemove.append(lastAtom)
                        lastAtom = value
                    else:
                        atomsToRemove.append(value)
                else:
                    bigger = value.ocupancy
                    lastAtom = value
        if len(atomsToRemove) > 0:
            print("Alternative positions found on residue %s %s. Keeping the positions with higher occupancy" % (self.id, self.num))
        for atom in atomsToRemove:
            self.childs.remove(atom)


class Atom(PDBase):
    def __init__(self, parent, ID, coordinates, ocupancy, Bfactor):
        PDBase.__init__(self, parent, ID)
        self.setHirearchy()
        self.coords = coordinates
        self.ocupancy = ocupancy
        if self.ocupancy:
            self.ocupancy = float(self.ocupancy)
        self.Bfactor = Bfactor

