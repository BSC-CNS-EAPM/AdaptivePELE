from __future__ import print_function
import os.path
import numpy as np

class PDBManager:
    """
    Class that loads and processes a given pdb file
    """
    def __init__(self, PDBtoLoad, resname):
        self.PDBtoLoad = PDBtoLoad
        self.resname = resname
        self.Protein = Structure(parent=None, ID="protein")
        self.Ligand = Structure(parent=None, ID=self.resname)
        self.Other = Structure(parent=None, ID="other")
        self.POSITIONS = {"DBREF": 0, "IDCODE": 1, "ATOMNAME": 2, "RESNAME": 3, "CHAINID": 4, "RESNUMBER": 5,
                          "COORDX": 6, "COORDY": 7, "COORDZ": 8, "OCUPANCY": 9, "BFACTOR": 10}
        self.bondedCYS = []
        self.loadPDB()

    def loadPDB(self):
        """
        Method that loads the pdb into memory
        """
        with open(self.PDBtoLoad, "r") as inp:
            currentStructure = None
            currentChain = None
            currentResidue = None
            currentStructureName = ''
            currentChainName = ''
            currentResidueName = ''
            currentResidueNumber = ''
            for line in inp:
                line = line.rstrip()
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    columns = [line[:6].strip(), line[6:11].strip(), line[12:16].strip(), line[17:20].strip(),
                               line[21].strip(), line[22:27].strip(), line[30:38].strip(), line[38:46].strip(),
                               line[46:54].strip(), line[54:60].strip(), line[60:66].strip()]

                    if columns[self.POSITIONS["DBREF"]] == "ATOM" and currentStructureName != "protein":
                        currentStructureName = "protein"
                        currentStructure = self.Protein
                        currentChainName = columns[self.POSITIONS["CHAINID"]]
                        currentChain = Chain(currentStructure, currentChainName)

                    elif columns[self.POSITIONS["DBREF"]] == "HETATM":
                        if columns[self.POSITIONS["RESNAME"]] == self.resname:
                            if currentStructureName != self.resname:
                                currentStructureName = self.resname
                                currentStructure = self.Ligand
                                currentChainName = columns[self.POSITIONS["CHAINID"]]
                                currentChain = Chain(currentStructure, currentChainName)
                        elif currentStructureName != "other":
                            currentStructure = self.Other
                            currentStructureName = "other"
                            currentChainName = columns[self.POSITIONS["CHAINID"]]
                            currentChain = Chain(currentStructure, currentChainName)

                    if currentChainName != columns[self.POSITIONS["CHAINID"]]:
                        currentChainName = columns[self.POSITIONS["CHAINID"]]
                        currentChain = Chain(currentStructure, currentChainName)

                    if currentResidueNumber != columns[self.POSITIONS["RESNUMBER"]]:
                        currentResidueNumber = columns[self.POSITIONS["RESNUMBER"]]
                        currentResidueName = columns[self.POSITIONS["RESNAME"]]
                        currentResidue = Residue(currentChain, currentResidueName, currentResidueNumber)
                    Atom(currentResidue, columns[self.POSITIONS["ATOMNAME"]],[columns[self.POSITIONS["COORDX"]],
                                                                              columns[self.POSITIONS["COORDY"]],
                                                                              columns[self.POSITIONS["COORDZ"]]],
                         columns[self.POSITIONS["OCUPANCY"]],
                         columns[self.POSITIONS["BFACTOR"]])

    def writePDB(self, finalpdb, *args):
        """
        function that writes a pdb of the selected structures
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
        finalpdb = os.path.join(outputpath, outputname)
        self.writePDB(finalpdb, self.Protein, self.Ligand, self.Other)
        return finalpdb

    def checkprotonation(self):
        print("Checking the Histidine protonation State for pdb %s" % self.PDBtoLoad)
        for chain in self.Protein:
            for residue in chain:
                residue.checkHISProtonationState()

    def renumber(self, starting_number=1):
        resnumber = starting_number
        for chain in self.Protein:
            for residue in chain:
                residue.renumber(resnumber)
                resnumber += 1

    def joinChains(self):
        # Tleap doesn't support chain identifiers
        # This method unifies all the chains of the protein in one single chain
        mainChain = self.Protein[0]
        for chain in self.Protein[1:]:
            mainChain.append(chain)
            self.Protein.remove(chain)

    def getDisulphideBondsforTleapTemplate(self):
        tleapString = "bond COMPLX.%s.SG COMPLX.%s.SG\n"
        bonds_to_return = []
        for cys_pair in self.bondedCYS:
            bonds_to_return.append(tleapString % (cys_pair[0].num, cys_pair[1].num))
        return "".join(bonds_to_return)

    def loadDisulphideBonds(self):
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
        print("%d disulphide bounds found" % len(self.bondedCYS))
        for cys_pair in self.bondedCYS:
            print("Disulphide bound between CYS number %s and CYS number %s" % (cys_pair[0].num, cys_pair[1].num))
            for cys in cys_pair:
                print("Cysteine number %s renamed to CYX" % cys.num)
                cys.rename("CYX")

    def preparePDBforMD(self):
        """
        Method that prepares the pdb to be used in adaptivePELE MD simulation

        """
        # Load the information of disulphidebonds
        self.loadDisulphideBonds()
        # check the protonation states of the histidines
        self.checkprotonation()
        # Make a unique chain for the protein to avoid problems with Tleap
        # Because Tleap doesn't support chain ids
        self.joinChains()
        # Renumber the pdb to renumber properly the new chain and to remove insertion codes
        self.renumber()


class PDBase:
    """
    Abstract class to hold the different pdb structures
    """
    def __init__(self, parent, ID):
        self.parent = parent
        self.id = ID
        self.childs = []

    def setHirearchy(self):
        self.parent.childs.append(self)

    def rename(self, newname):
        self.id = newname

    def __iter__(self):
        for child in self.childs:
            yield child

    def __getitem__(self, index):
        return self.childs[index]

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
        self.num = number
        self.setHirearchy()

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


class Atom(PDBase):
    def __init__(self, parent, ID, coordinates, ocupancy, Bfactor):
        PDBase.__init__(self, parent, ID)
        self.setHirearchy()
        self.coords = coordinates
        self.ocupancy = ocupancy
        self.Bfactor = Bfactor
