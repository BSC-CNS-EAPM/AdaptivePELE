import numpy as np


class RMSDCalculator:
    def __init__(self, symmetries=[]):
        """
            :param symmetries: List of dictionaries with gropus of symmetric atoms atomId:symmetricalAtomId corresponding with the symmetrical atoms
            :type symmetries: list
        """
        self.nonSymmetricalAtomsSet = None
        self.symmetries = symmetries

    def computeNonSymmAtoms(self, PDB):
        allAtomsSet = set(PDB.atoms.keys())
        for group in self.symmetries:
            symmetriesSet = set(group.keys()).union(set(group.values()))
            allAtomsSet -= symmetriesSet
        self.nonSymmetricalAtomsSet = allAtomsSet

    def computeRMSD2(self, PDB1, PDB2):
        """
            Compute the squared RMSD between two PDB

            :param PDB1: First PDB with which the RMSD will be calculated
            :type PDB1: PDB
            :param PDB2: First PDB with which the RMSD will be calculated
            :type PDB2: PDB
            :returns: float -- The squared RMSD between two PDB
        """
        rmsd = 0
        if self.nonSymmetricalAtomsSet is None:
            self.computeNonSymmAtoms(PDB1)

        for group in self.symmetries:
            d2 = 0
            d2sm = 0
            for atom1Id, atom2Id in group.iteritems():
                try:
                    atom11 = PDB1.getAtom(atom1Id)
                    atom12 = PDB1.getAtom(atom2Id)
                    atom21 = PDB2.getAtom(atom1Id)
                    atom22 = PDB2.getAtom(atom2Id)
                except KeyError as err:
                    raise KeyError("Atom %d not found in PDB" % err.message)
                d2 += atom11.squaredDistance(atom21) + atom12.squaredDistance(atom22)
                d2sm += atom12.squaredDistance(atom21) + atom11.squaredDistance(atom22)
            rmsd += min(d2, d2sm)
        for atomId in self.nonSymmetricalAtomsSet:
            try:
                atom1 = PDB1.getAtom(atomId)
                atom2 = PDB2.getAtom(atomId)
            except KeyError as err:
                raise KeyError("Atom %d not found in PDB" % err.message)
            rmsd += atom1.squaredDistance(atom2)
        n = len(PDB1.atoms.items())
        return rmsd/n

    def computeRMSD(self, PDB1, PDB2):
        """
            Compute the RMSD between two PDB

            :param PDB1: First PDB with which the RMSD will be calculated
            :type PDB1: PDB
            :param PDB2: First PDB with which the RMSD will be calculated
            :type PDB2: PDB
            :returns: float -- The squared RMSD between two PDB
        """
        return np.sqrt(self.computeRMSD2(PDB1, PDB2))
