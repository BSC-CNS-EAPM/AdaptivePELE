import atomset
import numpy as np


class SymmetryContactMapEvaluator:
    def __init__(self, symmetries=[]):
        """
            :param symmetries: List of dictionaries with gropus of symmetric atoms atomId:symmetricalAtomId corresponding with the symmetrical atoms
            :type symmetries: list
        """
        self.symmetricAtoms = set()
        for group in symmetries:
            self.symmetricAtoms.update(set(group.keys()).union(set(group.values())))
        self.symmetries = symmetries
        self.symToRowMap = {}

    def createContactMap(self, PDB, ligandResname, contactThresholdDistance):
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
        if self.symToRowMap:
            return self.buildContactMap(PDB, ligandResname, contactThresholdDistance)
        else:
            return self.buildContactMapWithRowMap(PDB, ligandResname, contactThresholdDistance)

    def buildContactMapWithRowMap(self, PDBobj, ligandResname, contactThresholdDistance):
        """
            Create a map that relates symmetric atoms in the ligand to row
            indices in the contactMap matrix and the contact map of the protein
            and ligand. The contact map is a boolean matrix that has as many
            rows as the number of ligand heavy atoms and as many columns as the
            number of alpha carbons. The value is one if the ligand atom and
            the alpha carbons are less than contactThresholdDistance Amstrongs
            away

            :param PDBobj: PDB object that contains the conformation's info
            :type PDBobj: PDB
            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: str
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: numpy.Array -- The contact map of the ligand and the protein
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        contactThresholdDistance2 = contactThresholdDistance**2

        ligandPDB = atomset.PDB()
        ligandPDB.initialise(PDBobj.pdb, resname=ligandResname,
                             heavyAtoms=True)

        alphaCarbonsPDB = atomset.PDB()
        alphaCarbonsPDB.initialise(PDBobj.pdb, type="PROTEIN", atomname="CA")

        # empty contact map, rows are atoms of the ligand, columns are protein
        # alpha carbons
        contactMap = np.zeros((len(ligandPDB.atomList),
                               len(alphaCarbonsPDB.atomList)), dtype=bool)
        contacts = set([])
        for rowind, ligandAtom in enumerate(ligandPDB.atomList):
            if ligandAtom in self.symmetricAtoms:
                self.symToRowMap[ligandAtom] = rowind
            for colind, proteinAtom in enumerate(alphaCarbonsPDB.atomList):
                dist2 = ligandPDB.atoms[ligandAtom].squaredDistance(alphaCarbonsPDB.atoms[proteinAtom])
                if dist2 < contactThresholdDistance2:
                    contactMap[rowind, colind] = True
                    contacts.update([proteinAtom])
        return contactMap, len(contacts)

    def buildContactMap(self, PDBobj, ligandResname, contactThresholdDistance):
        """
            Create the contact map of the protein and ligand. The contact map is
            a boolean matrix that has as many rows as the number of ligand heavy
            atoms and as many columns as the number of alpha carbons. The
            value is one if the ligand atom and the alpha carbons are less than
            contactThresholdDistance Amstrongs away

            :param PDBobj: PDB object that contains the conformation's info
            :type PDBobj: PDB
            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: str
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: numpy.Array -- The contact map of the ligand and the protein
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        contactThresholdDistance2 = contactThresholdDistance**2

        ligandPDB = atomset.PDB()
        ligandPDB.initialise(PDBobj.pdb, resname=ligandResname,
                             heavyAtoms=True)

        alphaCarbonsPDB = atomset.PDB()
        alphaCarbonsPDB.initialise(PDBobj.pdb, type="PROTEIN", atomname="CA")

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

    def evaluateJaccard(self, contactMap, cluster):
        permContactMap = self.buildOptimalPermutationContactMap(contactMap,
                                                                cluster)
        intersectContactMaps = (permContactMap == cluster.contactMap).sum()
        unionContactMaps = permContactMap.size + cluster.contactMap.size - intersectContactMaps
        similarity = float(intersectContactMaps)/unionContactMaps
        distance = 1-similarity
        return distance

    def evaluateCorrelation(self, contactMap, cluster):
        permContactMap = self.buildOptimalPermutationContactMap(contactMap,
                                                                cluster)
        similarity = calculateCorrelationContactMaps(permContactMap, cluster.contactMap)
        similarity += 1  # Necessary to omit negative correlations
        similarity /= 2.0  # Correlation values need to be higher now
        distance = 1-similarity
        return distance

    def evaluateDifferenceDistance(self, contactMap, cluster):
        permContactMap = self.buildOptimalPermutationContactMap(contactMap,
                                                                cluster)
        differenceContactMaps = np.abs(permContactMap-cluster.contactMap).sum()
        averageContacts = (0.5*(permContactMap.sum()+cluster.contactMap.sum()))
        if not averageContacts:
            # The only way the denominator can be 0 is if both contactMaps are
            # all zeros, thus being equal and belonging to the same cluster
            return True
        else:
            distance = differenceContactMaps/averageContacts
            return distance
    def buildOptimalPermutationContactMap(self, contactMap, cluster):
        """
            Build a permutated version of the contactMap which maximizes the
            similarity between the two contactMaps
            :param contactMap: contactMap of the conformation to compare
            :type contactMap: numpy.Array
            :param cluster: cluster object to which we are comparing
            :type cluster: Cluster
            :returns: numpy.Array -- The permutated contact map
        """
        if not self.symmetries and not self.symToRowMap:
            return contactMap
        permContactMap = np.array(contactMap)
        for group in self.symmetries:
            for atom1Id, atom2Id in group.iteritems():
                try:
                    atom1Row = self.symToRowMap[atom1Id]
                    atom2Row = self.symToRowMap[atom2Id]
                    line11 = cluster.contactMap[atom1Row]
                    line12 = cluster.contactMap[atom2Row]
                    line21 = contactMap[atom1Row]
                    line22 = contactMap[atom2Row]
                except KeyError as err:
                    raise KeyError("Atom %s not found in symmetries" % err.message)
                d2 = (line11 == line21).sum() + (line12 == line22).sum()
                d2sm = (line12 == line21).sum() + (line11 == line22).sum()
                if d2sm > d2:
                    # swap rows using advanced slicing from numpy
                    permContactMap[[atom1Row, atom2Row], :] = contactMap[[atom2Row, atom1Row], :]
        return permContactMap


def calculateCorrelationContactMaps(contactMap, clusterContactMap):
    """
        Calculate the correlation of two contactMaps
    """
    contactMap1 = contactMap.reshape((1, -1))
    # Reshape the matrix into a 1D array to use the numpy corrcoef function
    contactMap2 = clusterContactMap.reshape((1, -1))
    total1 = contactMap1.sum()
    total2 = contactMap2.sum()
    if not total1 or not total2:
        # if any array is all zeros the correlation will be NaN
        # if both are zero, correlation is perfect (i.e 1) else it is different
        # and will be in different clusters
        return total1 == total2
    return np.corrcoef(contactMap1, contactMap2)[0, 1]
