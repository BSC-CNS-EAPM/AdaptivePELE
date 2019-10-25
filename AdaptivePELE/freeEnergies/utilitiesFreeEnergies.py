from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from scipy import linalg


def getSortedEigen(T):
    eigenvals, eigenvectors = linalg.eig(T, left=True, right=False)
    sortedIndices = np.argsort(eigenvals)[::-1]
    return eigenvals[sortedIndices], eigenvectors[:, sortedIndices]


def getStationaryDistr(lowestEigenvector):
    absStationary = np.abs(lowestEigenvector)
    return absStationary / absStationary.sum()
