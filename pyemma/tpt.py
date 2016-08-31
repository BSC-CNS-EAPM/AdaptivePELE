import pyemma.msm as msm

def selectTPTSets(MSM_object, indexA, indexB):
    """ Extract from the sets of the PCCA clustering the sets that will serve
    as the extrems of the TPT
    """
    PCCAsets = MSM_object.metastable_sets
    SetA = PCCAsets[indexA]
    SetB = PCCAsets[indexB]
    return SetA, SetB

def createTPT(MSM_object, A, B):
    """ Calculate the reactive flux between sets A and B.
    Return a ReactiveFlux object"""
    return msm.tpt(MSM_object, A, B)
