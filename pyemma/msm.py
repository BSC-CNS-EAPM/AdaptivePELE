import pyemma.msm as msm
import pyemma.plots as mplt

def calculateITS(trajectories, lagtimes, its_plot_file=None):
    """ Calulate the implied time-scales at the given lagtimes"""
    its_object = msm.its(trajectories, lags=lagtimes)
    its_plot = mplt.plot_implied_timescales(its_object, outfile=its_plot_file)
    return its_object

def estimateMSM(trajectories,lagtime):
    """ Estimate a MSM from the trajectories using a provided lagtime that
    should be big enough so that the relevant processes have converged.
    Return a MaximumLikelihoodMSM object"""
    MSM_object = msm.estimate_markov_model(trajectories, lagtime)
    return MSM_object

def calculatePCCA(MSM_object, numPCCA):
    """ Coarse-cluster the MSM usin numPCCA clusters. 
    Return a PCCA object"""
    MSM_object.pccs(numPCCA)
    return MSM_object
