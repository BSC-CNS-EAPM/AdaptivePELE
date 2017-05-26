import os
import numpy as np
import shutil
import glob
import checkDetailedBalance
import ownBuildMSM
import sys
from AdaptivePELE.simulation import simulationrunner
from AdaptivePELE.pyemma_scripts import computeDeltaG


class Parameters:
    def __init__(self, ntrajs, length, lagtime, nclusters, nruns, useAllTrajInFirstRun, computeDetailedBalance, trajWildcard, folderWithTraj, lagtimes=[]):
        #If ntrajs/length = None, all trajs/lengths will be used
        self.trajWildcard = trajWildcard
        self.folderWithTraj = folderWithTraj
        self.ntrajs = ntrajs
        self.length = length
        self.lagtime = lagtime
        self.nclusters = nclusters
        self.nruns = nruns
        self.useAllTrajInFirstRun = useAllTrajInFirstRun
        self.computeDetailedBalance = computeDetailedBalance
        self.lagtimes = lagtimes

def _rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

def _rmFiles(trajWildcard):
    allfiles = glob.glob(trajWildcard)
    for f in allfiles:
        _rm(f)

def _prepareWorkingControlFile(lagtime, clusters, trajectoryFolder, trajectoryBasename, workingControlFile, lagtimes):
    """
    #Unused alternative #1, need of a templetized control file
    simulationParameters = simulationrunner.SimulationParameters()
    simulationParameters.templetizedControlFile = controlFile
    sr = simulationrunner.SimulationRunner(simulationParameters)

    controlFileDictionary = {"lagtime": lagtime, "clusters": clusters}
    sr.makeWorkingControlFile(workingControlFile, controlFileDictionary)
    """

    string = "{\"trajectoryFolder\":\"%s\", \"trajectoryBasename\":\"%s\", \"numClusters\":%d, \"lagtime\":%d, \"itsOutput\":\"its.png\", \"lagtimes\":%s}"%(trajectoryFolder, trajectoryBasename, clusters, lagtime, lagtimes)
    with open(workingControlFile, 'w') as f:
        f.write(string)

def _constructMSM(workingControlFile):
    ownBuildMSM.main(workingControlFile)

def _computeDG(trajWildcard):
    deltaGLine = computeDeltaG.main(trajWildcard)
    return deltaGLine

def _getDstName(bootstrap, i, trajFile):
    # Equiv lambda
    # getDstName = lambda bootstrap, i, trajFile: "traj_.%d.dat"%i if bootstrap else os.path.split(trajFile)[-1]
    if bootstrap:
        return "traj_.%d.dat"%i
    else:
        return os.path.split(trajFile)[-1]

def copyWorkingTrajectories(fileWildcard, length=None, ntrajs=None, bootstrap=True):
    """
        Function that copies trajectories that match "fileWildcard" into the current directory.
        It may copy a subset and a part of them (length)

        fileWildcard: wildcard to match original files
        length: trajectory length to consider, if None, the full trajectory will be considered
        ntrajs: number of trajs to consider. If None, a number equal to the total  will be considered
        bootstrap: bootstrap ntrajs from the original

        Returns: writenFiles, in order to ease a posterior cleanup

        Warning! If not using bootstrap, it uses all the trajectories
    """
    allFiles = glob.glob(fileWildcard)

    if bootstrap is False:
        trajFiles = allFiles
    else:
        if ntrajs is None:
            ntrajs = len(allFiles)
        trajFiles = np.random.choice(allFiles, ntrajs)

    writenFiles = []
    for i,trajFile in enumerate(trajFiles):
        dst = _getDstName(bootstrap, i, trajFile)
        writenFiles.append(dst)

        if length is None:
            shutil.copyfile(trajFile, dst)
        else:
            traj = np.loadtxt(trajFile)
            try:
                np.savetxt(dst, traj[:length+1,:], fmt="%d\t%.4f\t%.4f\t%.4f")
            except:
                sys.exit("There is a problem with %s"%trajFile)
    return writenFiles

def _cleanupFiles(trajWildcard):
    _rmFiles("clustering_object.pkl")
    _rmFiles("MSM_object.pkl")
    _rmFiles("discretized/traj_*")
    _rmFiles("discretized/clusterCenter*")
    _rmFiles(trajWildcard)

def _setVariablesForFirstIteration(useAllTrajInFirstRun, i, ntrajs):
    if useAllTrajInFirstRun and i == 0:
        bootstrap = False
        nWorkingTrajs = None # Not necessary, just to make it explicit that all of them are used
    else:
        bootstrap = True
        nWorkingTrajs = ntrajs

    return bootstrap, nWorkingTrajs

def _copyMSMDataFromRun(i):
    try: #it may not exist
        shutil.copyfile("its.png", "its_%d.png"%i)
    except IOError:
        pass
    shutil.copyfile("discretized/clusterCenters.dat", "clusterCenters_%d.dat"%i)
    shutil.copyfile("volumeOfClusters.dat", "volumeOfClusters_%d.dat"%i)
    shutil.copyfile("clusters.pdb", "clusters_%d.pdb" % i)
    shutil.copyfile("pmf_xyzg.dat", "pmf_xyzg_%d.dat" % i)
    if i == 0:
        try:
            shutil.copyfile("db_frobenius.eps", "db_frobenius_%d.eps"%i)
            shutil.copyfile("db_abs_diff.eps", "db_abs_diff_%d.eps"%i)
            shutil.copyfile("db_flux.eps", "db_flux_%d.eps"%i)
        except:
            pass

def _printList(l, label):
    print label
    print "====="
    for el in l:
        print el

def _getMeanAndStdFromList(l, accessFunction=lambda x:x):
    values = [float(accessFunction(element)) for element in l]
    return np.mean(values), np.std(values)


def estimateDG(parameters):
    """
        Estimates the absolute binding free energy using the parameters in the Parameters object.

        It copies the trajectory files from "folderWithTraj" into the current folder makes a certain number of iterations

        Documentation needs to be expanded, but the code style aims to help readability
    """

    workingControlFile = "control_MSM.conf"
    origFilesWildcard = os.path.join(parameters.folderWithTraj, parameters.trajWildcard)

    _prepareWorkingControlFile(parameters.lagtime, parameters.nclusters, parameters.folderWithTraj, parameters.trajWildcard, workingControlFile, parameters.lagtimes)

    deltaGs = []
    detailedBalance = []
    for i in range(parameters.nruns):
        _cleanupFiles(parameters.trajWildcard)

        bootstrap, nWorkingTrajs = _setVariablesForFirstIteration(parameters.useAllTrajInFirstRun, i, parameters.ntrajs)

        copiedFiles = copyWorkingTrajectories(origFilesWildcard, parameters.length, nWorkingTrajs, bootstrap)

        _constructMSM(workingControlFile)

        deltaG = _computeDG(parameters.trajWildcard)

        deltaGs.append(deltaG)

        if parameters.computeDetailedBalance:
            avgAsymmetricFlux = checkDetailedBalance.main(folder="discretized", countsThreshold=0, lagtime=parameters.lagtime, printFigs=False)
            detailedBalance.append(avgAsymmetricFlux)

        _copyMSMDataFromRun(i)

    _cleanupFiles(parameters.trajWildcard)

    #PLOT RESULTS
    #FIX TO WORK WITH NONES
    #print "clusters: %d, ntrajs: %d, trajLength: %d, lagtime: %d"%(parameters.nclusters, parameters.ntrajs, parameters.length, parameters.lagtime)
    _printList(deltaGs, "dG")
    meanDG, stdDG = _getMeanAndStdFromList(deltaGs, lambda element: element.split()[1])
    print "dG = %f +- %f"%(meanDG, stdDG)
    _printList(detailedBalance, "Asymmetric fluxes (see D.Lecina PhD thesis for more info)")
    meanDB, stdDB = _getMeanAndStdFromList(detailedBalance) #DB from detailed balance
    print "Asymmetric flux = %f +- %f"%(meanDB, stdDB)

    return meanDG, stdDG, meanDB, stdDB

if __name__ == "__main__":
    parameters = Parameters(ntrajs=None,
                            length=None,
                            lagtime=25,
                            nclusters=50,
                            nruns=10,
                            useAllTrajInFirstRun=True,
                            computeDetailedBalance=True,
                            trajWildcard="traj_*",
                            folderWithTraj="rawData",
                            lagtimes=[1,10,25,50,100,250,500,1000])
    estimateDG(parameters)
