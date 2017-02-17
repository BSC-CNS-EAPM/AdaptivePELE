import os
import numpy as np
import firstSnapshots
from simulation import simulationrunner
import shutil
import glob
from other import checkDetailedBalance
from pyemma_scripts import computeDeltaG
import ownBuildMSM

def rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

def rmTrajFiles(trajWildcard):
    allfiles = glob.glob(trajWildcard)
    for f in allfiles:
        rm(f)

def copyAllTrajectories(trajWildcard, folderWithTraj):
    allfiles = glob.glob(os.path.join(folderWithTraj, trajWildcard))
    for f in allfiles:
        shutil.copy(f, ".")

def computeDG(lagtime, clusters):
    simulationParameters = simulationrunner.SimulationParameters()
    controlFile = "templetized_control_MSM.conf"
    simulationParameters.templetizedControlFile = controlFile
    sr = simulationrunner.SimulationRunner(simulationParameters)

    controlFileDictionary = {"lagtime": lagtime, "clusters": clusters}
    sr.makeWorkingControlFile("control_MSM.conf", controlFileDictionary)
    ownBuildMSM.main("control_MSM.conf")
    deltaGLine = computeDeltaG.main("traj_*")
    return deltaGLine 

ntrajs = 512
length = 1001
lagtime = 200
nclusters = 100
nruns = 10
useAllTrajInFirstRun = True

computeDetailedBalance = True

controlFile = "templetized_control_MSM.conf" #unused, but remember it needs to be defined
trajWildcard = "traj*"
folderWithTraj = "rawData"

deltaGs = {}
detailedBalance = {}
for i in range(nruns):
    rmTrajFiles(trajWildcard)
    if useAllTrajInFirstRun and i == 0: #this uses full length trajectories
        rm("clustering_object.pkl") 
        rmTrajFiles("discretized/clusterCenter*")
        copyAllTrajectories(trajWildcard, folderWithTraj)
    else:
        firstSnapshots.main(length, ntrajs)
    rm("clustering_object.pkl") 
    rm("MSM_object.pkl") 
    rmTrajFiles("discretized/traj_*")
    rmTrajFiles("discretized/clusterCenter*")
    deltaG = computeDG(lagtime, nclusters)
    try:
        deltaGs[lagtime].append(deltaG)
    except KeyError:
        deltaGs[lagtime] = [deltaG]

    if computeDetailedBalance:
        frobeniusAvg, relativeEntropyT, unused = checkDetailedBalance.main("discretized", 0, lagtime)
        try:
            detailedBalance['frobenius'].append(frobeniusAvg)
            detailedBalance['relativeEntropy'].append(relativeEntropyT)
        except KeyError:
            detailedBalance['frobenius'] = [frobeniusAvg]
            detailedBalance['relativeEntropy'] = [relativeEntropyT]

    shutil.copyfile("its.png", "its_%d.png"%i)
    shutil.copyfile("volumeOfClusters.dat", "volumeOfClusters_%d.dat"%i)
    shutil.copyfile("clusters.pdb", "clusters_%d.pdb"%i)
    shutil.copyfile("pmf_xyzg.dat", "pmf_xyzg_%d.dat"%i)
    if i == 0:
        try:
            shutil.copyfile("db_frobenius.eps", "db_frobenius_%d.eps"%i)
            shutil.copyfile("db_abs_diff.eps", "db_abs_diff_%d.eps"%i)
            shutil.copyfile("db_flux.eps", "db_flux_%d.eps"%i)
        except:
            pass

#PLOT RESULTS
print "clusters: %d, ntrajs: %d, trajLength: %d"%(nclusters, ntrajs, length)
for key, val in deltaGs.iteritems():
    print key
    print "====="
    dGs = []
    for element in val:
        print element
        dG = element.split()[1]
        dGs.append(float(dG))
    print "dG = %f +- %f"%(np.mean(dGs), np.std(dGs))

if computeDetailedBalance:
    for key, val in detailedBalance.iteritems():
        print key
        print "====="
        values = []
        for element in val:
            print element
            values.append(float(element))
        print "%s = %f +- %f"%(key, np.mean(values), np.std(values))
