import os
import numpy as np
import firstSnapshots
from simulation import simulationrunner
import shutil
import glob
import buildMSM
import computeDeltaG

"""
    1) Change hardcoded params, e.g. ntrajs, length, lagtime, nclusters,...
    2) Needs of a templetized_control_MSM.conf with $clusters and $lagtime as template (if wanted), otherwise, create a file with the same name and no templetized attr
"""


"""
    Params TO BE CHANGED
"""
ntrajs = 511*2
length = 2000
lagtime = 200
nclusters = 0
nruns = 10
useAllTrajInFirstRun = True

controlFile = "templetized_control_MSM.conf" #unused, but remember it needs to be defined
trajWildcard = "traj*"
folderWithTraj = "rawData"


"""
    Functions
"""
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
    buildMSM.main("control_MSM.conf")
    deltaGLine = computeDeltaG.main("traj_*")
    return deltaGLine 

"""
    Main
"""
deltaGs = {}
for i in range(nruns):
    rmTrajFiles(trajWildcard)
    if useAllTrajInFirstRun and i == 0: #this uses full length trajectories
        rm("clustering_object.pkl") 
        copyAllTrajectories(trajWildcard, folderWithTraj)
    else:
        firstSnapshots.main(length, ntrajs)
    rm("clustering_object.pkl") 
    rm("MSM_object.pkl") 
    deltaG = computeDG(lagtime, nclusters)
    try:
        deltaGs[lagtime].append(deltaG)
    except KeyError:
        deltaGs[lagtime] = [deltaG]
    shutil.copyfile("its.png", "its_%d.png"%i)
    shutil.copyfile("volumeOfClusters.dat", "volumeOfClusters_%d.dat"%i)
    shutil.copyfile("clusters.pdb", "clusters_%d.pdb"%i)
    shutil.copyfile("pmf_xyzg.dat", "pmf_xyzg_%d.dat"%i)

#SUMMARY RESULTS
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
