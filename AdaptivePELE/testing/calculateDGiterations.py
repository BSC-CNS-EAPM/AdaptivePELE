import os
import shutil
from AdaptivePELE.testing import estimateDG, prepareMSMFolders

# list of tuples with format (lagtime, #clusters)
iterations = [(50, 200), (50, 400), (100, 200), (100, 400), (200, 200), (200, 400)]
runFolder = os.getcwd()
print "Running from " + runFolder
for tau, k in iterations:
    destFolder = "%d/%dcl" % (tau, k)
    prepareMSMFolders.main()
    os.makedirs(destFolder)
    shutil.move("MSM_0", destFolder)
    os.chdir(destFolder+"/MSM_0")
    print "***************"
    print "Estimating dG value in folder" + os.getcwd()
    parameters = estimateDG.Parameters(ntrajs=None,
                            length=None,
                            lagtime=tau,
                            nclusters=k,
                            nruns=10,
                            skipFirstSteps = 0,
                            useAllTrajInFirstRun=True,
                            computeDetailedBalance=True,
                            trajWildcard="traj_*",
                            folderWithTraj="rawData",
                            lagtimes=[1,10,25,50,100,250,500,1000],
                            clusterCountsThreshold=0)
    estimateDG.estimateDG(parameters, cleanupClusterCentersAtStart=True)
    os.chdir(runFolder)
