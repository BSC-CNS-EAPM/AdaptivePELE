import glob
import os
import estimateDG
import numpy as np

allFolders = np.array(glob.glob("MSM_*"))
epochs = [folder[4:] for folder in allFolders]
args = np.argsort(epochs)
sortedFolders = allFolders[args]
origDir = os.getcwd()
resultsFile=os.path.join(origDir, "results.txt")

with open(resultsFile, "a") as f:
    f.write("#Epoch DG StdDG Db StdDb\n")
    f.write("#=======================\n")

resultsEpoch = []
for epoch, folder in enumerate(sortedFolders):
    os.chdir(folder)
    parameters = estimateDG.Parameters(ntrajs=50*(epoch+1), 
                            length=None, 
                            lagtime=25, 
                            nclusters=100, 
                            nruns=10, 
                            useAllTrajInFirstRun=True, 
                            computeDetailedBalance=True, 
                            trajWildcard="traj_*", 
                            folderWithTraj="rawData")
    dG, stdDg, db, stdDb = estimateDG.estimateDG(parameters) 
    print "FINAL RESUTS EPOCH %d: dG: %f +- %f, asymmetric fluxes: %f +- %f" % (epoch, dG, stdDg, db, stdDb)
    resultsEpoch.append([dG, stdDg, db, stdDb])
    with open(resultsFile, "a") as f:
        f.write("%d %.3f %.3f %.3f %.3f\n"%(epoch, dG, stdDg, db, stdDb))
    os.chdir("..")

print "Results"
print "epoch, DG, stdDg, DB, stdDb"
print "====="
for i,results in enumerate(resultsEpoch):
    print i, results
