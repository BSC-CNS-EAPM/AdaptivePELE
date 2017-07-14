import os
import shutil
import glob
import numpy as np
from AdaptivePELE.testing import estimateDGAdaptive, prepareMSMFolders


def move(listFiles, dest):
    for element in listFiles:
        shutil.move(element, dest)

# list of tuples with format (lagtime(k), #clusters(tau))
iterations = [(50, 100), (50, 200), (50, 400), (100, 200), (100, 400)]
trajsPerEpoch = 50
runFolder = os.getcwd()
print "Running from " + runFolder
for tau, k in iterations:
    destFolder = "%d/%dcl" % (tau, k)
    prepareMSMFolders.main()
    os.makedirs(destFolder)
    foldersToMove = np.array(glob.glob("MSM_*"))
    epochs = [int(folder[4:]) for folder in foldersToMove]
    args = np.argsort(epochs)
    sortedFolders = foldersToMove[args]
    move(sortedFolders, destFolder)
    os.chdir(destFolder)
    print "***************"
    print "Estimating dG value in folder" + os.getcwd()
    estimateDGAdaptive.main(trajsPerEpoch, tau, k)
    os.chdir(runFolder)
