from __future__ import absolute_import, division, print_function, unicode_literals
from six import reraise as raise_
import sys
import os
import shutil
import glob
import numpy as np
from AdaptivePELE.freeEnergies import estimateDGAdaptive, prepareMSMFolders


def move(listFiles, dest):
    for element in listFiles:
        shutil.move(element, dest)

# list of tuples with format (lagtime(k), #clusters(tau))
iterations = [(25, 100), (50, 100), (100, 100), (200, 100), (400, 100),
              (25, 200), (50, 200), (100, 200), (200, 200), (400, 200),
              (25, 400), (50, 400), (100, 400), (200, 400), (400, 400)]
trajsPerEpoch = 239
nruns = 10
runFolder = os.getcwd()
print("Running from " + runFolder)
for tau, k in iterations:
    destFolder = "%dlag/%dcl" % (tau, k)
    prepareMSMFolders.main()
    if not os.path.exists(destFolder):
        os.makedirs(destFolder)
    print("***************")
    print("Estimating dG value in folder" + os.getcwd())
    try:
        estimateDGAdaptive.main(trajsPerEpoch, tau, k, nruns=nruns)
    except Exception as err:
        if "distribution contains entries smaller" in err.message:
            print("Caught exception in step with lag %d and k %d, moving to next iteration" % (tau, k))
            with open("error.txt", "w") as fe:
                fe.write("Caught exception in step with lag %d and k %d, moving to next iteration\n" % (tau, k))
        else:
            raise_(type(err), str(err), sys.exc_info()[2])

    foldersToMove = np.array(glob.glob("MSM_*"))
    epochs = [int(folder[4:]) for folder in foldersToMove]
    args = np.argsort(epochs)
    sortedFolders = foldersToMove[args]
    move(sortedFolders, destFolder)
    shutil.move("results.txt", destFolder)
