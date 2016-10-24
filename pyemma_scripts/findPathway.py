import MSMblocks
import numpy as np

control_file = "controlFile_pathway.conf"
params = MSMblocks.readParams(control_file)
trajectoryFolder = params["trajectoryFolder"]
trajectoryBasename = params["trajectoryBasename"]
numClusters = params["numClusters"]
# get microstates
prepareMSM = MSMblocks.PrepareMSM(numClusters, trajectoryFolder, trajectoryBasename)
cl = prepareMSM.getClusteringObject()

# obtain a connectivity matrix from microstates (pyemma?)


# use graph algorithm to establish a path
