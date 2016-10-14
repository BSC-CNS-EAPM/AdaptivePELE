import os
import socket
machine = socket.gethostname()
if "bsccv" in machine:
    PELE_EXECUTABLE  = "/data/EAPM/PELE/PELE++/bin/rev12025/Pele_rev12025_mpi"
    DATA_FOLDER  = "/data/EAPM/PELE/PELE++/data/rev12025/Data"
    DOCUMENTS_FOLDER  = "/data/EAPM/PELE/PELE++/Documents/rev12025"

    PYTHON = "/data2/apps/PYTHON/2.7.5/bin/python2.7"
else:
    PELE_EXECUTABLE  = "/gpfs/projects/bsc72/PELE++/bin/rev12025/Pele_rev12025_mpi"
    DATA_FOLDER  = "/gpfs/projects/bsc72/PELE++/data/rev12025/Data"
    DOCUMENTS_FOLDER  = "/gpfs/projects/bsc72/PELE++/Documents/rev12025"
    PYTHON = "python"

inputFileTemplate = "{ \"files\" : [ { \"path\" : \"%s\" } ] }"
trajectoryBasename = "*traj*"

class OutputPathConstants():
    """
        Class with constants that depend on the outputPath
    """
    def __init__(self, outputPath):
        self.buildOutputPathConstants(outputPath)

        self.tmpFolder = "tmp_" + outputPath.replace("/", "_")

        self.buildTmpFolderConstants(self.tmpFolder)

    def buildOutputPathConstants(self, outputPath):
        self.originalControlFile = os.path.join(outputPath, "originalControlFile.conf")
        self.outputPathTempletized = os.path.join(outputPath, "%d")
        self.clusteringOutputDir = os.path.join(self.outputPathTempletized, "clustering")
        self.clusteringOutputObject = os.path.join(self.clusteringOutputDir, "object.pkl")

    def buildTmpFolderConstants(self, tmpFolder):
        self.tmpInitialStructuresTemplate = tmpFolder+"/initial_%d_%d.pdb"
        self.tmpControlFilename = tmpFolder+"/controlFile%d.conf"

