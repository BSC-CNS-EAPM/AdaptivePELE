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
