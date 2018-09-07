import os
import time
import fcntl

try:
    ProcessLookupError
except NameError:
    ProcessLookupError = OSError


class ProcessesManager:
    """
        Object that sinchronizes multiple adaptivePELE instances, designed to
        be able to use multiple nodes of a gpu cluster
    """
    RUNNING = "RUNNING"
    WAITING = "WATING"
    INIT = "INIT"

    def __init__(self, output_path):
        self.lockFile = os.path.join(output_path, "syncFile.lock")
        self.pid = os.getpid()
        self.id = None
        self.lockInfo = {}
        self.status = self.INIT
        self.createLockFile()

    def __len__(self):
        # define the size of the ProcessesManager object as the number of
        # processes
        return len(self.lockInfo)

    def createLockFile(self):
        """
            Create the lock file and write the information for the current
            process
        """
        with open(self.lockFile, "w") as fw:
            # write something so that the file is created, it will be later
            # overwritten
            fw.write("0\n")
        file_lock = open(self.lockFile, "r+")
        while True:
            # loop until a lock can be aqcuired
            try:
                # get lock
                fcntl.lockf(file_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except (IOError, OSError):
                time.sleep(2)
        self.lockInfo = self.getLockInfo(file_lock)
        if len(self.lockInfo):
            self.id = len(self.lockInfo)
        else:
            self.id = 0
        self.lockInfo[self.pid] = (self.id, self.status)
        self.writeLockInfo(file_lock)
        fcntl.lockf(file_lock, fcntl.LOCK_UN)
        file_lock.close()

    def getLockInfo(self, file_descriptor):
        """
            Return the info stored in the lock file

            :param file_descriptor: File object of the lock file
            :type file_descriptor: file

            :returns: dict -- A dictonary containing the information of the
            different process managers initialized
        """
        info = {}
        file_descriptor.seek(0)
        for line in file_descriptor:
            if line == "0\n":
                break
            pid, id_num, label = line.rstrip().split(":")
            info[int(pid)] = (int(id_num), label)
        return info

    def writeLockInfo(self, file_descriptor):
        """
            Write the lock info to the lock file

            :param file_descriptor: File object of the lock file
            :type file_descriptor: file
        """
        file_descriptor.seek(0)
        for pid, data in self.lockInfo.items():
            file_descriptor.write("%d:%d:%s\n" % ((pid,)+data))
        file_descriptor.truncate()

    def isMaster(self):
        """
            Return wether the current process is the master process
        """
        return self.id == 0

    def setStatus(self, status):
        """
            Set the current status of the process

            :param status: Status of the process (INIT, WAITING or RUNNING)
            :type status: str
        """
        self.status = status
        self.lockInfo[self.pid] = (self.id, self.status)
        file_lock = open(self.lockFile, "r+")
        while True:
            # loop until a lock can be aqcuired
            try:
                # get lock
                fcntl.lockf(file_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except (IOError, OSError):
                time.sleep(2)
        self.writeLockInfo(file_lock)
        fcntl.lockf(file_lock, fcntl.LOCK_UN)
        file_lock.close()

    def getStatus(self):
        """
            Return the current status of the process
        """
        return self.status

    def isSynchronized(self, status):
        """
            Return wether all processes are synchronized, that is, they all
            have the same status

            :param status: Status of the process (INIT, WAITING or RUNNING)
            :type status: str

            :returns: bool -- Whether all processes are synchronized
        """
        assert len(self.lockInfo) > 0, "No processes found in lockInfo!!!"
        for pid in self.lockInfo:
            if self.lockInfo[pid][1] != status:
                return False
        return True

    def synchronize(self):
        """
            Create a barrier-like situation to wait for all processes to finish
        """
        while True:
            if self.isSynchronized(self.WAITING):
                return
            file_lock = open(self.lockFile, "r+")
            while True:
                # loop until a lock can be aqcuired
                try:
                    # get lock
                    fcntl.lockf(file_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except (IOError, OSError):
                    time.sleep(2)
            self.lockInfo = self.getLockInfo(file_lock)
            file_lock.close()

    def allRunning(self):
        """
            Check if all processes are still running
        """
        for pid in self.lockInfo:
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                return False
        return True

    def writeEquilibrationStructures(self, path, structures):
        """
            Write the equilibration structures for the current replica

            :param path: Path where to write the structures
            :type path: str
            :param structures: Filename with the structures
            :type structures: list
        """
        outfile = os.path.join(path, "structures_equilibration_%d.txt" % self.id)
        with open(outfile, "w") as fw:
            fw.write(",".join(structures))

    def readEquilibrationStructures(self, path):
        """
            Read the equilibration structures for all replicas

            :param path: Path from where to read the structures
            :type path: str
        """
        files = glob.glob(os.path.join(path, "structures_equilibration_*.txt"))
        assert len(files) == self.__len__(), "Missing files for some of the replicas"
        structures = []
        for i in range(self.__len__()):
            with open(os.path.join(path, "structures_equilibration_%d.txt" % i)) as fr:
                structure_partial = fr.read().rstrip().split(",")
            if structure_partial != [""]:
                structures.extend(structure_partial)
        return structures

    def getStructureListPerReplica(self, initialStructures, trajsPerReplica):
        """
            Filter the list of initial structures to select only the ones
            corresponding to the current replica

            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: list of str
            :param trajsPerReplica: Number of trajectories that each replica has to calculate
            :type trajsPerReplica: int

            :returns: list -- List with a tuple containing the initial structures of the replica and their indices
        """
        n_structures = len(initialStructures)
        end = min(n_structures, (self.id+1)*trajsPerReplica)
        return [(i, initialStructures[i]) for i in range(self.id*trajsPerReplica, end)]
