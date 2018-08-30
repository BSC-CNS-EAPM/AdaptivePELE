import os
import time
import fcntl


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
        self.createLockFile()
        self.status = self.INIT

    def createLockFile(self):
        """
            Create the lock file and write the information for the current
            process
        """
        if not os.path.exists(self.lockFile):
            # create file
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
            info[pid] = (id_num, label)
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
