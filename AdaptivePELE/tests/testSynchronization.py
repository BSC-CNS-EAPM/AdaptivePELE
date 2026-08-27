from __future__ import absolute_import, division, print_function, unicode_literals
import shutil
import tempfile
import unittest
from AdaptivePELE.utilities import synchronization


class recordingSleep(object):
    """
        Replacement for time.sleep that records the requested intervals
        instead of waiting for them
    """
    def __init__(self):
        self.intervals = []

    def __call__(self, interval):
        self.intervals.append(interval)


class synchronizationTest(unittest.TestCase):
    def setUp(self):
        self.tmpFolder = tempfile.mkdtemp()
        # a single replica is enough, the polling loop does not depend on the
        # number of replicas and this way the constructor does not block
        self.manager = synchronization.ProcessesManager(self.tmpFolder, 1)

    def tearDown(self):
        shutil.rmtree(self.tmpFolder, ignore_errors=True)

    def runSynchronizeWithSleepRecorded(self, notSynchronizedPolls):
        """
            Run synchronize, reporting the barrier as not reached for the first
            notSynchronizedPolls checks, and return the sleep intervals used
        """
        answers = [False]*notSynchronizedPolls + [True]
        self.manager.isSynchronized = lambda status: answers.pop(0)

        sleep = recordingSleep()
        originalSleep = synchronization.time.sleep
        synchronization.time.sleep = sleep
        try:
            self.manager.synchronize("WAITING-1")
        finally:
            synchronization.time.sleep = originalSleep
        return sleep.intervals

    def testSynchronizeWaitsBetweenPolls(self):
        # a process that reaches a barrier before the others has to wait
        # between two polls of the lock file, otherwise it keeps a core busy
        # and floods the shared filesystem with lock operations
        intervals = self.runSynchronizeWithSleepRecorded(3)

        self.assertEqual(len(intervals), 3)
        self.assertTrue(all(interval > 0 for interval in intervals))

    def testSynchronizeDoesNotWaitIfBarrierIsAlreadyReached(self):
        intervals = self.runSynchronizeWithSleepRecorded(0)

        self.assertEqual(intervals, [])

    def testSynchronizeBacksOffUpToSleepTime(self):
        intervals = self.runSynchronizeWithSleepRecorded(10)

        self.assertEqual(intervals[0], self.manager.MIN_SLEEP_TIME)
        self.assertEqual(intervals, sorted(intervals))
        self.assertEqual(max(intervals), self.manager.sleepTime)
