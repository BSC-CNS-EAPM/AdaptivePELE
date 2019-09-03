from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import numpy as np
from AdaptivePELE.utilities import utilities


def extendReportWithRmsd(reportFile, rmsds):
    """
        Extend a previous report file with corrected rmsd values

        :param reportFile: Report file to be corrected
        :type reportFile: np.ndarray
        :param rmsds: Rmsd corrected values
        :type rmsds: np.ndarray

        :returns: np.ndarray -- Extended report file with corrected rmsd values
    """
    newShape = reportFile.shape
    if len(rmsds.shape) < 2:
        add_cols = 1
    else:
        add_cols = rmsds.shape[1]
    fixedReport = np.zeros((newShape[0], newShape[1]+add_cols))
    fixedReport[:, :-add_cols] = reportFile
    fixedReport[:, -add_cols:] = rmsds
    return fixedReport


def process_folder(epoch, folder, trajName, reportName, output_filename, top, trajs_to_select=None):
    if epoch is None:
        allTrajs = glob.glob(os.path.join(folder, trajName))
        full_reportName = os.path.join(folder, reportName)
    else:
        allTrajs = glob.glob(os.path.join(folder, epoch, trajName))
        full_reportName = os.path.join(folder, epoch, reportName)
        epoch = int(epoch)

    allFiles = []
    for traj in allTrajs:
        trajNum = utilities.getTrajNum(traj)
        if trajs_to_select is not None and trajNum not in trajs_to_select:
            continue
        if top is not None:
            top_file = top.getTopologyFile(epoch, trajNum)
        else:
            top_file = None
        report_file = full_reportName % trajNum
        allFiles.append((traj, report_file, top_file, epoch, output_filename % trajNum))
    return allFiles
