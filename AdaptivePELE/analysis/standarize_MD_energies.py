from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import argparse
import numpy as np
import multiprocessing as mp
from sklearn import preprocessing
from AdaptivePELE.utilities import utilities
from AdaptivePELE.analysis import analysis_utils


def parseArguments():
    """
        Parse the command-line options

        :returns: object -- Object containing the options passed
    """
    desc = "Program that standarizes potential energy values from an MD simulation."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("energy_col", type=int, help="Ligand resname")
    parser.add_argument("--path", type=str, default=".", help="Path where the simulation is stored")
    parser.add_argument("--out_name", type=str, default="fixedReport", help="Name of the modified report files (default is fixedReport)")
    parser.add_argument("--out_folder", type=str, default=None, help="Path where to store the report files (default is fixedReport)")
    parser.add_argument("-n", type=int, default=1, help="Number of processors to parallelize")
    parser.add_argument("--fmt_str", type=str, default="%.4f", help="Format of the output file (default is .4f which means all floats with 4 decimal points)")
    parser.add_argument("--new_report", action="store_true", help="Whether to create new report files instead of modifying existing ones")
    parser.add_argument("--report_name", type=str, default=None, help="Name of report files, for example if reports are named 'report_1' use report")
    parser.add_argument("--traj_to_process", nargs="*", type=int, default=None, help="Number of the trajectories to filter, if not specified all of them will be processed")
    args = parser.parse_args()

    return args.energy_col, args.path, args.out_name, args.fmt_str, args.n, args.out_folder, args.new_report, args.report_name, args.traj_to_process


def process_file(report, outputFilename, format_out, new_report, epoch, energy_column):
    try:
        reportFilename = glob.glob(report)[0]
    except IndexError:
        raise IndexError("File %s not found" % report)

    with open(reportFilename) as f:
        header = f.readline().rstrip()
        if not header.startswith("#"):
            header = ""
        reportFile = utilities.loadtxtfile(f)
    energy_values = reportFile[:, energy_column]
    energy_values = preprocessing.scale(energy_values)

    if not new_report:
        if outputFilename != reportFilename:
            with open(outputFilename) as f:
                header = f.readline().rstrip()
                if not header.startswith("#"):
                    header = ""
            reportFile = utilities.loadtxtfile(outputFilename)
        fixedReport = analysis_utils.extendReportWithRmsd(reportFile, energy_values)
    else:
        header = ""
        indexes = np.array(range(energy_values.shape[0]))
        fixedReport = np.concatenate((indexes[:, None], energy_values[:, None]), axis=1)

    with open(outputFilename, "w") as fw:
        if header:
            fw.write("%s\tEnergy\n" % header)
        else:
            fw.write("# Step\tEnergy\n")
        np.savetxt(fw, fixedReport, fmt=format_out, delimiter="\t")


def main(col_energy, folder, out_report_name, format_out, nProcessors, output_folder, new_report, reportName, trajs_to_select):
    """
        Calculate the relative SASA values of the ligand

        :param col_energy: Column corresponding to the energy in the reports
        :type col_energy: int
        :param folder: Path the simulation
        :type folder: str
        :param out_report_name: Name of the output file
        :type out_report_name: str
        :param format_out: String with the format of the output
        :type format_out: str
        :param nProcessors: Number of processors to use
        :type nProcessors: int
        :param output_folder: Path where to store the new reports
        :type output_folder: str
        :param new_report: Whether to create new reports
        :type new_report: bool
    """
    # Constants
    if output_folder is not None:
        out_report_name = os.path.join(output_folder, out_report_name)
    outputFilename = "_".join([out_report_name, "%d"])
    trajName = "*traj*"
    if reportName is None:
        reportName = "report_%d"
    else:
        reportName += "_%d"
    if nProcessors is None:
        nProcessors = utilities.getCpuCount()
    nProcessors = max(1, nProcessors)
    print("Standarizing energy with %d processors" % nProcessors)
    epochs = utilities.get_epoch_folders(folder)
    files = []
    if not epochs:
        # path does not contain an adaptive simulation, we'll try to retrieve
        # trajectories from the specified path
        files = analysis_utils.process_folder(None, folder, trajName, reportName, os.path.join(folder, outputFilename), None, trajs_to_select)
    for epoch in epochs:
        print("Epoch", epoch)
        files.extend(analysis_utils.process_folder(epoch, folder, trajName, reportName, os.path.join(folder, epoch, outputFilename), None, trajs_to_select))
    pool = mp.Pool(nProcessors)
    results = [pool.apply_async(process_file, args=(info[1], info[4], format_out, new_report, info[3], col_energy)) for info in files]
    pool.close()
    pool.join()
    for res in results:
        res.get()

if __name__ == "__main__":
    energy_col, path, out_name, fmt_str, n_proc, out_folder, new_reports, name_report, traj_filter = parseArguments()
    if traj_filter is not None:
        traj_filter = set(traj_filter)
    main(energy_col, path, out_name, fmt_str, n_proc, out_folder, new_reports, name_report, traj_filter)
