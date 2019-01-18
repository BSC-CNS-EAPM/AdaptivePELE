from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import argparse
import numpy as np
import mdtraj as md
from AdaptivePELE.utilities import utilities
from AdaptivePELE.analysis import correctRMSD
from AdaptivePELE.freeEnergies.extractCoords import getTopologyObject


def parseArguments():
    """
        Parse the command-line options

        :returns: object -- Object containing the options passed
    """
    desc = "Program that caculates the relative SASA of a ligand."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("resname", type=str, help="Ligand resname")
    parser.add_argument("--path", type=str, default=".", help="Path where the simulation is stored")
    parser.add_argument("--top", type=str, default=None, help="Topology file for non-pdb trajectories or path to Adaptive topology object")
    parser.add_argument("--out_name", type=str, default="fixedReport", help="Name of the modified report files (default is fixedReport)")
    parser.add_argument("--fmt_str", type=str, default="%.4f", help="Format of the output file (default is .4f which means all floats with 4 decimal points)")
    args = parser.parse_args()

    return args.resname, args.path, args.top, args.out_name, args.fmt_str


def calculateSASA(trajectory, topology, res_name):
    """
        Calculate the SASA of a ligand in a trajectory

        :param trajectory: Name of the trajectory file
        :type trajectory: str
        :param topology: Topology of the trajectory (needed for non-pdb trajs)
        :type topology: str
        :param res_name: Ligand resname
        :type res_name: str
    """
    t = md.load(trajectory, top=topology)
    res_atoms = t.top.select("resname '%s'" % res_name)
    t2 = t.atom_slice(res_atoms)
    for atom in t2.top.atoms:
        # mdtraj complains if the ligand residue index is not 0 when
        # isolated
        atom.residue.index = 0
    res_index = t.top.atom(res_atoms[0]).residue.index
    sasa = md.shrake_rupley(t, mode="residue")
    sasa_empty = md.shrake_rupley(t2, mode="residue")
    return sasa[:, res_index]/sasa_empty[:, 0]


def main(resname, folder, top, out_report_name, format_out):
    """
        Calculate the relative SASA values of the ligand

        :param resname: Ligand resname
        :type resname: str
        :param folder: Path the simulation
        :type folder: str
        :param top: Path to the topology
        :type top: str
        :param out_report_name: Name of the output file
        :type out_report_name: str
        :param format_out: String with the format of the output
        :type format_out: str
    """
    # Constants
    outputFilename = "_".join([out_report_name, "%d"])
    trajName = "*traj*"
    reportName = "*report*_%d"

    epochs = utilities.get_epoch_folders(folder)
    if top is not None:
        top_obj = getTopologyObject(top)

    for epoch in epochs:
        print("Epoch", epoch)
        os.chdir(epoch)
        allTrajs = glob.glob(trajName)

        for traj in allTrajs:
            trajNum = utilities.getTrajNum(traj)
            if top is not None:
                top_file = os.path.abspath(top_obj.getTopologyFile(int(epoch), trajNum))
            else:
                top_file = None
            sasa_values = calculateSASA(traj, top_file, resname)
            try:
                reportFilename = glob.glob(reportName % trajNum)[0]
            except IndexError:
                raise IndexError("File %s not found in folder %s" % (reportName % trajNum, epoch))

            with open(reportFilename) as f:
                header = f.readline().rstrip()
                if not header.startswith("#"):
                    header = ""
                reportFile = utilities.loadtxtfile(f)

            fixedReport = correctRMSD.extendReportWithRmsd(reportFile, sasa_values)

            with open(outputFilename % trajNum, "w") as fw:
                if header:
                    fw.write("%s\tSASA\n" % header)
                np.savetxt(fw, fixedReport, fmt=format_out, delimiter="\t")

        os.chdir("..")

if __name__ == "__main__":
    lig_name, path, topology_path, out_name, fmt_str = parseArguments()
    main(lig_name, path, topology_path, out_name, fmt_str)
