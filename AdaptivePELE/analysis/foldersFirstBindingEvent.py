import argparse
import os
import findfirstbindingevent
import analyse

def parseArguments():
    """
        Parse the command-line options

        :returns: list, int, float, int, bool -- List of folders, column with
            binding event related metric, threshold for a binding event to be
            considered, number of steps per epoch to be consdidered, wether the
            simulation to analyse is and adaptive or sequential simulation
    """
    desc = "Program that computes the first binding event for a series of adaptive sampling runs"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("column", type=int, help="Column with binding event related metric")
    parser.add_argument("threshold", type=float, help="Threshold for a binding event to be considered")
    parser.add_argument("stepsPerEpoch", type=int, help="StepsPerEpoch")
    parser.add_argument("-seq", action='store_true', help="Use a sequential run, instead of adaptive")
    parser.add_argument("folders", nargs='+', default=".",  help="Folders with adaptive sampling runs")
    args = parser.parse_args()

    return  args.folders, args.column, args.threshold, args.stepsPerEpoch, args.seq


def main(folders, column, threshold, stepsPerEpoch, sequential):
    """
        Calculate first binding event statistics (mean, median, std)

        :param folders: List of folders
        :type folders: list
        :param column: Column with binding event related metric
        :type column: int
        :param threshold: Threshold for a binding event to be considered
        :type threshold: float
        :param stepsPerEpoch: Number of steps per epoch to be consdidered
        :type stepsPerEpoch: int
        :param sequential: Whether the simulation to analyse is and adaptive or
            sequential simulation
        :type sequential: bool
    """
    firstBE = []
    for folder in folders:
        cwd = os.getcwd()
        os.chdir(folder)
        if sequential:
            stepsFirstBindingEvent = findfirstbindingevent.findEpochFirstBindingEvent(threshold, column)
        else:
            stepsFirstBindingEvent = findfirstbindingevent.findFirstBindingEvent(stepsPerEpoch, column, threshold)

        if not stepsFirstBindingEvent:
            print "Didn't see any first binding event in folder:", folder
        else:
            firstBE.append(stepsFirstBindingEvent)

        os.chdir(cwd)

    if len(firstBE) > 1:
        print firstBE
        analyse.analyseData(firstBE)
    else:
        print firstBE[0]


if __name__ == "__main__":
    folders, column, threshold, stepsPerEpoch, seq = parseArguments()
    main(folders, column, threshold, stepsPerEpoch, seq)
