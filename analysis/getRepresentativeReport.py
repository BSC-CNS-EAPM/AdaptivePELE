import os
import argparse
import glob

OUTPUT_PATH_TEMPLETIZED = "%d"

REPORT_FILE = "*report_%d"

#cluster options
CLUSTER_CENTERS = "centers".upper()
CLUSTER_ELEMENTS = "all".upper()
CLUSTER_COMPRESSED = "compressed".upper()
OWN_CLUSTERING = "own".upper()

#folder options
PYPROCT_FOLDER = "pyproct"
PYPROCT_THIS_EPOCH_FOLDER = "pyproct_thisEpoch"
OWN_CLUSTERING_FOLDER = "clustering"

def getFileAndSnapshot(model):
    """
        The model number corresponds to the trajectory number in pele, but not with the report number (traj. starts counting with 1, whereas report starts counting with 0)
    """
    splitModel = model.split("\n")
    origFile = ""
    origSnapshot = 0
    for line in splitModel:
        if "REMARK source" in line:
            origFile = line.split(":")[1].strip()
        if "REMARK original model nr" in line:
            origSnapshot = int(line.split(":")[1])
            break
    return origFile, origSnapshot

def getPDBSnapshots(file, verbose=False):
    inputFile = open(file, "r")
    inputFileContent = inputFile.read()

    if "ENDMDL" in inputFileContent:
        snapshots = inputFileContent.split("ENDMDL")[:-1]
    else:
        snapshots = inputFileContent.split("ENDMDL")
    
    if not verbose: 
        return snapshots

    remarkInfo = "REMARK 000 File created using Prody and pyProCT\n\
                    REMARK source            : %s\n\
                    REMARK original model nr : %d\n"
    snapshotsWithInfo  = [remarkInfo%(file, i) + snapshots[i] for i in enumerate(snapshots)]
    return snapshotsWithInfo

def getReportFile(trajFilename):
    number = int(trajFilename[:-4].split("_")[-1])
    return REPORT_FILE%number

def getEpoch(trajFilename):
    """
        Expects as input XXX/epoch/trajectoryFile.pdb
        Returns epoch
    """
    return trajFilename.split("/")[-2]

def getReportFilename(pyProcTModel):
    origFile, origSnapshot = getFileAndSnapshot(pyProcTModel)
    epoch = getEpoch(origFile)
    reportFilename = getReportFile(origFile)
    return os.path.join(epoch, reportFilename), origSnapshot


def getReportLine(reportFilename, acceptedStep):
    reportFilename = glob.glob(reportFilename)[0]
    reportFile = open(reportFilename, 'r')
    return reportFile.readlines()[acceptedStep].strip('\n') #First line -> comment. However, starts counting with 0, whereas acceptedSteps starts counting in 1.



def getRepresentativesReportLine(representativesFile):
    """
        Expects the tree of folders: epochnum/simulationreports
        Ex:
            + 0/ --> trajectory.pdb
                 --> report

            + 1/ --> trajectory.pdb
                 --> report
            ...
    """
    models = getPDBSnapshots(representativesFile)

    reportLines = []
    for i, model in enumerate(models):
        reportFilename, acceptedStep = getReportFilename(model)
        reportLine = getReportLine(reportFilename, acceptedStep)
        reportLines.append({'file':reportFilename, 'reportLine':reportLine, 'acceptedStep': acceptedStep})
    return reportLines

def writeReportFile(outputReportFile, report):
    outputFile = open(outputReportFile, 'w')
    for reportLine in report:
        outputFile.write("%(reportLine)s #file: %(file)s, acceptedStep: %(acceptedStep)d\n"%reportLine)
    outputFile.close()
    print "Wrote ", outputReportFile



def parseArguments():
    desc = "Writes the report from the representatives in order to extract simulation information such as binding energy, rmsd...\n"\
            "It MUST be run from the root epoch folder (i.e., where it can find the folders 0/, 1/, 2/, ... lastEpoch/"
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-c", "--choice", default=CLUSTER_CENTERS, choices=[CLUSTER_CENTERS, CLUSTER_ELEMENTS, CLUSTER_COMPRESSED, OWN_CLUSTERING], type=str.upper, help="Either to get the report of the cluster centers or all elements")
    parser.add_argument("-f", "--folder", default=PYPROCT_FOLDER, choices=[PYPROCT_FOLDER, PYPROCT_THIS_EPOCH_FOLDER], help="")

    args = parser.parse_args()
    return args.choice, args.folder

def main():
    choice, pyproctFolder = parseArguments()

    PYPROCT_OUTPUT_DIR = os.path.join(OUTPUT_PATH_TEMPLETIZED, pyproctFolder)
    OWN_CLUSTERING_OUTPUT_DIR = os.path.join(OUTPUT_PATH_TEMPLETIZED, OWN_CLUSTERING_FOLDER)

    if choice == CLUSTER_CENTERS:
        PYPROCT_REPRESENTATIVE_OUTPUT = os.path.join(PYPROCT_OUTPUT_DIR, "results/representatives.pdb")
    elif choice == CLUSTER_ELEMENTS:
        PYPROCT_CLUSTERS_OUTPUT = os.path.join(PYPROCT_OUTPUT_DIR, "clusters/cluster_*.pdb")
        CLUSTERS_OUTPUT = os.path.join(PYPROCT_OUTPUT_DIR, "clusters/report_%d.pdb")
    elif choice == CLUSTER_COMPRESSED:
        PYPROCT_COMPRESSED_OUTPUT = os.path.join(PYPROCT_OUTPUT_DIR, "results/compressed.pdb")
    elif choice == OWN_CLUSTERING:
        OWN_CLUSTERING_REPRESENTATIVE_OUTPUT = os.path.join(OWN_CLUSTERING_OUTPUT_DIR, "cluster_*.pdb")

    allFolders = os.listdir('.')

    epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]

    numberOfEpochs=int(len(epochFolders))
    for i in range(numberOfEpochs): #Representatives are not obtained in the last round
        if choice == CLUSTER_CENTERS:
            representativesFile = PYPROCT_REPRESENTATIVE_OUTPUT%i
            reportLines = getRepresentativesReportLine(representativesFile)

            REPRESENTATIVES_REPORT=os.path.join(PYPROCT_OUTPUT_DIR, "reports")
            print reportLines
            writeReportFile(REPRESENTATIVES_REPORT%i, reportLines)

        elif choice == CLUSTER_ELEMENTS:
            clusterFiles = glob.glob(PYPROCT_CLUSTERS_OUTPUT%i)

            for clusterFile in clusterFiles:
                reportLines = getRepresentativesReportLine(clusterFile)
                clusterNumber = int(clusterFile.split("_")[-1][:-4])
                writeReportFile(CLUSTERS_OUTPUT%(i, clusterNumber), reportLines)

        elif choice == CLUSTER_COMPRESSED:
            representativesFile = PYPROCT_COMPRESSED_OUTPUT%i
            reportLines = getRepresentativesReportLine(representativesFile)

            REPRESENTATIVES_REPORT=os.path.join(PYPROCT_OUTPUT_DIR, "reports_compressed")
            writeReportFile(REPRESENTATIVES_REPORT%i, reportLines)

        elif choice == OWN_CLUSTERING:
            clusterFiles = glob.glob(OWN_CLUSTERING_REPRESENTATIVE_OUTPUT%i)

            reportLines = []
            for clusterFile in clusterFiles:
                reportLine = getRepresentativesReportLine(clusterFile)
                reportLines += reportLine

            REPRESENTATIVES_REPORT=os.path.join(OWN_CLUSTERING_OUTPUT_DIR, "reports")
            writeReportFile(REPRESENTATIVES_REPORT%i, reportLines)

if __name__ == "__main__":
    main()
