import os
import argparse
import glob
import types



def parseArguments():
    desc = "Generates output for gnuplot\n"\
            "It MUST be run from the root epoch folder (i.e., where it can find the folders 0/, 1/, 2/, ... lastEpoch/"\
            "To be run for example like: \">python generateGnuplotFile.py | gnuplot -persist\""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("steps", type=int, default=25, help="Pele steps per run") 
    args = parser.parse_args()
    return args.steps

def generateNestedString(gnuplotString, reportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps=False, replotFirst=False):
    """
        TotalNumberOfSteps -> not only considering steps in current epoch, but steps with all previous epochs
    """
    allFolders = os.listdir('.')
    epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]
    numberOfEpochs=int(len(epochFolders))


    dictionary = {'reportName':reportName, 'col2':column2, 'numberOfEpochs':numberOfEpochs, 'withLines':''}
    
    #runs of epoch 0, assumed constant
    numberOfRunsPerEpoch=len( glob.glob( os.path.join(str(0), reportName+"*") ) )
    dictionary['runsPerEpoch'] = numberOfRunsPerEpoch


    if printWithLines:
        dictionary['withLines'] = "w l"

    if type(column1) == types.IntType:
        if totalNumberOfSteps:
            dictionary['col1'] = "($" + str(column1) + "+ (%d*j))"%stepsPerRun #adds steps per runs, so that it mathes the total number of steps
        else:
            dictionary['col1'] = str(column1)
    elif type(column1) == types.LambdaType:
        dictionary['col1'] = "(" + str(column1(i)) + ")"

    print gnuplotString % dictionary + "\n"


def generateForLoopString(gnuplotString, reportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps=False, replotFirst=False):
    """
        TotalNumberOfSteps -> not only considering steps in current epoch, but steps with all previous epochs
    """
    allFolders = os.listdir('.')
    epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]
    numberOfEpochs=int(len(epochFolders))


    dictionary = {'reportName':reportName, 'col2':column2, 'numberOfEpochs':numberOfEpochs, 'withLines':''}
    plottingString = ""
    for i in range(numberOfEpochs):
        dictionary['epoch'] = i

        numberOfRunsPerEpoch=len( glob.glob( os.path.join(str(i), reportName+"*") ) )
        dictionary['runsPerEpoch'] = numberOfRunsPerEpoch

        if printWithLines:
            dictionary['withLines'] = "w l"

        if type(column1) == types.IntType:
            if totalNumberOfSteps:
                dictionary['col1'] = "($" + str(column1) + "+" + str(stepsPerRun*i) + ")" #adds steps per runs, so that it mathes the total number of steps
            else:
                dictionary['col1'] = str(column1)
        elif type(column1) == types.LambdaType:
            dictionary['col1'] = "(" + str(column1(i)) + ")"

        if i != 0 or replotFirst:
            plottingString += "re"

        plottingString += gnuplotString % dictionary + "\n"
        #plottingString += "pause 0.25\n"

    print plottingString



def main():
    stepsPerRun = parseArguments()


    #VARIABLES TO SET WHEN PRINTING
    kindOfPrint = "PRINT_RMSD_STEPS"
    #kindOfPrint = "PRINT_BE_RMSD"

    ownClustering = True

    if kindOfPrint == "PRINT_RMSD_STEPS":
        printWithLines = True
        column1 = 2 #steps
        #column2 = 7 #rmsd

        column2 = 8 #rmsd

        #PR
        #column2 = 7 #sasa

        totalNumberOfSteps=True
    elif kindOfPrint == "PRINT_BE_RMSD":
        printWithLines = False
        column1 = 7 #rmsd
        column2 = 6 #binding energy

        column1 = 5 #rmsd
        column2 = 6 #binding energy

        column1 = 8 #rmsd
        column2 = 7 #binding energy

        #PR
        #column1 = 7 #sasa
        #column2 = 5 #binding energy
        totalNumberOfSteps=False

    if ownClustering:
        representativeReportName='clustering/reports'
    else:
        representativeReportName='pyproct/reports'



    reportName='run_report_'
    #reportName='report_'
    """
    gnuplotString = "plot for [i=1:%(runsPerEpoch)d] \"%(epoch)d/%(reportName)s\".i u %(col1)s:%(col2)d lt 6 lc palette frac %(epoch)d/%(numberOfEpochs)d. notitle %(withLines)s"
    generateForLoopString(gnuplotString, reportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps, False)


    """
    gnuplotString = "plot for [j=0:%(numberOfEpochs)d-1] for [i=1:%(runsPerEpoch)d] \"\".j.\"/%(reportName)s\".i u %(col1)s:%(col2)d lt 6 lc palette frac j/%(numberOfEpochs)d. notitle %(withLines)s"
    generateNestedString(gnuplotString, reportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps, False)


    printWithLines = False
    
    if kindOfPrint == "PRINT_RMSD_STEPS":
        column1 = lambda x: (x+1)*stepsPerRun

    #Comment if you don't want to plot representative info
    #refactor
    gnuplotString = "plot for [j=0:%(numberOfEpochs)d-1] \"\".j.\"/%(reportName)s\" u %(col1)s:%(col2)d lt 6 lc palette frac j/%(numberOfEpochs)d. notitle %(withLines)s"
    #generateForLoopString(gnuplotString, representativeReportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps, replotFirst=True)

    #os.system(plottingString + "\| gnuplot -persist")

if __name__ == "__main__":
    main()
