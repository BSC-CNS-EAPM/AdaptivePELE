import os
import shutil 

def cleanup(tmpFolder):
    if os.path.exists(tmpFolder):
        shutil.rmtree(tmpFolder)

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
