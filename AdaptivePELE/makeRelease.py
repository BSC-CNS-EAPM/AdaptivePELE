"""
    Only for developers.
    Change releaseName to build a new release.
"""

import glob
import shutil
import os

releaseName = "v1.4"
# releaseFolder = "/gpfs/projects/bsc72/adaptiveSampling/bin"
releaseFolder = "/data2/bsc72/AdaptiveSampling/bin"
toOmit = ["tests", "runAllTests.py", "os", "sys", "TODO.txt", "Data", "Documents", "DataLocal", "epsilon_values.txt", "makeRelease.py", ".git", ".gitignore"]


files = glob.glob("*")

destFolder = os.path.join(releaseFolder, releaseName, "AdaptivePELE", "%s")
for filename in files:
    if filename in toOmit or filename.startswith(".") or filename.endswith("pyc"):
        continue
    print "Copying", filename
    try:
        if not os.path.exists(destFolder % filename):
            shutil.copytree(filename, destFolder % filename)
    except OSError:
        shutil.copyfile(filename, destFolder % filename)

print "Done with release %s!" % releaseName
