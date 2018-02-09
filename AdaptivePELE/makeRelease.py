"""
    Only for developers.
    Change releaseName to build a new release.
"""
import socket
import glob
import shutil
import os

def copy_ignore(src, names):
    return filter(lambda x: x.endswith(".c") or x.endswith(".so"), names)


machine = socket.gethostname()
if "bsccv" in machine:
    releaseFolder = "/data2/bsc72/AdaptiveSampling/bin"
elif 'login' in machine:
    name = os.getenv("BSC_MACHINE")
    if name == "mn4":
        releaseFolder = "/gpfs/projects/bsc72/adaptiveSampling/bin"
    elif name == "nord3":
        releaseFolder = "/gpfs/projects/bsc72/adaptiveSampling/bin_nord"

releaseName = "v1.4.1"
toOmit = ["tests", "runAllTests.py", "os", "sys", "TODO.txt", "Data", "Documents", "DataLocal", "epsilon_values.txt", "makeRelease.py", ".git", ".gitignore"]


files = glob.glob("*")

destFolder = os.path.join(releaseFolder, releaseName, "AdaptivePELE", "%s")
for filename in files:
    if filename in toOmit or filename.startswith(".") or filename.endswith("pyc"):
        continue
    print "Copying", filename
    try:
        if not os.path.exists(destFolder % filename):
            shutil.copytree(filename, destFolder % filename, ignore=copy_ignore)
    except (IOError, OSError):
        shutil.copyfile(filename, destFolder % filename)

extraFiles = ["../README.rst", "../setup.py"]
for filename in extraFiles:
    if not os.path.exists(destFolder % filename):
        shutil.copyfile(filename, destFolder % filename)
        print "Copying", os.path.split(filename)[1]
print "Done with release %s!" % releaseName
