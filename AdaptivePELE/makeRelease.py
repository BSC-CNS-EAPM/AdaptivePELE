"""
    Only for developers.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import socket
import shutil
import argparse
import subprocess
from datetime import datetime
import AdaptivePELE as a
from AdaptivePELE.utilities import utilities


def parseArgs():
    parser = argparse.ArgumentParser(description="Helper script to automate the process of releasing a new version")
    parser.add_argument('--name', type=str, default=None)
    arg = parser.parse_args()
    return arg.name


def copy_ignore(src, names):
    return [x for x in names if x.endswith(".c") or x.endswith(".so")]


def main(releaseName):
    machine = socket.gethostname()
    if "bsccv" in machine:
        releaseFolder = "/data2/bsc72/AdaptiveSampling/bin"
    elif 'login' in machine:
        name = os.getenv("BSC_MACHINE")
        if name == "mn4":
            releaseFolder = "/gpfs/projects/bsc72/adaptiveSampling/bin"
        elif name == "nord3":
            releaseFolder = "/gpfs/projects/bsc72/adaptiveSampling/bin_nord"
        elif name == "nvidia":
            releaseFolder = "/gpfs/projects/bsc72/adaptiveSampling/bin_mt"
        elif name == "power":
            releaseFolder = "/gpfs/projects/bsc72/adaptiveSampling/bin_cte"

    if releaseName is None:
        releaseName = "v%s" % a.__version__
    toOmit = ["tests", "runAllTests.py", "os", "sys", "TODO.txt", "Data", "Documents", "DataLocal", "epsilon_values.txt", "makeRelease.py", ".git", ".gitignore"]
    toOmit += ['runTestsCuda.sl', 'runMDTest.sl', 'runAllTests.sl', 'runAllTests_nord.sl', 'runTestsCuda_CTE.sl', 'AdaptiveTest_CUDA.err', 'AdaptiveTest_CUDA.out']


    files = glob.glob("*")
    if os.path.exists(os.path.join(releaseFolder, releaseName)):
        raise ValueError("Installation already found! Check that the version of the user-provided release name is correct")
    utilities.makeFolder(os.path.join(releaseFolder, releaseName))
    destFolder = os.path.join(releaseFolder, releaseName, "AdaptivePELE", "%s")
    for filename in files:
        if filename in toOmit or filename.startswith(".") or filename.endswith("pyc"):
            continue
        try:
            if not os.path.exists(destFolder % filename):
                print("Copying", filename)
                shutil.copytree(filename, destFolder % filename, ignore=copy_ignore)
        except (IOError, OSError):
            shutil.copyfile(filename, destFolder % filename)

    extraFiles = ["../README.rst", "../setup.py"]
    for filename in extraFiles:
        if not os.path.exists(destFolder % filename):
            shutil.copyfile(filename, destFolder % filename)
            print("Copying", os.path.split(filename)[1])

    print("Compiling cython extensions")
    os.chdir(destFolder % "..")
    subprocess.call(['python', 'setup.py', 'build_ext', '--inplace'])
    with open(os.path.join(releaseFolder, releaseName, "installation_info.txt"), "w") as fw:
        fw.write("Python used in installation: ")
        version = subprocess.check_output(["python", "--version"], universal_newlines=True)
        fw.write(version)

        fw.write("Compiler used in installation: ")
        fw.write(os.getenv("CC")+"\n")

        fw.write("Installed on %s\n" % str(datetime.now()))

        fw.write("Modules loaded in installation:\n")

        modules = subprocess.check_output("module list", universal_newlines=True, shell=True, stderr=subprocess.STDOUT)
        fw.write(modules)

    print("Done with release %s!" % releaseName)

if __name__ == "__main__":
    name_install = parseArgs()
    main(name_install)
