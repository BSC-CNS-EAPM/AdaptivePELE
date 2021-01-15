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
from AdaptivePELE.constants import environments


def parseArgs():
    parser = argparse.ArgumentParser(description="Helper script to automate the process of releasing a new version")
    parser.add_argument('--name', type=str, default=None)
    arg = parser.parse_args()
    return arg.name


def copy_ignore(src, names):
    return [x for x in names if x.endswith(".so")]

def join_cmds(prefix, suffix):
    if not prefix:
        return suffix
    else:
        return prefix+"; "+suffix

def log_install(file_descriptor, prepare_str):
    file_descriptor.write("Python used in installation: ")
    version = subprocess.check_output(join_cmds(prepare_str, "python --version"), shell=True,
                                      universal_newlines=True, stderr=subprocess.STDOUT)
    file_descriptor.write(version)

    file_descriptor.write("Compiler used in installation: ")
    compiler = subprocess.check_output(join_cmds(prepare_str, "echo $CC"), shell=True,
                                       universal_newlines=True, stderr=subprocess.STDOUT)
    file_descriptor.write(compiler)

    file_descriptor.write("Installed on %s\n" % str(datetime.now()))
    file_descriptor.write("Modules loaded in installation:\n")

    modules = subprocess.check_output(join_cmds(prepare_str, "module list"), shell=True,
                                       universal_newlines=True, stderr=subprocess.STDOUT)
    file_descriptor.write(modules)

def build_extensions(name, releaseFolder, releaseName):
    all_envs = environments.get(name)
    # force recompiles everything even if no changes are detected, needed
    # to recompile with different versions
    compile_cmd = 'python setup.py build_ext --inplace --force'
    with open(os.path.join(releaseFolder, releaseName, "installation_info.txt"), "w") as fw:
        if all_envs is None:
            # if no particular environment is specified just rely on whathever is
            # set when calling the script
            subprocess.call(['python', 'setup.py', 'build_ext', '--inplace'])
            log_install(fw, "")
        else:
            for env_str in all_envs:
                prepare_str = "; ".join(["module purge 2> /dev/null", env_str])
                # call all commands in the same shell, so the module changes
                # take effect
                print(join_cmds(prepare_str, compile_cmd))
                subprocess.call(join_cmds(prepare_str, "python --version"), universal_newlines=True, shell=True)
                subprocess.call(join_cmds(prepare_str, compile_cmd), universal_newlines=True, shell=True)
                log_install(fw, prepare_str)
                fw.write("\n")

def main(releaseName):
    machine = socket.gethostname()
    name = machine
    if "bsccv" in machine:
        name = "life"
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
    build_extensions(name, releaseFolder, releaseName)

    print("Done with release %s!" % releaseName)

if __name__ == "__main__":
    name_install = parseArgs()
    main(name_install)
