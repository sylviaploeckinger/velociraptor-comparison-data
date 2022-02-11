#! /usr/bin/env python3
"""
Usage:
  convert.py [--nproc NPROC]

where NPROC is the number of parallel processes used to run conversion scripts.

This script will run all the conversion scripts found in the data directory:
  data/*/conversion/*.py
Every script is run as follows:
  python3 $script_name ../../../cosmology.py
using the script's conversion directory as working directory.

To allow for efficient parallel execution, the script first creates a list of
all the scripts that need to be run. It then creates NPROC slots and launches
a subprocess running a script in each slot (or until there are no more scripts
to run). As soon as a subprocess slot finishes, a next script is launched in
the same slot, until all scripts have been processed.

Based on an earlier bash script that processed scripts serially.
"""

import subprocess
import glob
import os
import argparse

# Parse the optional command line argument
argparser = argparse.ArgumentParser()
argparser.add_argument("--nproc", "-n", type=int, default=1)
args = argparser.parse_args()

# List all the scripts
conversion_scripts = sorted(glob.glob("data/*/conversion/*.py"))
# Convert every script into a (command, working directory) pair
cmds = []
for script_path in conversion_scripts:
    wdir, script = os.path.split(script_path)
    cmd = f"python3 {script} ../../../cosmology.py"
    cmds.append((cmd, wdir))

# Create NPROC empty slots
slots = [None] * args.nproc
# Fill the slots with the first NPROC scripts
icmd = 0
while icmd < len(cmds) and icmd < args.nproc:
    cmd, wdir = cmds[icmd]
    # setting 'cwd' runs the script as if it was run from that directory
    slots[icmd] = (cmd, wdir, subprocess.Popen(cmd, cwd=wdir, shell=True))
    icmd += 1

# Keep track of failed processes. As soon as one process fails, we need to
# exit the script with a non-zero return code
general_return_code = 0
# Keep processing slots until all commands have been processed
while icmd < len(cmds):
    # Loop over slots
    for islot in range(args.nproc):
        cmd, wdir, handle = slots[islot]
        # poll() is None as long as the subprocess is running
        if not handle.poll() is None:
            # Check the return code of the process
            if handle.returncode != 0:
                print(
                    f'Script {wdir} -> "{cmd}" failed.\nReturn code {handle.returncode}.'
                )
            # Update the general return code
            general_return_code = max(general_return_code, handle.returncode)
            # Launch the next script
            cmd, wdir = cmds[icmd]
            slots[islot] = (cmd, wdir, subprocess.Popen(cmd, cwd=wdir, shell=True))
            icmd += 1
            # Jump out of the loop if there are no more processes
            if icmd == len(cmds):
                break

# Now wait for the final subprocesses to come home
for islot in range(args.nproc):
    if not slots[islot] is None:
        cmd, wdir, handle = slots[islot]
        # Do a blocking wait(), since we have nothing else to do
        handle.wait()
        if handle.returncode != 0:
            print(f'Script {wdir} -> "{cmd}" failed.\nReturn code {handle.returncode}.')
        # Update the general return code
        general_return_code = max(general_return_code, handle.returncode)

exit(general_return_code)
