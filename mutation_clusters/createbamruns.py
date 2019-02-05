#! /usr/bin/env python2

# This program takes in a list of S3 "paths" and creates, for each entry
# in that list, a subdirectory in which to do a bamscore analysis of the
# corresponding BAM file.  It copies bamscore.py into the directory and
# writes a runbam.sh script with the correct inputs.  It does NOT actually
# run bamscore.py.

# This will make a directory structure starting where you run it.  It does
# DESTROY a previous structure of that name; if you want to run it
# again, create a sandbox each time!

bybranch = True

rundir = "bamruns"

import os
import sys

pathfile = "curated_bamlist"


if bybranch:
    rundir += "_bybranch/"
else:
    rundir += "/"

if len(sys.argv) == 2:
    pathfile = sys.argv[1]

elif len(sys.argv) > 2:
  print("USAGE:  python createbamruns.py pathfile")
  exit()


# read path list
paths = []
for line in open(pathfile,"r"):
  paths.append(line.rstrip())

if os.path.isdir(rundir):
    print("Please move or remove the directory", rundir, "before running this program: it's designed to create it from scratch.")
    exit()

os.mkdir(rundir)

shname = "run_bams.sh"
if bybranch:
    shname = "run_bams_bybranch.sh"
runall = open(shname, "w")


for path in paths:
  items = path.split("/")
  items = items[-1].split("-")
  print("path = ",path)
  pid = items[0]
  sid = items[1]
  print("pid", pid, "sid", sid)
  myname = pid + "-" + sid
  dirname = rundir + "run" + myname
  if not(os.path.isdir(dirname)):
      os.mkdir(dirname)
  outname = rundir + "run" + myname + "/runbam.sh"
  outfile = open(outname,"w")
  outline = "#!/bin/bash\n"
  outfile.write(outline)
  outline = "ml Python/2.7.14-foss-2016b-fh1\n"
  outfile.write(outline)
  outline = "sbatch -M beagle -n 4 bamscore.py "
  outline += "s3://fh-pi-reid-b/"
  outline += path + "\n"
  outfile.write(outline)
  outfile.close()
  os.system("chmod a+x " + outname)
  os.system("cp bamscore.py " + dirname)
  print("processed",pid,sid)
  runall.write("./" + dirname + "/runbam.sh\n")
print("All setup completed")
runall.close()

