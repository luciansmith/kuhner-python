
# collate ASCAT data from all challenge patients; optimized for use of ncpus CPUs
# (change ncpus to run on Mac)
# assumes existance of directories run1, run2, ... runNCPUS

# DEBUG ncpus = 6
ncpus = 6

import os.path
import p_rawsegs_module   # routine to collate one ascat rawsegs file

samples = []
for line in open("samples","r"):
  samples.append(line.rstrip())

# divide samples into patients
pats = {}
for sample in samples:
  pid = sample.split("_")[0]
  if pid in pats:
    pats[pid].append(sample)
  else:
    pats[pid] = [sample,]

# REMOVE all patients with less than 3 samples
for key in pats.keys():
  if len(pats[key]) < 3:
    print "Deleting patient",key,"due to having only",len(pats[key]),"samples"
    del pats[key]

pids = pats.keys()
pids.sort()

npat = len(pids)
if npat < ncpus:
  print "ERROR:  trying to divide",npat,"runs across",ncpus,"processors"
  exit()

nperdir = int(float(npat)/ncpus) + 1     # +1 handles the remainder

startpat = 0
endpat = nperdir

# create outfile and write header to it
outfile = open("all_lesions.txt","w")
outline = "patient\tbiopsy\tchrom\tsegstart\tsegend\t"
outline += "rawA\trawB\tintA\tintB\n"
outfile.write(outline)

# make probe location dictionary (done just once for speed!)
markerlocations = p_rawsegs_module.makeprobedict()

# collate data from each file into "outfile"
for n in range(1,ncpus+1):
  dirname = "run"+str(n)+"/"
  # locate each rawsegs file
  for x in range(startpat, endpat):
    mypid = pids[x]
    for sample in pats[mypid]:
      pid, sid = sample.split("_")[0:2]
      print "collating",dirname,pid,sid
      goodresult = p_rawsegs_module.collate(dirname, pid, sid, outfile, 
        markerlocations)
      if not goodresult:
        print "P-ASCAT failed on",pid,sid
  startpat = endpat
  endpat = min(endpat + nperdir,npat)

outfile.close()
