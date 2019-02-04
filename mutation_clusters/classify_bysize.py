# classify.py                    Mary Kuhner and Jon Yamato 2018/04/02

# This program receives the output of bamscore.py 
# to answer the question "For two mutations close enough together
# that they might be in the same read, how often are they actually both
# present on a read spanning both their positions?" 

# NOTE:  BAM is zero based.  VCF is one based.  Subtract 1 from your VCF
# mutation position before looking it up in BAM!

# WARNING:  We wrote this code against a BAM file where paired reads had the
# same name, but this is not standardized:  it will NOT WORK RIGHT if paired
# reads append /1, /2 or _1, _2 or _a, _b or anything like that!

import os
from os import path
from os import mkdir
import matplotlib.pyplot as plt
import lucianSNPLibrary as lsl

read_cutoff = 5

# how many non ref/alt bases can we tolerate at a mutation pair
# before throwing it out?
max_unexpected_bases = 5

# lowest acceptable mapping quality for a read
min_quality = 25

resultdir = "results/"
outdir = "analysis/"


if not(path.isdir(outdir)):
    mkdir(outdir)

# functions

def hasbin(result,bin):
  if result[bin] > read_cutoff:
    return True
  return False

def score_read(base1,base2,mut1,mut2,scorearray):
  chr1,pos1,ref1,alt1 = mut1
  chr2,pos2,ref2,alt2 = mut2

  # score the mutation
  if base1 == ref1 and base2 == ref2:  
    scorearray[0] += 1
  if base1 == ref1 and base2 == alt2:  
    scorearray[1] += 1
  if base1 == alt1 and base2 == ref2:  
    scorearray[2] += 1
  if base1 == alt1 and base2 == alt2:  
    scorearray[3] += 1
  scorearray[4] += 1

def summarize(pairlist, cisdists, transdists, nesteddists):
# taxonomy:
#   less than 12 total reads:  toosmall
#   just bin0:  wt
#   just bin1 or just bin2:  missed
#   all three bins:  fourgamete
#   bins 1 and 2 but not 3:  trans
#   bin 3 alone:  cis
#   bin 3 with bin 1 or 2 but not both:  nested
  results = {}
  results["noreads"] = 0
  results["anomaly"] = 0
  results["toosmall"] = 0
  results["wt"] = 0
  results["missed"] = 0
  results["fourgamete"] = 0
  results["cis"] = 0
  results["trans"] = 0
  results["nested"] = 0

  for tally in pairlist:

    pos1 = int(tally[6])
    pos2 = int(tally[7])
    dist = pos2 - pos1
    # noreads
    if tally[4] == 0:
      results["noreads"] += 1
      continue

    # anomaly
    unexpected_bases = tally[4] - sum(tally[0:4])
    if unexpected_bases > max_unexpected_bases:
      results["anomaly"] += 1
      continue

    # toosmall
    if tally[4] < 12:
      results["toosmall"] += 1
      continue

    # wt
    if not hasbin(tally,1) and not hasbin(tally,2) and not hasbin(tally,3):
      results["wt"] += 1
      continue
    
    # missed
    if hasbin(tally,1) and not hasbin(tally,2) and not hasbin(tally,3):
      results["missed"] += 1
      continue
    if not hasbin(tally,1) and hasbin(tally,2) and not hasbin(tally,3):
      results["missed"] += 1
      continue

    # fourgamete
    if hasbin(tally,1) and hasbin(tally,2) and hasbin(tally,3):
      results["fourgamete"] += 1
      continue

    # trans
    if hasbin(tally,1) and hasbin(tally,2) and not hasbin(tally,3):
      results["trans"] += 1
      transdists.append(dist)
      continue

    # cis
    if not hasbin(tally,1) and not hasbin(tally,2) and hasbin(tally,3):
      results["cis"] += 1
      cisdists.append(dist)
      continue

    # nested
    if (hasbin(tally,1) or hasbin(tally,2)) and hasbin(tally,3):
      results["nested"] += 1
      nesteddists.append(dist)
      continue

    print("found anomaly:",tally)
    assert False

  return results


########################################################################
# main program

rfiles = []
for root, dirs, files in os.walk(resultdir):
  for file in files:
    if file.endswith("_results.txt"):
      rfiles.append(file)

cisdists = []
transdists = []
nesteddists = []
for rfile in rfiles:
    data = []
    (pid, sid, A, B) = rfile.split("_")[0:4]
    if not(A=='1' and B=='1'):
        continue
    for line in open(resultdir + rfile,"r"):
      line = line.rstrip().split()
    
      datum = []
      for entry in line:
        entry = int(entry)
        datum.append(entry)
      data.append(datum)
    
    # classify by distances
    
    ############################################################
    # write reports
    
    results = summarize(data, cisdists, transdists, nesteddists)

lsl.createPrintAndSaveHistogram(cisdists, "Cis Distances", 0.5, xdata="distance", axis=(-20, 500, 0))
lsl.createPrintAndSaveHistogram(transdists, "Trans Distances", 0.5, xdata="distance", axis=(-20, 500, 0))
lsl.createPrintAndSaveHistogram(nesteddists, "Nested Distances", 0.5, xdata="distance", axis=(-20, 500, 0))

plt.ylim(0, 400)
plt.hist(cisdists, 100)
plt.hist(transdists, 100)
plt.hist(nesteddists, 100)
plt.show()
plt.close()
