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

import pickle
import os
from os import path
from os import mkdir

read_cutoff = 5

# how many non ref/alt bases can we tolerate at a mutation pair
# before throwing it out?
max_unexpected_bases = 5

# lowest acceptable mapping quality for a read
min_quality = 25

resultdir = "results/"
bybranchdir = "results_bybranch/"
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

def classify(tally):
# taxonomy:
#   less than 12 total reads:  toosmall
#   just bin0:  wt
#   just bin1 or just bin2:  missed
#   all three bins:  fourgamete
#   bins 1 and 2 but not 3:  trans
#   bin 3 alone:  cis
#   bin 3 with bin 1 or 2 but not both:  nested

    # noreads
    if tally[4] == 0:
        return "noreads"

    # anomaly
    unexpected_bases = tally[4] - sum(tally[0:4])
    if unexpected_bases > max_unexpected_bases:
        return "anomaly"

    # toosmall
    if tally[4] < 12:
        return "toosmall"

    # wt
    if not hasbin(tally,1) and not hasbin(tally,2) and not hasbin(tally,3):
        return "wt"
    
    # missed
    if hasbin(tally,1) and not hasbin(tally,2) and not hasbin(tally,3):
        return "missed"
    if not hasbin(tally,1) and hasbin(tally,2) and not hasbin(tally,3):
        return "missed"

    # fourgamete
    if hasbin(tally,1) and hasbin(tally,2) and hasbin(tally,3):
        return "fourgamete"

    # trans
    if hasbin(tally,1) and hasbin(tally,2) and not hasbin(tally,3):
        return "trans"

    # cis
    if not hasbin(tally,1) and not hasbin(tally,2) and hasbin(tally,3):
        return "cis"

    # nested
    if (hasbin(tally,1) or hasbin(tally,2)) and hasbin(tally,3):
        return "nested"

    return "unknown"

########################################################################
# main program

rfiles = []
for root, dirs, files in os.walk(resultdir):
  for file in files:
    if file.endswith("_results.txt"):
      rfiles.append((resultdir + file, "all"))

for root, dirs, files in os.walk(bybranchdir):
  for file in files:
    if file.endswith("_results.txt"):
      rfiles.append((bybranchdir + file, "bybranch"))


allmuts = {}
for (rfile, allorbranch) in rfiles:
    data = []
    (patient, sample, A, B) = rfile.split("/")[1].split("_")[0:4]
    AB = (A, B)
    if patient not in allmuts:
        allmuts[patient] = {}
    if sample not in allmuts[patient]:
        allmuts[patient][sample] = {}
    if AB not in allmuts[patient][sample]:
        allmuts[patient][sample][AB] = {}
    
    for line in open(rfile,"r"):
        line = line.rstrip().split()

        tally = []
        for entry in line:
            entry = int(entry)
            tally.append(entry)
        (chr, pos1, pos2)  = tally[5:8]
        dist = pos2 - pos1
        pos12 = (pos1, pos2)
        assert dist > 0
        if chr not in allmuts[patient][sample][AB]:
            allmuts[patient][sample][AB][chr] = {}
        if pos12 not in allmuts[patient][sample][AB][chr]:
            allmuts[patient][sample][AB][chr][pos12] = {}
        allmuts[patient][sample][AB][chr][pos12][allorbranch] = classify(tally)


############################################################
# write reports

pickle.dump(allmuts, open("all_mutations.bin", "wb"))

labels = "Patient"
labels += "\tSample"
labels += "\tA"
labels += "\tB"
labels += "\tpos1"
labels += "\tpos2"
labels += "\tall"
labels += "\tpos2"
labels += "\n"

outfile = open("all_mutations.tsv", "w")
outfile.write(labels)

for patient in allmuts:
    for sample in allmuts[patient]:
        for AB in allmuts[patient][sample]:
            (A, B) = AB
            for chr in allmuts[patient][sample][AB]:
                for pos12 in allmuts[patient][sample][AB][pos12]:
                    (pos1, pos2) = pos12
                    outline = patient
                    outline += "\t" + sample
                    outline += "\t" + A
                    outline += "\t" + B
                    outline += "\t" + str(pos1)
                    outline += "\t" + str(pos2)
                    if "all" in allmuts[patient][sample][AB][pos12]:
                        outline += "\t" + allmuts[patient][sample][AB][pos12]["all"]
                    else:
                        outline += "\tmissing"
                    if "bybranch" in allmuts[patient][sample][AB][pos12]:
                        outline += "\t" + allmuts[patient][sample][AB][pos12]["bybranch"]
                    else:
                        outline += "\tmissing"
                    outline += "\n"
                    outfile.write(outline)
outfile.close()