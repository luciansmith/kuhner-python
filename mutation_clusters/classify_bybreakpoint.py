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
import csv
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

mutation_file = "../snv_plus_indels.twoPlus.20181030.csv"


onlysomepatients = False
somepatients = ["1005"]

breakError = 10000

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
  
def readAllMuts():
    mutations = {}
    with open(mutation_file, 'r') as csvfile:
        for lvec in csv.reader(csvfile):
            if "DNANum" in lvec[0]:
                continue
            (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
            if (is_snv=="f"):
                continue
            if (is_2p=="f"):
                continue
    #        if ("N" in sample):
    #            continue
            if sample not in mutations:
                mutations[sample] = {}
            if chr not in mutations[sample]:
                mutations[sample][chr] = {}
            pos = int(pos)
            mutations[sample][chr][pos] = (ref, alt)
    return mutations

def readBreakpoints():
    breakpoints = {}
    bdir = "../gamma_test_output/pASCAT_input_g500/"
    bfiles = []
    for root, dirs, files in os.walk(bdir):
      for file in files:
        if file.endswith("copynumber_segments.txt"):
          bfiles.append(file)
    for bfilename in bfiles:
        patient = bfilename.split('_')[0]
        if onlysomepatients and patient not in somepatients:
            continue
        breakpoints[patient] = {}
        prevchr = "0"
        prevend = 0
        bfile = open(bdir+bfilename, "r")
        for line in bfile:
            if "Chr" in line:
                continue
            (chr, start, end, __, __) = line.rstrip().split()
            start = int(start)
            end = int(end)
            if (chr==prevchr):
                if chr not in breakpoints[patient]:
                    breakpoints[patient][chr] = []
                breakpoints[patient][chr].append((prevend, start))
            prevchr = chr
            prevend = end
    return breakpoints

def isNearBreak(chr, pos1, pos2, breakpoints):
    if chr not in breakpoints:
        return "Far"
    for (start, end) in breakpoints[chr]:
        if start > pos2+breakError:
            continue
        if end < pos1-breakError:
            continue
        return "Near"
    return "Far"

def summarize(pairlist, cisdists, transdists, nesteddists, breakpoints):
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
  results["neartotal"] = 0
  results["fartotal"] = 0

  for tally in pairlist:

    chr = str(tally[5])
    pos1 = int(tally[6])
    pos2 = int(tally[7])
    dist = pos2 - pos1
    nearbreak = isNearBreak(chr, pos1, pos2, breakpoints)
    for posa in (pos1, pos2):
        if isNearBreak(chr, posa, posa, breakpoints) == "Near":
            results["neartotal"] += 1
        else:
            results["fartotal"] += 1

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
      transdists[nearbreak].append(dist)
      continue

    # cis
    if not hasbin(tally,1) and not hasbin(tally,2) and hasbin(tally,3):
      results["cis"] += 1
      cisdists[nearbreak].append(dist)
      continue

    # nested
    if (hasbin(tally,1) or hasbin(tally,2)) and hasbin(tally,3):
      results["nested"] += 1
      nesteddists[nearbreak].append(dist)
      continue

    print("found anomaly:",tally)
    assert False

  return results

def getPatientSampleMap():
    patientSampleMap = {}
    samplePatientMap = {}
    callfile = open("../calling_evidence.tsv", "r")
    for line in callfile:
        if "Patient" in line:
            continue
        (patient, sample) = line.rstrip().split()[0:2]
        patientSampleMap[sample] = patient
        if patient not in samplePatientMap:
            samplePatientMap[patient] = []
        samplePatientMap[patient].append(sample)
    return patientSampleMap, samplePatientMap

########################################################################
# main program

rfiles = []
for root, dirs, files in os.walk(resultdir):
  for file in files:
    if file.endswith("_results.txt"):
      rfiles.append(file)

breakpoints = readBreakpoints()
cisdists = {}
transdists = {}
nesteddists = {}
overallnear = 0
overallfar = 0

for nb in ("Near", "Far"):
    cisdists[nb] = []
    transdists[nb] = []
    nesteddists[nb] = []

for rfile in rfiles:
    data = []
    (pid, sid, A, B) = rfile.split("_")[0:4]
    if onlysomepatients and pid not in somepatients:
        continue
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
    
    results = summarize(data, cisdists, transdists, nesteddists, breakpoints[pid])
    overallnear += results["neartotal"]
    overallfar += results["fartotal"]

for nb in ("Near", "Far"):
    print("Near break =", nb)
    if nb=="Far":
        plt.ylim(0, 400)
    plt.hist(cisdists[nb], 100)
    plt.hist(transdists[nb], 100)
    plt.hist(nesteddists[nb], 100)
    plt.show()
    plt.close()

nearfar = open("conclusions/all_distances.tsv", "w")
nearfar.write("Distance\tType\tNearFar\n")
for nb in ("Near", "Far"):
    for dist in cisdists[nb]:
        nearfar.write(str(dist));
        nearfar.write("\t" + "cis")
        nearfar.write("\t" + nb)
        nearfar.write("\n")
    for dist in transdists[nb]:
        nearfar.write(str(dist));
        nearfar.write("\t" + "trans")
        nearfar.write("\t" + nb)
        nearfar.write("\n")
    for dist in nesteddists[nb]:
        nearfar.write(str(dist));
        nearfar.write("\t" + "nested")
        nearfar.write("\t" + nb)
        nearfar.write("\n")
nearfar.close()

(patientSampleMap, samplePatientMap) = getPatientSampleMap()
mutations = readAllMuts()
allnearcount = 0
allfarcount = 0
for sample in mutations:
    for chr in mutations[sample]:
        for pos in mutations[sample][chr]:
            if isNearBreak(chr, pos, pos, breakpoints[patientSampleMap[sample]]) == "Near":
                allnearcount += 1
            else:
                allfarcount += 1

print("Allowed distance from breakpoint:", str(breakError))

print("Total 'near' mutations:", str(allnearcount))
print("Total 'far' mutations:", str(allfarcount))
print("Fraction near mutations =", str(allnearcount/(allnearcount+allfarcount)))

print("Total paired 'near' mutations:", str(overallnear))
print("Total paired 'far' mutations:", str(overallfar))
print("Fraction near mutations =", str(overallnear/(overallnear+overallfar)))


