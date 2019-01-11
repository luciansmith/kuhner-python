# mutpairhistogram.py           Mary Kuhner and Jon Yamato 2018/10/12

# This program receives the output of bamscore.py 

# NOTE:  BAM is zero based.  VCF is one based.  Subtract 1 from your VCF
# mutation position before looking it up in BAM!

# WARNING:  We wrote this code against a BAM file where paired reads had the
# same name, but this is not standardized:  it will NOT WORK RIGHT if paired
# reads append /1, /2 or _1, _2 or _a, _b or anything like that!

# This version reads peak numbers from tmscore results file and can
# stratify by peaks.

import pysam

read_cutoff = 5

# how many non ref/alt bases can we tolerate at a mutation pair
# before throwing it out?
max_unexpected_bases = 5

# lowest acceptable mapping quality for a read
min_quality = 25

# functions

def hasbin(result,bin):
  if result[bin] > read_cutoff:
    return True
  return False

def summarize(pairlist):
# taxonomy:
#   less than 12 total reads:  toosmall
#   just bin0:  wt
#   just bin1 or just bin2:  missed
#   all three bins:  fourgamete
#   bins 1 and 2 but not 3:  trans
#   bin 3 alone:  cis
#   bin 3 with bin 1 or 2 but not both:  nested
  results = {}
  results["noreads"] = []
  results["anomaly"] = []
  results["toosmall"] = []
  results["wt"] = []
  results["missed"] = []
  results["fourgamete"] = []
  results["cis"] = []
  results["trans"] = []
  results["nested"] = []
  results["nocategory"] = []

  for tally in pairlist:

    distance = tally[9] - tally[8]

    # noreads
    if tally[4] == 0:
      results["noreads"].append(distance)
      continue

    # anomaly
    unexpected_bases = tally[4] - sum(tally[0:4])
    if unexpected_bases > max_unexpected_bases:
      results["anomaly"].append(distance)
      continue

    # toosmall
    if tally[4] < 12:
      results["toosmall"].append(distance)
      continue

    # wt
    if not hasbin(tally,1) and not hasbin(tally,2) and not hasbin(tally,3):
      results["wt"].append(distance)
      continue
    
    # missed
    if hasbin(tally,1) and not hasbin(tally,2) and not hasbin(tally,3):
      results["missed"].append(distance)
      continue
    if not hasbin(tally,1) and hasbin(tally,2) and not hasbin(tally,3):
      results["missed"].append(distance)
      continue

    # fourgamete
    if hasbin(tally,1) and hasbin(tally,2) and hasbin(tally,3):
      results["fourgamete"].append(distance)
      continue

    # trans
    if hasbin(tally,1) and hasbin(tally,2) and not hasbin(tally,3):
      results["trans"].append(distance)
      continue

    # cis
    if not hasbin(tally,1) and not hasbin(tally,2) and hasbin(tally,3):
      results["cis"].append(distance)
      continue

    # nested
    if (hasbin(tally,1) or hasbin(tally,2)) and hasbin(tally,3):
      results["nested"].append(distance)
      continue

    results["nocategory"].append(distance)

  return results

########################################################################
# main program

import sys
import os

if len(sys.argv) != 2:
  print("USAGE:  classify.py root_resultsdir")
  exit()

resultshere = sys.argv[1]

resultfiles = []
pidsused = []
for rt,ds,fs in os.walk(resultshere):
  for dr in ds:
    if dr.startswith("run"):
      pid = dr.split("-")[0].split("n")[1]
      if pid in pidsused:
        continue
      wdir = resultshere+"/"+dr
      for root,dirs,files in os.walk(wdir):
        for file in files:
          if file.endswith("results.txt"):
            resultfiles.append(wdir+"/"+file)
            pidsused.append(pid)

print("found",len(resultfiles),"patients")

allresults = {}
for resultfile in resultfiles:
  items = resultfile.split("/")
  items = items[-1].split("-")
  pid = items[0]
  sid = items[1]

  data = []
  peaks = []
  for line in open(resultfile,"r"):
    line = line.rstrip().split()
  
    # process peak header
    if line[0] == "Peaks:":   # header giving peak data
      for entry in line[1:]:
        vaf,peakstart,peakend = entry.split(",")
        vaf = float(vaf[1:])   # lose starting paren
        peakstart = float(peakstart)
        peakend = float(peakend[:-1])
        peaks.append([vaf,peakstart,peakend])  # lose ending paren
      numpeaks = len(peaks)
      continue
        
    # otherwise, process mutation entry
    datum = []
    for entry in line:
      entry = int(entry)
      datum.append(entry)
    data.append(datum)
  
  # classify by peaks
  peakdata = [[[] for x in range(numpeaks)] for x in range(numpeaks)]
  for datum in data:
    mypeak = [datum[5],datum[6]]
    mypeak.sort()
    peakdata[mypeak[0]][mypeak[1]].append(datum)
  
  # classify by distances
  results = summarize(data)
  for key in results:
    if key in allresults:
      allresults[key].extend(results[key])
    else:
      allresults[key] = results[key]

for key in allresults:
  print(key, len(allresults[key]))

import matplotlib.pyplot as plt

cis = allresults["cis"][:]
trans_nest = allresults["trans"][:]
trans_nest.extend(allresults["nested"])

cis_by_size = [0.0 for x in range(500)]
trans_by_size = [0.0 for x in range(500)]

for item in cis:
  cis_by_size[item] += 1.0

for item in trans_nest:
  trans_by_size[item] += 1.0

distances = []
ratios = []

for x in range(500):
  if cis_by_size[x] + trans_by_size[x] > 100:
    ratio = cis_by_size[x]/(cis_by_size[x] + trans_by_size[x])
    distances.append(x)
    ratios.append(ratio)

plt.figure(0)
plt.plot(distances,ratios,"ko")

plt.figure(1)
distances = []
ratios = []
for x in range(0,500,10):

  cissum = sum(cis_by_size[x:x+10])
  transsum = sum(trans_by_size[x:x+10])
  if cissum + transsum > 100:
    ratio = cissum / (cissum + transsum)
    distances.append(x)
    ratios.append(ratio)
    
plt.plot(distances,ratios,"ko")


#n_bins = 20

#plt.figure(1)
#plt.hist(cis,n_bins)
#plt.title("Cis")

#plt.figure(2)
#plt.hist(trans_nest,n_bins)
#plt.title("Trans + nested")
plt.show()

