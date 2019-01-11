# classify.py                    Mary Kuhner and Jon Yamato 2018/04/02

# This program receives the output of bamscore.py 
# to answer the question "For two mutations close enough together
# that they might be in the same read, how often are they actually both
# present on a read spanning both their positions?" 

# Currently findmutations.py filters out any position that was not
# called BOTH 1/1 in diploid and 2/2 in tetraploid, if both solutions
# existed.  This means that there will not be mutation files for all
# samples.

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

def unpack(line):
  mut1 = [line[0],int(line[1]),line[2],line[3],float(line[4])]
  mut2 = [line[5],int(line[6]),line[7],line[8],float(line[9])]
  peakinfo = [int(line[10]),int(line[11])]
  return [mut1,mut2,peakinfo]

def hasbin(result,bin):
  if result[bin] > read_cutoff:
    return True
  return False

def score_read(base1,base2,mut1,mut2,scorearray):
  chr1,pos1,ref1,alt1,vaf1 = mut1
  chr2,pos2,ref2,alt2,vaf2 = mut2

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
      continue

    # cis
    if not hasbin(tally,1) and not hasbin(tally,2) and hasbin(tally,3):
      results["cis"] += 1
      continue

    # nested
    if (hasbin(tally,1) or hasbin(tally,2)) and hasbin(tally,3):
      results["nested"] += 1
      continue

    print "found anomaly:",tally
    assert False

  return results

########################################################################
# main program


import sys
if len(sys.argv) != 2:
  print "USAGE:  classify.py *-results.txt"
  exit()

resultfile = sys.argv[1]

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
peakdata = [[[] for x in xrange(numpeaks)] for x in xrange(numpeaks)]
for datum in data:
  mypeak = [datum[5],datum[6]]
  mypeak.sort()
  peakdata[mypeak[0]][mypeak[1]].append(datum)

# classify by distances
bin5 = [[] for x in xrange(101)]
bin20 = [[] for x in xrange(26)]
bin50 = [[] for x in xrange(11)]
for datum in data:
  pos1 = int(datum[8])
  pos2 = int(datum[9])
  dist = pos2 - pos1
  assert dist > 0
  ind5 = int(dist/5.0)
  bin5[ind5].append(datum)
  ind20 = int(dist/20.0)
  bin20[ind20].append(datum)
  ind50 = int(dist/50.0)
  bin50[ind50].append(datum)

############################################################
# write reports

outfilename = pid + "-" + sid + "-analysis.txt" 
outfile = open(outfilename,"w")

results = summarize(data)
scorable = results["cis"] + results["trans"] + results["nested"]
scorable = float(scorable)
cistrans = (results["cis"] + results["nested"])/scorable
cistrans = cistrans * 100.0

outline = "Overall results:\n"
outfile.write(outline)
outline = "total " + str(len(data)) + "\n"
outfile.write(outline)
outline = "cis/trans ratio " + str(cistrans) + "\n"
outfile.write(outline)
outline = "cis " + str(results["cis"]) + "\n"
outfile.write(outline)
outline = "nested " + str(results["nested"]) + "\n"
outfile.write(outline)
outline = "trans " + str(results["trans"]) + "\n"
outfile.write(outline)
outline = "missed " + str(results["missed"]) + "\n"
outfile.write(outline)
outline = "wt " + str(results["wt"]) + "\n"
outfile.write(outline)
outline = "fourgamete " + str(results["fourgamete"]) + "\n"
outfile.write(outline)
outline = "anomaly " + str(results["anomaly"]) + "\n"
outfile.write(outline)
outline = "noreads " + str(results["noreads"]) + "\n"
outfile.write(outline)
outline = "toosmall " + str(results["toosmall"]) + "\n"
outfile.write(outline)

total = 0
for categ in results.keys():
  total += results[categ]
if total != len(data):
  print "WARNING:  not all pairs accounted for in total pairs!"


# report on peaks
outline = "\nPer-peak results\n"
outfile.write(outline)
peakno = 1
for vaf,peakstart,peakend in peaks:
  outline = ""
  outline += "Peak #" + str(peakno) + ": "
  peakno += 1
  outline += "Mean VAF " + str(vaf)
  outline += ", lower bound " + str(peakstart)
  outline += ", upper bound " + str(peakend) + "\n"
  outfile.write(outline)
outfile.write("\n")

outline = "      "
for col in range(0,numpeaks):
  outline += "Peak" + str(col+1) + "      "
outline += "\n"
outfile.write(outline)
outfile.write("\n")
  
for i in range(0,numpeaks):
  outline = "Peak" + str(i+1) + " "
  for j in range(0,numpeaks):
    results = summarize(peakdata[i][j])
    nresults = len(peakdata[i][j])
    scorable = float(results["cis"] + results["nested"] + results["trans"])
    if scorable > 0:
      cistrans = 100.0 * (results["cis"] + results["nested"])/scorable
      cistrans = "{0:2.1f}".format(cistrans)
    else:
      cistrans = "  NA"
    outline += str(cistrans) + " "
    counter = "{0:2d}".format(nresults)
    outline += "[" + counter + "] "
  outline += "\n"
  outfile.write(outline)

outline = "bin5 "
for entry in bin5:
  results = summarize(entry)
  scorable = float(results["cis"] + results["nested"] + results["trans"])
  if scorable > 0:
    cistrans = 100.0 * (results["cis"] + results["nested"])/scorable
    cistrans = "{0:2.1f}".format(cistrans)
  else:
    cistrans = "  NA"
  outline +=  str(cistrans) + " "
outline += "\n"
outfile.write(outline)

outline = "bin20 "
for entry in bin20:
  results = summarize(entry)
  scorable = float(results["cis"] + results["nested"] + results["trans"])
  if scorable > 0:
    cistrans = 100.0 * (results["cis"] + results["nested"])/scorable
    cistrans = "{0:2.1f}".format(cistrans)
  else:
    cistrans = "  NA"
  outline +=  str(cistrans) + " "
outline += "\n"
outfile.write(outline)

outline = "bin50 "
for entry in bin50:
  results = summarize(entry)
  scorable = float(results["cis"] + results["nested"] + results["trans"])
  if scorable > 0:
    cistrans = 100.0 * (results["cis"] + results["nested"])/scorable
    cistrans = "{0:2.1f}".format(cistrans)
  else:
    cistrans = "  NA"
  outline +=  str(cistrans) + " "
outline += "\n"
outfile.write(outline)




outfile.close()
