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

    print("found anomaly:",tally)
    assert False

  return results

########################################################################
# main program


import sys
if len(sys.argv) != 2:
  print("USAGE:  classify.py *-results.txt")
  exit()

resultfile = sys.argv[1]

items = resultfile.split("/")
items = items[-1].split("_")
(pid, sid, A, B) = items[0:4]

data = []
for line in open(resultfile,"r"):
  line = line.rstrip().split()

  datum = []
  for entry in line:
    entry = int(entry)
    datum.append(entry)
  data.append(datum)

# classify by distances
bin5 = [[] for x in range(101)]
bin20 = [[] for x in range(26)]
bin50 = [[] for x in range(11)]
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

outfilename = pid + "_" + sid + "_" + A + "_" + B + "_analysis.txt" 
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
  print("WARNING:  not all pairs accounted for in total pairs!")


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
