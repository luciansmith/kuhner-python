# construct a trinucleotide input data matrix for use in the R package
# deconstructsigs

# NB:  This uses the REF trinucleotide context, not the one from the patient's
# normal sample, alas.

# NB:  We assume that the input mutations are sorted by chromosome and position

inputdir = "/home/mkkuhner/seagate/data/wgs_dna/trinucleotide/"
cutoff = 500     # how close do two mutations need to be to count?

####################################################
# functions

# this is a quasi-function that reverse compliments a base
rev = {}
rev["A"] = "T"
rev["C"] = "G"
rev["G"] = "C"
rev["T"] = "A"

# we only use trinucleotides where the middle base is C or T.  if
# the trinucleotide is not like that, we reverse complement it.

def canonical(trinuc,alt):
  if trinuc[1] == "C" or trinuc[1] == "T":    # already okay
    return trinuc,alt

  # reverse the order and change its base to its partner
  return rev[trinuc[2]] + rev[trinuc[1]] + rev[trinuc[0]], rev[alt]

def makelinewithtabs(columncontents):
  outline = str(columncontents[0])
  for stuff in columncontents[1:]:
    outline += "\t" + str(stuff) 
  outline += "\n"
  return outline

####################################################
# main program

import os
import matplotlib.pyplot as plt
infiles = []
for root,dirs,files in os.walk(inputdir):
  for file in files:
    if file.endswith("trinuc.txt"):
      infiles.append(file)

allscores = {}      # dictionary by sample
bases = ["A","C","G","T"]

outfile = open("deconstruct_input.tsv","w")
writeheader = True

totscorefile = open("total_clustered.csv","w")

# process mutation files
for file in infiles:
  pid,sid = file.split("-")[0:2]
  id = pid + "_" + sid
  if id.endswith("N"):
    print "\tSkipping a presumed normal sample:",id
    continue
  print "Processing",id

  # initialize storage for trinucs
  scores = {}
  for b1 in bases:
    for b2 in ["C","T"]:
      for b3 in bases:
        if b2 == "C":
          others = ["A","G","T"]
        else:
          others = ["A","C","G"]
        for b4 in others:
          key = b1+b2+b3+b4
          scores[key] = 0

  # score mutations by trinuc
  prevchr = "0"
  prevpos = -10000
  prevtri = "NA"
  prevref = "NA"
  prevalt = "NA"
  scored = False
  totalscore = 0
  for line in open(inputdir + file,"r"):
    line = line.rstrip().split()
    # process header line
    if line[0].startswith("CHROM"):  
      chromind = line.index("CHROM")
      posind = line.index("POS")
      trinucind = line.index("TRINUC")
      refind = line.index("REF")
      altind = line.index("ALT")
      continue

    # process non header line
    chr = line[chromind]
    pos = int(line[posind])
    tri = line[trinucind]
    alt = line[altind]
    ref = line[refind]
    # reverse complement if necessary
    tri,alt = canonical(tri,alt)
    if prevchr == chr and pos - prevpos <= cutoff:   # two close-together mutations!
      if not scored:
        scores[prevtri+prevalt] += 1
        totalscore += 1
      scores[tri+alt] += 1
      totalscore += 1
      scored = True
    else:
      scored = False
    prevchr = chr
    prevpos = pos
    prevtri = tri
    prevref = ref
    prevalt = alt
  allscores[id] = scores
  totline = id + "," + str(totalscore) + "\n"
  totscorefile.write(totline)


# must be done in this order?
# AXA AXC AXG AXT CXA CXC CXG CXT GXA GXC GXG GXT TXA TXC TXG TXT
# X = C->A,C->G,C->T,T->A,T->C,T->G
  myscores = allscores[id]
  normscores = [id,]
  header = ["\t",]
  for b2 in ["C","T"]:
    if b2 == "C":
      others = ["A","G","T"]
    else:
      others = ["A","C","G"]
    for b4 in others:
      for b1 in bases:
        for b3 in bases:
          normscores.append(myscores[b1+b2+b3+b4] / float(totalscore))
          header.append('"' + b1+'['+b2+'>'+b4+']'+b3 + '"')

  if writeheader:
    outline = makelinewithtabs(header)
    outfile.write(outline)
    writeheader = False
    
  outline = makelinewithtabs(normscores)
  outfile.write(outline)
     
outfile.close()
totscorefile.close()
