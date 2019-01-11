# compute unnormalized mutation signatures from trinucleotide context
# data, 2+ callers, from NYGC

# NB:  This uses the REF trinucleotide context, not the one from the patient's
# normal sample, alas.

# NB:  We assume that the input mutations are sorted by chromosome and position

inputdir = "/home/mkkuhner/seagate/data/wgs_dna/trinucleotide/"
cutoff = 5000000000     # how close do two mutations need to be to count?

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

outfile = open("mutfile_all.tsv","w")
outline = "sample.id\tchr\tpos\tref\talt\n"
outfile.write(outline)

figno = 0
figdir = "mutsigs/"

# process mutation files
for file in infiles:
  pid,sid = file.split("-")[0:2]
  id = pid + "_" + sid
  print("Processing",id)

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
    if prevpos != -10000:
      old_outline = id + "\t" + "chr" + prevchr + "\t" + str(prevpos) + "\t" 
      old_outline += prevref + "\t" + prevalt + "\n"
    outline = id + "\t" + "chr" + chr + "\t" + str(pos) + "\t"
    outline += ref + "\t" + alt + "\n"
    # reverse complement if necessary
    tri,alt = canonical(tri,alt)
    if prevchr == chr and pos - prevpos <= cutoff:   # two close-together mutations!
      if not scored:
        outfile.write(old_outline)
        scores[prevtri+prevalt] += 1
      scores[tri+alt] += 1
      outfile.write(outline)
      scored = True
    else:
      scored = False
    prevchr = chr
    prevpos = pos
    prevtri = tri
    prevref = ref
    prevalt = alt
  allscores[id] = scores


  myscores = allscores[id]
  orderedscores = []
  pointlabels = []
  colors = []
  for b2 in ["C","T"]:
    if b2 == "C":
      others = ["A","G","T"]
    else:
      others = ["A","C","G"]
    for b4 in others:
      for b1 in bases:
        for b3 in bases:
          label = b1+b2+b3+b4
          orderedscores.append(myscores[label])
          pointlabels.append(label)
          if b2 == "C":
            if b4 == "A": colors.append("blue")
            if b4 == "G": colors.append("black")
            if b4 == "T": colors.append("red")
          else:
            if b4 == "A": colors.append("lightgray")
            if b4 == "C": colors.append("lime")
            if b4 == "G": colors.append("pink")

  if False:
    plt.figure(figno)
    plt.title(id)
    barpos = [x for x in range(1,len(orderedscores)+1)]
    plt.bar(barpos,orderedscores,color=colors)
    figno += 1
    plt.savefig(figdir + id+".jpg")
    plt.close(figno)

outfile.close()
