#! /usr/bin/env python2

# bamscore.py                    Mary Kuhner and Jon Yamato 2018/04/02

# This program receives the output of findmutations.py and combines it with a BAM
# file to answer the question "For two mutations close enough together
# that they might be in the same read, how often are they actually both
# present on a read spanning both their positions?"  Its input is a mutation-pair
# file and an indexed BAM and BAI file for a single sample.  Its output is
# a file with scores (ref/ref, ref/alt, etc) for each mutation pair; we
# save postprocessing for classify.py to get it off the cloud and allow
# easier reruns.

# Currently findmutations.py filters out any position that was not
# called BOTH 1/1 in diploid and 2/2 in tetraploid, if both solutions
# existed.  This means that there will not be mutation files for all
# samples.

# NOTE:  BAM is zero based.  VCF is one based.  Subtract 1 from your VCF
# mutation position before looking it up in BAM!

# NOTE:  For patient 391 sample 23521, use BAM file containing "REI_12051"
# in the name (main pipeline) and not the other one (pilot pipeline).

# WARNING:  We wrote this code against a BAM file where paired reads had the
# same name, but this is not standardized:  it will NOT WORK RIGHT if paired
# reads append /1, /2 or _1, _2 or _a, _b or anything like that!

# This version reads peak numbers from the mutation pairs file and can
# stratify by peaks.


########################################################################
# constants

import pysam

# lowest acceptable mapping quality for a read
min_quality = 25

mutdir = "/home/mkuhner/mutationfiles/"

##########################################################################

# functions

def unpack(line):
  mut1 = [line[0],int(line[1]),line[2],line[3],float(line[4])]
  mut2 = [line[5],int(line[6]),line[7],line[8],float(line[9])]
  peakinfo = [int(line[10]),int(line[11])]
  return [mut1,mut2,peakinfo]

def hasbin(result,bin):
  if result[bin] > 1:
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

########################################################################
# main program

import sys
if len(sys.argv) != 2:
  print("USAGE:  bamscore.py bamfileurl")
  print("This version takes its mutation files from /mutationfiles,")
  print("its BAM/BAI from the S3 cloud, and writes to ." )
  exit()

bamurl = sys.argv[1]
items = bamurl.split("/")
pid,sid,dna,level = items[-1].split("-")

# read in mutation pairs based on VCF and copy-number calls
mutationfile = mutdir + pid + "-" + sid + "-mutations.txt"

mutpairs = []
peaks = []
for line in open(mutationfile,"r"):
  if line.startswith("Peaks"):   # header giving peak data
    peakheader = line
    line = line.rstrip().split()
    assert line[0] == "Peaks:"
    for entry in line[1:]:
      vaf,peakstart,peakend = entry.split(",")
      vaf = float(vaf[1:])   # lose starting paren
      peakstart = float(peakstart)
      peakend = float(peakend[:-1])
      peaks.append([vaf,peakstart,peakend])  # lose ending paren
    numpeaks = len(peaks)
    continue
      
  line = line.rstrip().split()   # body lines giving mutation pairs
  mutpair = unpack(line)
  mutpairs.append(mutpair)

if len(mutpairs) == 0:
  print("No mutation pairs found; bailing out now")
  exit()
print("Assessing",len(mutpairs),"mutation pairs")

##########################################################################3

# test mutation pairs in bam file

baiurl = bamurl[0:-1] + "i"

bamfile = pysam.AlignmentFile(bamurl,"rb", index_filename=baiurl)

pairresults = []
pairresults_same = []
peakresults = [[[] for x in xrange(numpeaks)] for x in xrange(numpeaks)]

for mut1,mut2,peakinfo in mutpairs:
  chr1,pos1,ref1,alt1,vaf1 = mut1
  chr2,pos2,ref2,alt2,vaf2 = mut2
  peak1,peak2 = peakinfo

  # correct for zero versus one based
  pos1 -= 1
  pos2 -= 1

  assert chr1 == chr2
  myresult = [0,0,0,0,0,peak1,peak2]

  # pull reads for mutation position 1 into a list
  reads1 = list(bamfile.fetch(chr1,pos1,pos1+1))
  reads2 = list(bamfile.fetch(chr2,pos2,pos2+1))

  already_scored = []
  for read1 in reads1:
    if read1.mapping_quality < min_quality:  
      continue    # a bad read
    if read1.query_name in already_scored:  # this read-pair has already been scored
      continue
    # find out which, if any, mutation positions are present
    aligned_pairs = read1.get_aligned_pairs(matches_only=True)
    found1 = False
    found2 = False
    for mypair in aligned_pairs:
      if mypair[1] == pos1:
        found1 = True
        base1 = read1.query_sequence[mypair[0]]
      if mypair[1] == pos2:
        found2 = True
        base2 = read1.query_sequence[mypair[0]]

    # mutation position 1 is missing:  skip this read
    if not found1:  
      continue

    # both mutation positions 1 and 2 are present:  same-end score
    if found1 and found2:
      score_read(base1,base2,mut1,mut2,myresult)
      already_scored.append(read1.query_name)
      continue
 
    # only mutation position 1 is present:  find the mate and try an
    # opposite-end score
    if found1 and not found2:
      name1 = read1.query_name
      for read2 in reads2:
        if read2.mapping_quality < min_quality:
          continue
        found2 = False
        name2 = read2.query_name
        if name1 != name2:
          continue
        if (read1.is_read1 == read2.is_read1):   # these are the same read!
          continue
        found2 = True
        break

      if not found2:
        continue

      aligned_pairs2 = read2.get_aligned_pairs(matches_only=True)
      found2 = False
      for mypair in aligned_pairs2:
        if mypair[1] == pos2:
          found2 = True
          base2 = read2.query_sequence[mypair[0]]
          break
    if found1 and found2:    # the mate works, score it
      score_read(base1,base2,mut1,mut2,myresult)
      already_scored.append(read1.query_name)
      continue
   
  pairresults.append(myresult)

############################################################
# write reports


outfilename = pid + "-" + sid + "-results.txt"
outfile = open(outfilename,"w")
outfile.write(peakheader)
for myresult in pairresults:
  outline = ""
  for item in myresult:
    outline += str(item) + " "
  outline += "\n"
  outfile.write(outline)
outfile.close()
