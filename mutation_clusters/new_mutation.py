# This program answers the question "For two mutations close enough together
# that they might be in the same read, how often are they actually both
# present on a read spanning both their positions?"  Its input is a filtered
# Strelka output file and an indexed BAM and BAI file for a single sample.

# NOTE:  BAM is zero based.  VCF is one based.  Subtract 1 from your VCF
# mutation position before looking it up in BAM!

# WARNING:  We wrote this code against a BAM file where paired reads had the
# same name, but this is not standardized:  it will NOT WORK RIGHT if paired
# reads append /1, /2 or _1, _2 or _a, _b or anything like that!

import pysam

strelkafile = "/home/mkkuhner/seagate/data/wgs_dna/filtered_strelka/163-23740-14O-37--163-25059N-78R-BLD.snv.strelka.filtered.vcf"

samfile = pysam.AlignmentFile("163-23740-14O-37.final.bam","rb",
  index_filename="163-23740-14O-37.final.bai",require_index=True)

read_cutoff = 1

# how many non ref/alt bases can we tolerate at a mutation pair
# before throwing it out?
max_unexpected_bases = 5

# lowest acceptable mapping quality for a read
min_quality = 25

# functions

def tier1(entry):
  entry = entry.split(",")
  entry = int(entry[0])
  return entry

def tier1_cutoff(entry, cutoff):
  val = tier1(entry)
  assert val >= 0
  if val >= cutoff:  return val
  else:  return 0

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

#DEBUG
  anomindexes = []
  index = -1

  for tally in pairlist:

#DEBUG
    index += 1

    # noreads
    if tally[4] == 0:
      results["noreads"] += 1
      continue

    # anomaly
    if tally[4] - max_unexpected_bases >= sum(tally[0:4]):
      results["anomaly"] += 1
#DEBUG
      anomindexes.append(index)
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

  return results,anomindexes

########################################################################
# main program

# parse the vcf to find eligible mutation pairs
prevch = 0
prevpos = -500
mutpairs = []
mut1 = None
mut2 = None
for line in open(strelkafile,"r"):
  if line.startswith("#"):  continue     # skip headers
# find a mutation
  chr,pos,id,ref,alt,qual,filter,info,format,normal,tumor = line.rstrip().split()
  pos = int(pos)

  # compute its VAF
  # parse the NORMAL and TUMOR fields to get raw read counts
  trimmed_normal = []
  trimmed_tumor = []
  norm = normal.split(":")[4:]
  for entry in norm:
    trimmed_normal.append(tier1_cutoff(entry,read_cutoff))
  tumor = tumor.split(":")[4:]
  for entry in tumor:
    trimmed_tumor.append(tier1_cutoff(entry,read_cutoff))
  assert sum(trimmed_normal) > 0
  assert sum(trimmed_tumor) > 0

  # compute VAF based on trimmed read-counts
  numer = 0.0
  denom = 0.0
  for norm, tum in zip(trimmed_normal, trimmed_tumor):
    if norm == 0 and tum != 0:
      numer += tum
    denom += tum
  vaf = numer/denom
  if vaf > 0.5:  continue

  # pair the mutations
  muttype = info.split(";")[0]
  if muttype != "NT=ref":  continue      # heterozygous sites are too confusing for now
  if chr == "X":  continue               # no X chromosomes please
  if chr != prevch:     # new chromosome, so no pair
    prevch = chr
    prevpos = pos
    mut1 = [chr,pos,ref,alt,vaf]
    mut2 = None
    continue
  if pos - prevpos > 500:  # too far away, so no pair
    prevpos = pos
    mut1 = [chr,pos,ref,alt,vaf]
    mut2 = None
    continue
  # if we have a previous mutation in hand, this is a pair!
  if mut1 != None:
    mut2 = [chr,pos,ref,alt,vaf]
    mutpairs.append([mut1,mut2])
    mut1 = mut2 = None
    prevpos = -500
  else:     # we already used up mut1 in a different pair, so we skip
    mut1 = [chr,pos,ref,alt,vaf]
    mut2 = None
    prevpos = pos 
    

##########################################################################3

# test mutation pairs in bam file
pairresults = []
pairresults_same = []

for mut1,mut2 in mutpairs:
  # mut1 = ['21', 16264255, 'C', 'T', 0.15492957746478872]
  # mut2 = ['21', 16264256, 'G', 'A', 0.09859154929577464]
  # print("trying a mut-pair")
  chr1,pos1,ref1,alt1,vaf1 = mut1
  chr2,pos2,ref2,alt2,vaf2 = mut2

  # correct for zero versus one based
  pos1 -= 1
  pos2 -= 1

#  #DEBUG
#  if pos1 != 34310880:
#    continue

  assert chr1 == chr2
  myresult = [0,0,0,0,0]
  myresult_sameonly = [0,0,0,0,0]

  # pull reads for mutation position 1 into a list
  reads1 = list(samfile.fetch(chr1,pos1,pos1+1))
  # print("Number of reads",len(reads1))
  reads2 = list(samfile.fetch(chr2,pos2,pos2+1))

  already_scored = []
  for read1 in reads1:
    if read1.mapping_quality < min_quality:  
      #print("poor mapping quality")
      continue    # a bad read
    if read1.query_name in already_scored:  # this read-pair has already been scored
      #print("rejected due to overlap")
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
      #print("does not contain mutations")
      continue

    # both mutation positions 1 and 2 are present:  same-end score
    if found1 and found2:
      score_read(base1,base2,mut1,mut2,myresult_sameonly)
      score_read(base1,base2,mut1,mut2,myresult)
      already_scored.append(read1.query_name)
     # print("scored")
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
      # NO! score_read(base1,base2,mut1,mut2,myresult_sameonly)
      score_read(base1,base2,mut1,mut2,myresult)
      already_scored.append(read1.query_name)
      continue
   
  #myresult_sameonly.append(mut1)  #DEBUG
  #myresult_sameonly.append(mut2)  #DEBUG
  # print(chr1,pos1,pos2,myresult_sameonly)
  pairresults.append(myresult)
  pairresults_same.append(myresult_sameonly)

print("Total pairs",len(mutpairs))
print

#DEBUG added aindex
results,aindex = summarize(pairresults)
print("All reads:")
print("cis",results["cis"],)
print("trans",results["trans"])
print("nested",results["nested"],)
print("toosmall",results["toosmall"],)
print("wt",results["wt"],)
print("missed",results["missed"],)
print("fourgamete",results["fourgamete"])
print("noreads",results["noreads"],)
print("anomaly",results["anomaly"])
print()

total = 0
for categ in results.keys():
  total += results[categ]
if total != len(mutpairs):
  print("WARNING:  not all pairs accounted for in total pairs!")

#DEBUG added aindex
results,aindex = summarize(pairresults_same)
print("Same-end reads:")
print("cis",results["cis"],)
print("trans",results["trans"],)
print("nested",results["nested"])
print("toosmall",results["toosmall"],)
print("wt",results["wt"],)
print("missed",results["missed"],)
print("fourgamete",results["fourgamete"])
print("noreads",results["noreads"],)
print("anomaly",results["anomaly"])

total = 0
for categ in results.keys():
  total += results[categ]
if total != len(mutpairs):
  print("WARNING:  not all pairs accounted for in same-end pairs!")

#DEBUG
#for ind in aindex:
#  print(pairresults_same[ind][5],pairresults_same[ind][6])
