# This program answers the question "For two mutations close enough together
# that they might be in the same read, how often are they actually both
# present on a read spanning both their positions?"  Its input is a filtered
# Strelka output file and an indexed BAM and BAI file for a single sample.

# NOTE:  BAM is zero based.  VCF is one based.  Subtract 1 from your VCF
# mutation position before looking it up in BAM!

import pysam

strelkafile = "/home/mkkuhner/seagate/data/wgs_dna/filtered_strelka/163-23740-14O-37--163-25059N-78R-BLD.snv.strelka.filtered.vcf"

samfile = pysam.AlignmentFile("163-23740-14O-37.final.bam","rb",
  index_filename="163-23740-14O-37.final.bai",require_index=True)

read_cutoff = 1


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

def nooverlap(read1,read2):
  nbases = read1.get_overlap(read2.reference_start,read2.reference_end)
  if nbases == None: # if there is no comparison possible, return True
    return True
  if nbases == 0:
    return True
  return False

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

################################################################################
    
# test mutation pairs in bam file
pairresults = []
anomalypairs = 0    # too many non-ref/alt bases
noreadspairs = 0      # no usable reads


for mut1,mut2 in mutpairs:
  #mut1 = ['21', 16264255, 'C', 'T', 0.15492957746478872]
  #mut2 = ['21', 16264256, 'G', 'A', 0.09859154929577464]

  chr1,pos1,ref1,alt1,vaf1 = mut1
  chr2,pos2,ref2,alt2,vaf2 = mut2


  # correct for zero versus one based
  pos1 -= 1
  pos2 -= 1

# DEBUG
#  if pos1 != 34310880:
#    continue

  assert chr1 == chr2
  myresult = [0,0,0,0,0]

  found = False

  # we don't want to score a read-pair twice if due to overlap both ends
  # could be scored for the same mutation pair, so once we score an end,
  # we put its name in scored_names and refuse to score the other end.
  scored_names = []
  scored_reads = {}

  count = 0
  reads1 = list(samfile.fetch(chr1,pos1,pos1+1))

  for read in reads1:
    count += 1
    if read.mapping_quality < 25:
      # print "poor mapping quality"
      continue     # throw out poor reads
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    found1 = False
    found2 = False
    for mypair in aligned_pairs:
      if mypair[1] == pos1:   # found the first position
        found1 = True
        base1 = read.query_sequence[mypair[0]]
        pair1 = mypair
      if mypair[1] == pos2:   # found the second position
        found2 = True
        base2 = read.query_sequence[mypair[0]]
        pair2 = mypair
    if not (found1 and found2):  
      # print "does not contain mutations"
      continue     # one of the positions is missing

    if read.query_name in scored_names:  
     #  print "rejected due to overlap"
# check that scored_names are in fact overlapping...
      oread = scored_reads[read.query_name]
      if (read.is_read1 and oread.is_read2) or (read.is_read2 and oread.is_read1):
        if nooverlap(read,oread):
          print "found no overlap!"
      continue   # we already scored this read-pair
    scored_names.append(read.query_name)
    scored_reads[read.query_name] = read

    assert (ref1 != alt1)
    assert (ref2 != alt2)

    if base1 == ref1 and base2 == ref2:  myresult[0] += 1
    if base1 == ref1 and base2 == alt2:  myresult[1] += 1
    if base1 == alt1 and base2 == ref2:  myresult[2] += 1
    if base1 == alt1 and base2 == alt2:  myresult[3] += 1
    myresult[4] += 1
    # print "scored"
     
  if myresult[4] == 0:
    noreadspairs += 1
    print chr1,pos1,pos2,myresult
    continue

  if myresult[4] - 5 >= sum(myresult[0:4]):
    # print "Found 5+ non-ref/alt bases:",myresult
    # print "At mutation pair",mut1,mut2
    anomalypairs += 1
    print chr1,pos1,pos2,myresult
    continue

  print chr1,pos1,pos2,myresult
  pairresults.append(myresult)

# print "Number of reads",count

# Note that we do not test for how many non-ref/non-alt bases we saw!

#print pairresults

# summarize the results

# taxonomy:
#   less than 12 total reads:  toosmall
#   just bin0:  wt
#   just bin1 or just bin2:  missed
#   all three bins:  fourgamete
#   bins 1 and 2 but not 3:  trans
#   bin 3 alone:  cis
#   bin 3 with bin 1 or 2 but not both:  nested

toosmall = 0
wt = 0
missed = 0
fourgamete = 0
trans = 0
cis = 0
nested = 0

def hasbin(result,bin):
  if result[bin] > 1:
    return True
  return False

for result in pairresults:

  # toosmall
  if result[4] < 12:
    toosmall += 1
    continue

  # wt
  if not hasbin(result,1) and not hasbin(result,2) and not hasbin(result,3):
    wt += 1
    continue

  # missed
  if hasbin(result,1) and not hasbin(result,2) and not hasbin(result,3):
    missed += 1
    continue
  if not hasbin(result,1) and hasbin(result,2) and not hasbin(result,3):
    missed += 1
    continue

  # fourgamete
  if hasbin(result,1) and hasbin(result,2) and hasbin(result,3):
    fourgamete += 1
    continue

  # trans
  if hasbin(result,1) and hasbin(result,2) and not hasbin(result,3):
    trans += 1
    continue

  # cis
  if not hasbin(result,1) and not hasbin(result,2) and hasbin(result,3):
    cis += 1
    continue

  # nested
  if (hasbin(result,1) or hasbin(result,2)) and hasbin(result,3):
    nested += 1
    continue

  print "found unclassifiable mutpair:",result

print "Total pairs",len(mutpairs)
print "cis",cis,"trans",trans,"nested",nested
print "toosmall",toosmall,"wt",wt,"missed",missed,"fourgamete",fourgamete
print "noreads",noreadspairs,"anomaly",anomalypairs
total = cis + trans + nested + toosmall + wt + missed + fourgamete
total += noreadspairs + anomalypairs
if total != len(mutpairs):
  print "WARNING:  not all pairs accounted for!"
