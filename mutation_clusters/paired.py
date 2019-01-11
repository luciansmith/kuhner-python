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
    
print(len(mutpairs))

# test mutation pairs in bam file
pairresults = []
allerrors = []
numbases = 20
count = 0
for mut1,mut2 in mutpairs:
  chr1,pos1,ref1,alt1,vaf1 = mut1
  chr2,pos2,ref2,alt2,vaf2 = mut2

  # correct for zero versus one based
  pos1 -= 1
  pos2 -= 1

  assert chr1 == chr2
  myresult = [0,0,0,0,0]
  errors = [0,0,None,None]

  found = False
  q_seq = []
  qa_seq = []
  namedict = {}
  print("Reads for mutation 1")
  for read in samfile.fetch(chr1,pos1,pos1+1):
    print(read.query_name, read.is_paired, read.is_read1, read.is_read2)
    otherread = samfile.mate(read)
    print(otherread.query_name, otherread.is_paired, otherread.is_read1, otherread.is_read2)

    exit()
    
#    print(read.query_name, read.is_paired, read.is_read1, read.is_read2)
#    if read.query_name not in namedict:
#      namedict[read.query_name] = 0
#    namedict[read.query_name] += 1
#
#  print("Reads for mutation 2")
#  for read in samfile.fetch(chr2,pos2,pos2+1):
#    print(read.query_name, read.is_paired, read.is_read1, read.is_read2)
#    if read.query_name not in namedict:
#      namedict[read.query_name] = 0
#    namedict[read.query_name] += 1
#
#
#  for key in namedict.keys():
#    print(key, namedict[key])
#  exit()

# fetch reads for mutation1 and store
# fetch reads for mutation2 and store
# filter reads
# for each read in store1:
#   check for same-name read in store2 (may be same or mate)
#   if you find one:
#     do indexing for each read (ignore that they may be same!)
#     score cis/trans
