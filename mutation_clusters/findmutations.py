# findmutations.py            Mary Kuhner and Jon Yamato 2018/04/02

# This program plus bamscore.py work together to answer the question 
# "For two mutations close enough together
# that they might be in the same read, how often are they actually both
# present on a read spanning both their positions?" This program gets
# eligible pairs of mutations and filters out those in areas of
# copy-number abnormality using Lucian's "2v4 intersection" files,
# which list only SNVs in regions where the diploid solution was 1/1
# AND the tetraploid solution was 2/2.  Note that if only one solution
# was produced 2v4 constrains only on that solution, and that there are
# a few cases where the 2v4 is EMPTY except for its header as no
# such intersection positions existed.  We skip those cases.

# NOTE:  BAM is zero based.  VCF is one based.  Subtract 1 from your VCF
# mutation position before looking it up in BAM!

# WARNING:  We wrote this code against a BAM file where paired reads had the
# same name, but this is not standardized:  it will NOT WORK RIGHT if paired
# reads append /1, /2 or _1, _2 or _a, _b or anything like that!

# This version takes in an additional file giving peak and valley
# locations for Lucian's VAF histograms, 1/1 and 2/2 regions only.
# These are used to add a VAF classification to each mutation pair.

import pysam
import os

read_cutoff = 1


########################################################################
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

def whichpeak(vaf,peaks):
  for i in xrange(numpeaks):
    if vaf >= peaks[i][1] and vaf < peaks[i][2]:  return i
  # fell off the end?
  assert vaf == 1.0
  return numpeaks - 1

def findcnfilename(pid,sid,filelist):
  for file in filelist:
    entry = file.split("_")
    if entry[0] == pid and entry[1] == sid:
      return file
  return None
    

########################################################################
# main program

#strelkafile = "/home/mkkuhner/seagate/data/wgs_dna/filtered_strelka/163-23740-14O-37--163-25059N-78R-BLD.snv.strelka.filtered.vcf"
#cnfile = "/home/mkkuhner/seagate/data/wgs_array/noninteger_processed_CNs/163_23740_g1000_tetraploid_nonint_CNs.txt"

peakfile = "/home/mkkuhner/seagate/data/wgs_dna/bam/2v4_peaks.tsv"

strelkalocation = "/home/mkkuhner/seagate/data/wgs_dna/filtered_strelka/"
strelkafiles = []
for root, dirs, files in os.walk(strelkalocation):
  for file in files:
    if file.endswith(".filtered.vcf"):
      strelkafiles.append(file)

cnlocation = "/home/mkkuhner/seagate/data/wgs_joint/2v4_intersection/"
cnfiles = []
for root, dirs, files in os.walk(cnlocation):
  for file in files:
    if file.endswith("_positions.tsv"):
      cnfiles.append(file)

for strelkafile in strelkafiles:   
   prefix = strelkafile.split("/")[-1]
   prefix = prefix.split("-")
   pid = prefix[0]
   sid = prefix[1]
   prefix = pid + "-" + sid + "-"
   print "Analyzing",pid,sid
 
   # note that if the cnfile analysis couldn't be done, it does NOT
   # fail to generate a file, it generates a header-only file.  Failure
   # to find a file here is therefore a fatal error.
   cnfile = findcnfilename(pid,sid,cnfiles)
   if cnfile == None:
     print "FAILED to find a copy number file for patient",pid,"sample",sid
     print "in directory",cnlocation
     exit()

   #################################
   # read the peaks-and-valleys file
   # this consists of alternate peaks and valleys, starting and ending with
   # a peak--there are implicit valleys at 0 and 1.
   # we will convert to a format of [vaf,start,end] where start and end
   # are the valleys surrounding the peak with that vaf.
   
   peaks = []
   peakstart = 0.0
   peaktop = None
   peakend = None
   found = False
   for line in open(peakfile,"r"):
     if line.startswith("Patient"):  continue    # skip header
     mypid,mysid,vaf,prob,kind,nvafs = line.rstrip().split()
     vaf = float(vaf)
     if mypid != pid:  continue
     if mysid != sid:  continue
     # we have an entry for this patient-sample combo
     found = True
     if kind == "Valley":
       peakend = vaf
       mypeak = [peaktop,peakstart,peakend]
       peaks.append(mypeak)
       peakstart = vaf
       continue
     if kind == "Peak":
       peaktop = vaf
       continue

   if not found:
     print "No peaks-and-valleys info for",pid,sid,"--skipping"
     continue
   
   # catch the last peak
   mypeak = [peaktop,peakstart,1.0]
   peaks.append(mypeak)
   numpeaks = len(peaks)

   #################################
   # parse the CN file to find usable mutation positions (those
   # where diploid gave 1/1 and tetraploid gave 2/2)
   goodpositions = []
   for line in open(cnlocation + cnfile,"r"):
     if line.startswith("Patient"):      # header
       continue
     mypid,mysid,mychr,mypos = line.rstrip().split()
     assert mypid == pid and mysid == sid
     goodpositions.append([mychr,int(mypos)])
   if len(goodpositions) == 0:
     print "Skipping",pid,sid,"as no copy-number information was available"
     continue
   
   #################################
   # parse the vcf to find eligible mutation pairs
   prevch = 0
   prevpos = -500
   mutpairs = []
   
   mut1 = None
   mut2 = None
   for line in open(strelkalocation + strelkafile,"r"):
     if line.startswith("#"):  continue     # skip headers
   # find a mutation
     chr,pos,id,ref,alt,qual,filter,info,format,normal,tumor = line.rstrip().split()
     pos = int(pos)
     if [chr,pos] not in goodpositions:   # is this position in a 1/1:2/2 region?  If not, bail
       continue
   
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
   
     # pair the mutations
     muttype = info.split(";")[0]
     if muttype != "NT=ref":  continue      # heterozygous sites are too confusing for now
     if chr == "X":  continue               # no X chromosomes please
     if chr != prevch:     # new chromosome, so no pair
       prevch = chr
       prevpos = pos
       mut1 = [chr,pos,ref,alt,vaf]
       peak1 = whichpeak(vaf,peaks)
       mut2 = None
       continue
     if pos - prevpos > 500:  # too far away, so no pair
       prevpos = pos
       mut1 = [chr,pos,ref,alt,vaf]
       peak1 = whichpeak(vaf,peaks)
       mut2 = None
       continue
     # if we have a previous mutation in hand, this is a pair!
     if mut1 != None:
       mut2 = [chr,pos,ref,alt,vaf]
       peak2 = whichpeak(vaf,peaks)
       peakinfo = [peak1,peak2]
       mutpairs.append([mut1,mut2,peakinfo])
       mut1 = mut2 = None
       prevpos = -500
     else:     # we already used up mut1 in a different pair, so we skip
       mut1 = [chr,pos,ref,alt,vaf]
       peak1 = whichpeak(vaf,peaks)
       mut2 = None
       prevpos = pos 
       
       
   ###########################################################################
   # write the output files
   
   # make peak header
   header = "Peaks:\t"
   for peak in peaks:
     header += "(" + str(peak[0]) + "," + str(peak[1]) + "," + str(peak[2]) + ")  "
   header += "\n"
   
   outname = prefix + "mutations.txt"

   outfile = open(outname,"w")
   outfile.write(header)
   for mut1,mut2,peakinfo in mutpairs:
     outline = mut1[0]
     for item in mut1[1:]:
       outline += "\t" + str(item) 
     for item in mut2:
       outline += "\t" + str(item)
     for item in peakinfo:
       outline += "\t" + str(item)
     outline += "\n"
     outfile.write(outline)
   outfile.close()
