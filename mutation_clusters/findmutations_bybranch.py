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

import os
import csv
from os import path
from os import mkdir

somepatientsonly = False
somepatients = ["1005"]


########################################################################
# functions

def findcnfilename(pid,sid,ploidy,filelist):
  #cnfile = "noninteger_processed_CNs/163_23740_g500_tetraploid_nonint_CNs.txt"
  for file in filelist:
    entry = file.split("_")
    if entry[0] == pid and entry[1] == sid and entry[3] == ploidy:
      return file
  return None

def getCNCallFor(chr, pos, cncalls):
    unknown = (-1, -1)
    if chr not in cncalls:
        return unknown
    prevcall = unknown
    for (start, end, A, B) in cncalls[chr]:
        if pos>= start and pos<= end:
            return (A, B)
        if pos < start and (A, B) == prevcall:
            return (A, B)
        if pos < start:
            return unknown
        prevcall = (A, B)
    return unknown
    

########################################################################
# main program


mutation_file = "../snv_plus_indels.twoPlus.20181030.csv"
dipvtet_file = "../calling_evidence_odds.tsv"
cndir = "../noninteger_processed_CNs/"
mutdir = "mutationfiles/"

if not(path.isdir(mutdir)):
    mkdir(mutdir)

mutations = {}
with open(mutation_file, 'r') as csvfile:
    for lvec in csv.reader(csvfile):
        if "DNANum" in lvec[0]:
            continue
        (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
        if (is_snv=="f"):
            continue
        if (is_2p=="f"):
            continue
#        if ("N" in sample):
#            continue
        if sample not in mutations:
            mutations[sample] = {}
        if chr not in mutations[sample]:
            mutations[sample][chr] = {}
        pos = int(pos)
        mutations[sample][chr][pos] = (ref, alt)


cnfiles = []
for root, dirs, files in os.walk(cndir):
  for file in files:
    if file.endswith("_CNs.txt"):
      cnfiles.append(file)

for line in open(dipvtet_file, "r"):
    if "Patient" in line:
        continue
    lvec = line.rstrip().split();
    ploidy = lvec[-1].lower()
    (pid, sid) = lvec[0:2]
    if ploidy=="unknown" and "N" in sid:
        ploidy = "diploid"
    if ploidy=="unknown":
        print("Unknown ploidy for patient", pid, "sample", sid)
        continue
#    if "N" in sid:
#        print("Skipping gastric sample", pid, sid)
#        continue
    if somepatientsonly and pid not in somepatients:
        continue
    print("Analyzing",pid,sid)

    #If there's no CN file, something has gone wrong.
    cnfile = findcnfilename(pid, sid, ploidy, cnfiles)
    if cnfile == None:
      print("FAILED to find a copy number file for patient",pid,"sample",sid)
      print("in directory",cndir)
      exit()

    #################################
    # parse the CN file to find the int calls, and save them.
    cncalls = {}
    for line in open(cndir + cnfile,"r"):
      if "atient" in line:
          # header
         continue
      (mypid,mysid,chr, start, end, rawA, rawB, intA, intB) = line.rstrip().split()
      assert mypid == pid and mysid == sid
      if chr not in cncalls:
          cncalls[chr] = []
      try:
          A = int(intA)
          B = int(intB)
          start = int(start)
          end = int(end)
      except:
          continue
      if B<A:
          tmp = A
          A = B
          B = tmp
      cncalls[chr].append((start, end, A, B))

    #################################
    # parse the two-plus caller file to find eligible mutation pairs
    
    mut1 = None
    mut2 = None
    mutpairs = {}
    for chr in mutations[sid]:
        if chr == "X" or chr=="23" or chr=="24":
            continue    # no sex chromosomes please
        prevpos = -500
        prevAB = (-1, -1)
        for pos in mutations[sid][chr]:
            (ref, alt) = mutations[sid][chr][pos]
            ABpair = getCNCallFor(chr, pos, cncalls)
            if ABpair == (-1, -1):
                prevAB = ABpair
                mut1 = None
                mut2 = None
                prevpos = -500
                continue
            
            #This is where we'd save VAF information if we wanted it.
    
            # pair the mutations
            if pos - prevpos > 500 or ABpair != prevAB:  # too far away or different CN values.
                prevpos = pos
                mut1 = [chr,pos,ref,alt]
                mut2 = None
                prevAB = ABpair
                continue
            # if we have a previous mutation in hand, this is a pair!
            if mut1 != None:
                mut2 = [chr,pos,ref,alt]
                if ABpair not in mutpairs:
                    mutpairs[ABpair] = []
                mutpairs[ABpair].append([mut1,mut2])
                mut1 = mut2 = None
                prevpos = -500
                prevAB = ABpair
            else:   
                # we already used up mut1 in a different pair, so we skip
                mut1 = [chr,pos,ref,alt]
                mut2 = None
                prevpos = pos 
                prevAB = ABpair
           
           
    ###########################################################################
    # write the output files
    

    for ABpair in mutpairs:
        outname = mutdir + pid + "_" + sid + "_" + str(ABpair[0]) + "_" + str(ABpair[1]) + "_mutations.txt"
        outfile = open(outname,"w")
        outfile.write("Chr\tpos1\tref1\talt1\tpos2\tref2\talt2\n")
        for mut1,mut2 in mutpairs[ABpair]:
          outline = mut1[0]
          for item in mut1[1:]:
             outline += "\t" + str(item) 
          for item in mut2:
             outline += "\t" + str(item)
          outline += "\n"
          outfile.write(outline)
        outfile.close()
