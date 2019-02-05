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
somechrsonly = False
somechrs = ['1']

mutation_file = "../snv_plus_indels.twoPlus.20181030.csv"
dipvtet_file = "../calling_evidence_odds.tsv"
cndir = "../noninteger_processed_CNs/"
mutdir = "mutationfiles_bybranch/"

if not(path.isdir(mutdir)):
    mkdir(mutdir)


########################################################################
# functions

def findcnfilename(pid,sid,ploidy,filelist):
  #cnfile = "noninteger_processed_CNs/163_23740_g500_tetraploid_nonint_CNs.txt"
  ploidy = ploidy.lower()
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


def getPatientSampleMap():
    patientSampleMap = {}
#    samplePatientMap = {}
    callfile = open(dipvtet_file, "r")
    for line in callfile:
        if "Patient" in line:
            continue
        lvec = line.rstrip().split()
        (patient, sample) = lvec[0:2]
        ploidy = lvec[-1]
        if ploidy=="Unknown":
            odds = float(lvec[-2])
            if odds > 0.5:
                ploidy = "Diploid"
            else:
                ploidy = "Tetraploid"
        patientSampleMap[sample] = (patient, ploidy)
#        if patient not in samplePatientMap:
#            samplePatientMap[patient] = []
#        samplePatientMap[patient].append(sample)
    return patientSampleMap#, samplePatientMap


def readAllMuts(patientSampleMap):
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
            pos = int(pos)

            (patient, __) = patientSampleMap[sample]
            if patient not in mutations:
                mutations[patient] = {}
            if chr not in mutations[patient]:
                mutations[patient][chr] = {}
            if pos not in mutations[patient][chr]:
                mutations[patient][chr][pos] = {}
            if (ref, alt) not in mutations[patient][chr][pos]:
                mutations[patient][chr][pos][(ref, alt)] = set()
            mutations[patient][chr][pos][(ref, alt)].add(sample)
    return mutations

def convertMapping(mutations):
    altmuts = {}
    for patient in mutations:
        altmuts[patient] = {}
        for chr in mutations[patient]:
            for pos in mutations[patient][chr]:
                for (ref, alt) in mutations[patient][chr][pos]:
                    samples = mutations[patient][chr][pos][(ref, alt)]
                    sl = list(samples)
                    sl.sort()
                    samples = tuple(sl)
                    if samples not in altmuts[patient]:
                        altmuts[patient][samples] = {}
                    if chr not in altmuts[patient][samples]:
                        altmuts[patient][samples][chr] = {}
                    if pos not in altmuts[patient][samples][chr]:
                        altmuts[patient][samples][chr][pos] = (ref, alt)
    return altmuts

def getAllCNCalls(patientSampleMap):
    cncalls = {}
    cnfiles = []
    for root, dirs, files in os.walk(cndir):
      for file in files:
        if file.endswith("_CNs.txt"):
          cnfiles.append(file)
    for sample in patientSampleMap:
        cncalls[sample] = {}
        (patient, ploidy) = patientSampleMap[sample] 
        cnfile = findcnfilename(patient, sample, ploidy, cnfiles)
        if cnfile == None:
            print("FAILED to find a copy number file for patient",patient,"sample",sample)
            print("in directory",cndir)
            exit()

        #################################
        # parse the CN file to find the int calls, and save them.
        for line in open(cndir + cnfile,"r"):
          if "atient" in line:
              # header
             continue
          (mypid,mysid,chr, start, end, rawA, rawB, intA, intB) = line.rstrip().split()
          assert mypid == patient and mysid == sample
          if chr not in cncalls[sample]:
              cncalls[sample][chr] = []
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
          cncalls[sample][chr].append((start, end, A, B))
    return cncalls


########################################################################
# main program


patientSampleMap = getPatientSampleMap()
mutations = readAllMuts(patientSampleMap)
altmuts = convertMapping(mutations)
allCNCalls = getAllCNCalls(patientSampleMap)

for patient in altmuts:
    if somepatientsonly and patient not in somepatients:
        continue
    print("Processing patient", str(patient))
    mutpairs = {}
    for samples in altmuts[patient]:
        if len(samples)==0:
            continue
        mut1 = None
        mut2 = None
        for sample in samples:
#            print(sample)
            for chr in altmuts[patient][samples]:
                if chr == "X" or chr=="23" or chr=="24":
                    continue    # no sex chromosomes please
                if somechrsonly and chr not in somechrs:
                    continue
                prevpos = -500
                prevAB = (-1, -1)
                positions = list(altmuts[patient][samples][chr].keys())
                positions.sort()
                for pos in positions:
                    (ref, alt) = altmuts[patient][samples][chr][pos]
                    ABpair = getCNCallFor(chr, pos, allCNCalls[sample])
                    if sample not in mutpairs:
                        mutpairs[sample] = {}
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
                        if ABpair not in mutpairs[sample]:
                            mutpairs[sample][ABpair] = []
                        mutpairs[sample][ABpair].append([mut1,mut2])
#                        print("Found", sample, str(mut1), str(mut2))
                        assert(mut1[1] < mut2[1])
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
    

    for sample in mutpairs:
        for ABpair in mutpairs[sample]:
            outname = mutdir + patient + "_" + sample + "_" + str(ABpair[0]) + "_" + str(ABpair[1]) + "_mutations.txt"
            outfile = open(outname,"w")
            outfile.write("Chr\tpos1\tref1\talt1\tpos2\tref2\talt2\n")
            mutpairs[sample][ABpair].sort()
            for mut1,mut2 in mutpairs[sample][ABpair]:
              outline = mut1[0]
              for item in mut1[1:]:
                 outline += "\t" + str(item) 
              for item in mut2:
                 outline += "\t" + str(item)
              outline += "\n"
              outfile.write(outline)
            outfile.close()
