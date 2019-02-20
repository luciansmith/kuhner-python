#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Cproreated on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

from __future__ import division
from os import walk
from os import path
from os import mkdir
from os.path import isfile
import random

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

tag = "_g500_better_ploidy/"

BAF_input = ["BAF_filtered_data_25M_only_40_65/", "BAF_filtered_data_1M_only_40_65/", "BAF_filtered_data_Pilot_40_65/"]
CN_input = "noninteger_processed_CNs/"
dipvtet_file = "calling_evidence_challenge_inc_odds.tsv"
BEAST_output = "BEAST_newinput" + tag


somepatientsonly = True
#somepatients = ["141", "163", "194", "512", "954"]
somepatients = ["891"]

use_nonints = True


if not(path.isdir(BEAST_output + "/")):
    mkdir(BEAST_output + "/")


CNlist = []
for (_, _, f) in walk(CN_input):
    CNlist += f

BAFlist = []
for bafin in BAF_input:
    for (_, _, f) in walk(bafin):
        BAFlist += f


def validateSegments(BAFs_by_sample, validation_output, id, failfile, summaryfile):
    #assumes that any double deletion is already removed.
    output = {}
    samples = list(BAFs_by_sample.keys())
    sampairs = set()
    summatches = 0
    sumantimatches = 0
    sumfails = 0
    samplefails = {}
    for sample in samples:
        samplefails[sample] = [0, 0, 0]
    for n in range(0,len(samples)-1):
        sample1 = samples[n]
        for m in range(n+1, len(samples)):
            sample2 = samples[m]
            matches = []
            for segname in BAFs_by_sample[sample1]:
                if not(segname in BAFs_by_sample[sample2]):
                    continue
                match = 0
                antimatch = 0
                for pos in BAFs_by_sample[sample1][segname]:
                    if not(pos in BAFs_by_sample[sample2][segname]):
                        continue
                    baf1 = BAFs_by_sample[sample1][segname][pos]
                    baf2 = BAFs_by_sample[sample2][segname][pos]
#                    if 0.4 < baf1 and baf1 < 0.6:
#                        continue
#                    if 0.4 < baf2 and baf2 < 0.6:
#                        continue
                    if baf1 < 0.5 and baf2 < 0.5:
                        match += 1
                    elif baf1 > 0.5 and baf2 > 0.5:
                        match += 1
                    elif baf1 < 0.5 and baf2 > 0.5:
                        antimatch += 1
                    elif baf1 > 0.5 and baf2 < 0.5:
                        antimatch += 1
                    #print(baf1, baf2, match, antimatch)
                matches.append((match, antimatch))
                nmax = match + antimatch
                if nmax >= 2:
                    if match/nmax > 0.9:
                        samplefails[sample1][0] += 1
                        samplefails[sample2][0] += 1
                    elif antimatch/nmax > 0.9:
                        samplefails[sample1][1] += 1
                        samplefails[sample2][1] += 1
                    else:
                        samplefails[sample1][2] += 1
                        samplefails[sample2][2] += 1
                if not(segname in output):
                    output[segname] = {}
                output[segname][sample1+"_"+sample2] = (match, antimatch)
                sampairs.add(sample1 + "_" + sample2)
    for sample in samples:
        summaryfile.write(id + "\t" + sample + "\t" + str(samplefails[sample][0]) + "\t" + str(samplefails[sample][1]) + "\t" +  str(samplefails[sample][2])  + "\n")
    outfile = open(validation_output + id + "_segvalidate.txt", "w")
    outfile.write("patient\tchr\tstart\tend\tmax_nBAFs\tmatches\tantimatches\tfails")
    sampairs = sorted(sampairs)
    for sampair in sampairs:
        outfile.write("\t" + sampair + " match\t" + sampair + " anti-match")
    outfile.write("\n")
    for segname in output:
        outfile.write(id + "\t" + segname[0] +"\t" + str(segname[1]) +"\t" + str(segname[2]))
        outline = ""
        nummatches = 0
        numantimatches = 0
        numfails = 0
        maxBAFs = 0
        for sampair in sampairs:
            if sampair in output[segname]:
                match = output[segname][sampair][0]
                antimatch = output[segname][sampair][1]
                outline += "\t" + str(match) + "\t" + str(antimatch)
                numBAFs = match+antimatch
                if maxBAFs<numBAFs:
                    maxBAFs = numBAFs
                if numBAFs<=1:
                    continue
                elif (match/numBAFs) > .9:
                    #print("Match", match, antimatch, numBAFs, match/numBAFs)
                    nummatches += 1
                elif (antimatch/numBAFs) > .9:
                    numantimatches += 1
                    #print("Antimatch", match, antimatch, numBAFs, match/numBAFs)
                else:
                    numfails += 1
                    (sample1, sample2) = sampair.split("_")
                    failfile.write(id + "\t" + segname[0] + "\t" + str(segname[1]) +"\t" + str(segname[2]) + "\t" + sample1)
                    for pos in BAFs_by_sample[sample1][segname]:
                        failfile.write("\t" + str(BAFs_by_sample[sample1][segname][pos]))
                    failfile.write("\n")
                    failfile.write(id + "\t" + segname[0] + "\t" + str(segname[1]) +"\t" + str(segname[2]) + "\t" + sample2)
                    for pos in BAFs_by_sample[sample2][segname]:
                        failfile.write("\t" + str(BAFs_by_sample[sample2][segname][pos]))
                    failfile.write("\n")
                    #print("Fail", match, antimatch, numBAFs, match/numBAFs)
            else:
                outline += "\t--\t--"
        outfile.write("\t" + str(maxBAFs) + "\t" + str(nummatches) + "\t" + str(numantimatches) + "\t" + str(numfails) + outline + "\n")
        summatches += nummatches
        sumantimatches += numantimatches
        sumfails += numfails
    outfile.close()
    summaryfile.write(id + "\tall\t" + str(summatches) + "\t" + str(sumantimatches) + "\t" + str(sumfails) + "\n")


def findCNfilename(pid,sid,ploidy,filelist):
  #cnfile = "noninteger_processed_CNs/163_23740_g500_tetraploid_nonint_CNs.txt"
  ploidy = ploidy.lower()
  for file in filelist:
    entry = file.split("_")
    if entry[0] == pid and entry[1] == sid and entry[3] == ploidy:
      return file
  for file in filelist:
    entry = file.split("_")
    if entry[0] == pid and entry[1] == sid:
      return file
  return None

def getSegmentCalls(patient, sampleList, s2p, CNlist):
    segments = {}
    for sample in sampleList:
        ploidy = s2p[sample][1]
        assert(patient == s2p[sample][0])
        CNfile = findCNfilename(patient, sample, ploidy, CNlist)
        cnf = open(CN_input + CNfile, "r")
        for line in cnf:
            if (line.find("chr") != -1):
                continue
            (__, __, chr, start, end, __, __, intA, intB) = line.rstrip().split()
            if (chr == "23"):
                continue
            if (chr == "24"):
                continue
            chr = int(chr)
            try:
                intA = int(intA)
                intB = int(intB)
                if intB<intA:
                    temp = intA
                    intA = intB
                    intB = temp
            except:
                assert(intA == intB)
                pass #Probably '?'s
            if not(chr in segments):
                segments[chr] = {}
            segpair = (int(start), int(end))
            if not (segpair in segments[chr]):
                segments[chr][segpair] = {}
            segments[chr][segpair][sample] = (intA, intB)
    return segments

def writeHeader(file, samplelist):
    file.write("Chr")
    file.write("\tStart")
    file.write("\tEnd")
    for sample in samplelist:
        file.write("\t" + sample)
    file.write("\n")

def writeLine(file, vector):
    file.write(str(vector[0]))
    for n in range(1,len(vector)):
        file.write("\t")
        file.write(str(vector[n]))
    file.write("\n")

def randomizeFiles(allA, allB):
    if(random.choice([True, False])):
        return (allA, allB)
    else:
        return (allB, allA)

def getMatchVsAntimatch(samp1, samp2, chr, segpair):
    return random.choice(["match", "antimatch"])

def validateMatchgrid(matchgrid):
    set1 = set()
    set2 = set()
    for samp1 in matchgrid:
        if matchgrid[samp1] == "balanced":
            continue
        for samp2 in matchgrid[samp1]:
            if matchgrid[samp1][samp2] == "balanced":
                continue
            if matchgrid[samp1][samp2] == "match":
                if samp1 in set2:
                    set2.add(samp2)
                else:
                    set1.add(samp1)
                    set1.add(samp2)
            elif matchgrid[samp1][samp2] == "antimatch":
                if samp1 in set1:
                    set2.add(samp2)
                elif samp1 in set2:
                    set1.add(samp2)
                elif samp2 in set1:
                    set2.add(samp1)
                elif samp2 in set2:
                    set1.add(samp1)
                else:
                    set1.add(samp1)
                    set2.add(samp2)
    for sample in set1:
        if sample in set2:
            print("Invalid match grid:\n", matchgrid)
            return False
    return True

def finalizeMatches(matchgrid, samples, segments, prevN, prevS):
    Nvec = [-1] * len(samples)
    Svec = [-1] * len(samples)

    # First pass:  any balanced entries, and any entries that match the previous row
    for s in range(len(samples)):
        sample = samples[s]
        if matchgrid[sample] == "balanced":
            Nvec[s] = segments[sample][0]
            Svec[s] = segments[sample][1]
            continue
        if segments[sample][0] == prevN[s] and segments[sample][1] == prevS[s]:
            Nvec[s] = prevN[s]
            Svec[s] = prevS[s]
            continue
        if segments[sample][0] == prevS[s] and segments[sample][1] == prevN[s]:
            Nvec[s] = prevN[s]
            Svec[s] = prevS[s]
            continue

    # Second pass:  anything that matched a pulled-down entry, put that in.
    for s1 in range(len(samples)):
        if Nvec[s1] != -1:
            continue
        samp1 = samples[s1]
        for s2 in range(len(samples)):
            if Nvec[s2] != Svec[s2]:
                samp2 = samples[s2]
                match = matchgrid[samp1][samp2]
                if (Nvec[s2] < Svec[s2] and match == "match") or (Nvec[s2] > Svec[s2] and match == "antimatch"):
                    Nvec[s1] = segments[samp1][0]
                    Svec[s1] = segments[samp1][1]
                else:
                    Nvec[s1] = segments[samp1][1]
                    Svec[s1] = segments[samp1][0]
                break
    
    #Third pass: make a new decision, and then match everything else to it:
    for s1 in range(len(samples)):
        if Nvec[s1] != -1:
            continue
        samp1 = samples[s1]
        Nvec[s1] = segments[samp1][0]
        Svec[s1] = segments[samp1][1]
        for s2 in range(s1+1, len(samples)):
            if Nvec[s2] != -1:
                continue
            samp2 = samples[s2]
            match = matchgrid[samples[s1]][samples[s2]]
            if match == "match":
                Nvec[s2] = segments[samp2][0]
                Svec[s2] = segments[samp2][1]
            else:
                assert(match=="antimatch")
                Nvec[s2] = segments[samp2][1]
                Svec[s2] = segments[samp2][0]
            
        

    return (Nvec, Svec)
    
def getSortedCalls(patient, samples, chr, segpair, segments, prevN, prevS):
    matchgrid = {}
    allBalanced = True
    for s1 in range(len(samples)):
        samp1 = samples[s1]
        (s1A, s1B) = segments[samp1]
        if s1A==s1B:
            matchgrid[samp1] = "balanced"
        else:
            allBalanced = False
            if samp1 not in matchgrid:
                matchgrid[samp1] = {}
            for s2 in range(s1+1, len(samples)):
                samp2 = samples[s2]
                (s2A, s2B) = segments[samp2]
                if s2A!=s2B:
                    matchornot = getMatchVsAntimatch(samp1, samp2, chr, segpair)
                    matchgrid[samp1][samp2] = matchornot
                    if samp2 not in matchgrid:
                        matchgrid[samp2] = {}
                    matchgrid[samp2][samp1] = matchornot
    (Nvec, Svec) = finalizeMatches(matchgrid, samples, segments, prevN, prevS)
    if not(validateMatchgrid(matchgrid)):
        print("Invalid match grid for patient", patient, "samples", str(samples), "at chr", str(chr), str(segpair))
    return (Nvec, Svec, allBalanced)
        
(s2p, p2s) = lsl.getPatientSampleMap(dipvtet_file)
for patient in p2s:
    samples = p2s[patient]
    samples.sort()
    segments = getSegmentCalls(patient, samples, s2p, CNlist)
    allA = open(BEAST_output + patient + "_allA.txt", "w")
    allB = open(BEAST_output + patient + "_allB.txt", "w")
    writeHeader(allA, samples)
    writeHeader(allB, samples)
    chrs = list(segments.keys());
    chrs.sort()
    for chr in chrs:
        shouldSwitch = True
        prevN = [-1]*len(samples)
        prevS = [-1]*len(samples)
        segpairs = list(segments[chr].keys())
        segpairs.sort()
        for segpair in segpairs:
            if (shouldSwitch):
                (Ns, Ss) = randomizeFiles(allA, allB)
            (Nvec, Svec, shouldSwitch) = getSortedCalls(patient, samples, chr, segpair, segments[chr][segpair], prevN, prevS)
            writeLine(Ns, Nvec)
            writeLine(Ss, Svec)
            prevN = Nvec
            prevS = Svec
    

segments = {}
for CNfile in CNlist:
    if use_nonints:
        (patient, sample, gamma, ploidy, __, __) = CNfile.split("_")
        if ploidy != lsl.getBestPloidyFor(patient, sample):
            continue
    else:
        (patient, sample, __) = CNfile.split("_")
    if somepatientsonly and patient not in somepatients:
        continue
    BAFname = patient + "_" + sample + "_BAF.txt"
    BAFname = BAFname.replace('b', '')
    if not(BAFname in BAFlist):
        print("Couldn't find expected BAF file", BAFname, "from CN file", CNfile)
        continue
    if not(patient in segments):
        segments[patient] = {}
    segments[patient][sample] = {}

failfile = open(validation_output + "failures.txt", "w")
summaryfile = open(validation_output + "summary.txt", "w")
failfile.write("patient\tchr\tstart\tend\tsample\tBAFs-.5\n")
summaryfile.write("patient\tsample\tmatches\tantimatches\tfails\n")

for patient in segments:
    if somepatientsonly and patient not in somepatients:
        continue
    BAFs_by_sample = {}
    for sample in segments[patient]:
        #if sample != "18992":
        #    continue
        if not(sample in BAFs_by_sample):
            BAFs_by_sample[sample] = {}
        for chr in segments[patient][sample]:
            for seg in segments[patient][sample][chr]:
                segname = (chr, seg[0], seg[1])
                if not(segname in BAFs_by_sample[sample]):
                    BAFs_by_sample[sample][segname] = {}
        BAFname = patient + "_" + sample + "_BAF.txt"
        BAFname = BAFname.replace('b', '')
        baffilename = BAF_input[0] + BAFname
        if not isfile(baffilename):
            baffilename = BAF_input[1] + BAFname
        if not isfile(baffilename):
            baffilename = BAF_input[2] + BAFname
        baffile = open(baffilename, "r")
        for line in baffile:
            if (line.find("BAF") != -1):
                continue
            (snpid, chr, pos, baf) = line.rstrip().split()
            if (baf=="?"):
                continue
            if (chr == "23"):
                continue
            if (chr == "24"):
                continue
            baf = float(baf)
            pos = int(pos)
            if chr not in segments[patient][sample]:
                continue
            for seg in segments[patient][sample][chr]:
                #These segments come from ASCAT, which is a both-end-inclusive caller (i.e. calls "5-10", "11-24", etc.)
                if pos < seg[0]:
                    continue
                if pos > seg[1]:
                    continue
                segname = (chr, seg[0], seg[1])
                pos = str(pos)
                if pos in BAFs_by_sample[sample][segname]:
                    pos = pos + "_b"
                BAFs_by_sample[sample][segname][pos] = baf
                break
    print("Processing patient", patient)
    validateSegments(BAFs_by_sample, validation_output, patient, failfile, summaryfile)


failfile.close()
summaryfile.close()



