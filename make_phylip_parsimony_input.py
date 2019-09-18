# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:43:41 2018

@author: Lucian
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt
import csv

import lucianSNPLibrary as lsl

highVAF = False
lowVAF = False
pilot = True

tag = ""
if highVAF:
    tag = "_highVAF"
if lowVAF:
    tag = "_lowVAF"
if not pilot:
    tag = "_noPilot"

CNdir = "noninteger_processed_CNs/"
mutation_file = "snv_plus_indels.twoPlus.20181030.csv"
ploidy_file = "calling_evidence_odds.tsv"
outdir = "phylip_input" + tag + "/"
samplefile = "20181031_SampleCodeLucianTrees_agesDiffRounding_withPilot.txt"

if not path.isdir(outdir):
    mkdir(outdir)


def sampleToCode():
    sampleCodeMap = {}
    s2c = open(samplefile, "r")
    levels = {}
    for line in s2c:
        if "Patient" in line:
            continue
        (patient, sample, __, age, level, __, stype, lcode, __, prog, gej, os) = line.rstrip().split('\t')
        if stype=="T1" or stype=="T2" or stype=="Index":
            level = int(level)
            gej = int(gej)
            dist = gej-level
            
            if patient not in levels:
                levels[patient] = []
            levels[patient].append(dist)
    s2c.close()
    
    for patient in levels:
        levels[patient].sort()
        
    s2c = open(samplefile, "r")
    for line in s2c:
        if "Patient" in line:
            continue
        (patient, sample, __, age, level, __, __, lcode, __, prog, gej, os) = line.rstrip().split('\t')
        if lcode=="G":
            scode = sample + "_" + "G" + age
            sampleCodeMap[sample] = scode
            continue
        level = int(level)
        gej = int(gej)
        dist = gej-level
        lcode = "S" #'stomach'
        if dist - levels[patient][0] > levels[patient][-1] - dist:
            lcode = "T" #'teeth'
        if dist - levels[patient][0] == levels[patient][-1] - dist:
            #print("Equal distances:", dist, levels[patient])
            assert(len(levels[patient]) == 4)
            if dist==levels[patient][2]:
                lcode = "T"
        #print(patient, ":", lcode, str(dist), str(levels[patient]))
        scode = sample + "_" + lcode + age
        
        sampleCodeMap[sample] = scode
    return sampleCodeMap

def includeList(pilot = False, T3 = True):
    useable_samples = {}
    for line in open(samplefile, "r"):
        if "Patient" in line:
            continue
        lvec = line.rstrip().split()
        patient = lvec[0]
        sample = lvec[1]
        stype = lvec[6]
        if pilot and stype=="Pilot":
            useable_samples[sample] = patient
        elif T3 and stype=="LongFollowUp":
            useable_samples[sample] = patient
        elif stype=="T1" or stype=="T2" or stype=="Index":
            useable_samples[sample] = patient
        elif stype=="GastricSample":
            continue
        else:
            print("Unknown stype", stype)
    return useable_samples

def readProgressorOrNot():
    progressorMap = {}
    pfile = open("Patient_status.txt", "r")
    for line in pfile:
        if "Patient" in line:
            continue
        (patient, pstat, sex) = line.rstrip().split()
        progressorMap[patient] = pstat
    return progressorMap

def loadDeletions(useable_samples):
    deletions = {}
    for sample in useable_samples:
        patient = useable_samples[sample]
        ploidy = lsl.getBestPloidyFor(patient, sample)
        filename = CNdir + patient + "_" + sample + "_g500_" + ploidy + "_nonint_CNs.txt"
        if not isfile(filename):
            filename = CNdir + patient + "_" + sample + "_g550_" + ploidy + "_nonint_CNs.txt"
        for line in open(filename, "r"):
            if "patient" in line:
                continue
            lvec = line.rstrip().split()
            if lvec[7] == "0" or lvec[8] == "0":
                (chrom, start, end) = lvec[2:5]
                start = int(start)
                end = int(end)
                if patient not in deletions:
                    deletions[patient] = {}
                if sample not in deletions[patient]:
                    deletions[patient][sample] = {}
                if chrom not in deletions[patient][sample]:
                    deletions[patient][sample][chrom] = []
                deletions[patient][sample][chrom].append((start, end))
    return deletions

def isDeleted(patient, sample, chrom, pos, deletions):
    if patient not in deletions:
        return False
    if sample not in deletions[patient]:
        return False
    if chrom not in deletions[patient][sample]:
        return False
    for (start, end) in deletions[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return True
    return False

useable_samples = includeList(pilot = pilot)
deletions = loadDeletions(useable_samples)

mutations = {}
with open(mutation_file, 'r') as csvfile:
    for lvec in csv.reader(csvfile):
        if "DNANum" in lvec[0]:
            continue
        (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
        if sample not in useable_samples:
            continue
        refcnt = int(lvec[47])
        mutcnt = int(lvec[48])
        VAF = mutcnt/(refcnt+mutcnt)
        #print(VAF)
        if highVAF and VAF < 0.25:
            continue
        if lowVAF and VAF >= 0.25:
            continue
        if (is_snv=="f"):
            continue
        if (is_2p=="f"):
            continue
    #    if ("N" in sample):
    #        continue
        if sample not in mutations:
            mutations[sample] = {}
        if chr not in mutations[sample]:
            mutations[sample][chr] = {}
        pos = int(pos)
        mutations[sample][chr][pos] = (ref, alt)

patients = {}
for line in open(ploidy_file, "r"):
    if "Patient" in line:
        continue
    lvec = line.rstrip().split()
    (patient, sample) = lvec[0:2]
    hutchcall = lvec[-1]
    if patient not in patients:
        patients[patient] = {}
    patients[patient][sample] = hutchcall

batchfile = open(outdir + "runphylip.bat", "w")

sampleCodeMap = sampleToCode()
progressorMap = readProgressorOrNot()

for patient in patients:
    origsamples = list(patients[patient].keys())
    samples = []
    for sample in origsamples:
        if sample in mutations:
            samples.append(sample)
        else:
            print("Unknown sample", sample)
    samples.sort()
    allsamples = {}
    allref = {}
    poscount = 0
    for snum in range(0,len(samples)):
        sample = samples[snum]
        for chr in mutations[sample]:
            if chr not in allsamples:
                allsamples[chr] = {}
                allref[chr] = {}
            for pos in mutations[sample][chr]:
                (ref, alt) = mutations[sample][chr][pos]
                allref[chr][pos] = ref
                if pos not in allsamples[chr]:
                    allsamples[chr][pos] = [ref]*(len(samples)+1)
                    poscount += 1
                allsamples[chr][pos][snum] = alt
    for chr in allsamples:
        for pos in allsamples[chr]:
            for snum in range(len(samples)):
                if isDeleted(patient, sample, chr, pos, deletions) and allsamples[chr][pos][snum] == allref[chr][pos]:
                    allsamples[chr][pos][snum] = "?"
    outfile = open(outdir + patient + "_phylip_in.txt", "w")
    outfile.write("  " + str(len(samples)+1) + "  " + str(poscount) + "\n")
    outstrings = {}
    for snum in range(0,len(samples)):
        outstrings[snum] = sampleCodeMap[samples[snum]] + " "
        if (len(outstrings[snum])) > 10:
            outstrings[snum] = outstrings[snum][0:10]
    
    outstrings[snum+1] = "blood"
    if progressorMap[patient] == "NP":
        outstrings[snum+1] += "_non "
    else:
        outstrings[snum+1] += "_prog"
    stringnum = 10
    for chr in allsamples:
        for pos in allsamples[chr]:
            for snum in range(0,len(samples)+1):
                outstrings[snum] += allsamples[chr][pos][snum]
                if stringnum==80:
                    outstrings[snum] += "\n"
            if stringnum==80:
                stringnum = 0
            else:
                stringnum += 1
    for snum in range(0,len(samples)+1):
        outfile.write(outstrings[snum])
        if not(outstrings[snum][-1] == '\n'):
            outfile.write("\n")
    outfile.close()

    infile= open(outdir + patient + "_infile.txt", "w")
    infile.write(patient + "_phylip_in.txt\n")
    infile.write("I\n")
    infile.write("O\n")
    infile.write(str(len(samples)+1) + "\n")
    infile.write("Y\n")
    infile.close()
    
    batchfile.write("./dnapars < " + patient + "_infile.txt\n")
    batchfile.write("mv outfile " + patient + "_outfile.txt\n")
    batchfile.write("mv outtree " + patient + "_outtree.txt\n")

batchfile.close()