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

mutation_file = "snv_plus_indels.20180919.csv"
ploidy_file = "calling_evidence_odds.tsv"
outdir = "phylip_input/"
samplefile = "20181019_SampleCodeLucianTrees.txt"

def sampleToCode():
    sampleCodeMap = {}
    s2c = open(outdir + samplefile, "r")
    for line in s2c:
        if "Patient" in line:
            continue
        lvec = line.rstrip().split('\t')
        sample = lvec[1]
        scode = lvec[8]
        sampleCodeMap[sample] = scode
    return sampleCodeMap


if not path.isdir(outdir):
    mkdir(outdir)

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
    poscount = 0
    for snum in range(0,len(samples)):
        sample = samples[snum]
        for chr in mutations[sample]:
            if chr not in allsamples:
                allsamples[chr] = {}
            for pos in mutations[sample][chr]:
                (ref, alt) = mutations[sample][chr][pos]
                if pos not in allsamples[chr]:
                    allsamples[chr][pos] = [ref]*(len(samples)+1)
                    poscount += 1
                allsamples[chr][pos][snum] = alt
    outfile = open(outdir + patient + "_phylip_in.txt", "w")
    outfile.write("  " + str(len(samples)+1) + "  " + str(poscount) + "\n")
    outstrings = {}
    for snum in range(0,len(samples)):
        outstrings[snum] = sampleCodeMap[samples[snum]] + " "
    outstrings[snum+1] = "blood     "
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