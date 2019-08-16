#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

from __future__ import division
from os import walk
from os import path
from os import mkdir
from random import shuffle

import lucianSNPLibrary as lsl

tag = "_g500_better_ploidy/"

CN_input = "noninteger_processed_CNs/"
validation_input = "segmentation_validation" + tag

tree_output = "joint_processed_segmentation" + tag

somepatientsonly = False
somepatients = ["1001"]

somechromsonly = False
somechroms = ["9"]

if not(path.isdir(tree_output + "/")):
    mkdir(tree_output + "/")

CNlist = []
for (_, _, f) in walk(CN_input):
    CNlist += f

vallist = []
for (_, _, f) in walk(validation_input):
    vallist += f


def isAllWT(segment):
    for sample in segment:
        (intA, intB) = segment[sample][0:2]
        if (intA != 1 or intB != 1):
            return False
    return True

def isAllEven(segment):
    for sample in segment:
        (intA, intB) = segment[sample][0:2]
        if (intA != intB):
            return False
    return True


def getMatchFrom(currentout, whichmatch, validated_labels, thissample, N, S, intA, intB):
    thisN = N
    thisS = S
    for sample in currentout:
        if sample=="label":
            continue
        if currentout[sample][0][1] == currentout[sample][1][1]:
            continue
        #otherwise, it's also uneven, and we need to find out if we match or antimatch
        nmatch = -1
        nantimatch = -1
        for l in range(len(validated_labels)):
            label = validated_labels[l]
            if label.find(sample) == -1:
                continue
            if label.find(thissample) == -1:
                continue
            num = whichmatch[l]
            if num=="--":
                num=0
            else:
                num = int(num)
            if label.find("anti-match"):
                nantimatch = num
            else:
                nmatch = num
        allmatch = nmatch + nantimatch
        if allmatch > 0:
            if nmatch/allmatch > 0.9:
                thisN = currentout[sample][0][0]
                thisS = currentout[sample][1][0]
            elif nantimatch/allmatch > 0.9:
                thisN = currentout[sample][1][0]
                thisS = currentout[sample][0][0]
            elif nmatch >= nantimatch:
                thisN = currentout[sample][0][0] + "?"
                thisS = currentout[sample][1][0] + "?"
            else:
                thisN = currentout[sample][1][0] + "?"
                thisS = currentout[sample][0][0] + "?"
        elif allmatch == 0:
            thisN = thisN + "?"
            thisS = thisS + "?"
    return(thisN, thisS)


def writeOneSet(chr, oneset, n, samples, output, labels):
    N = "N" + str(n)
    S = "S" + str(n)
    current = len(output)
    for r in range(len(oneset)):
        (segpair, row) = oneset[r]
        output.append({})
        output[current+r]["label"] = str(chr) + "\t" + str(segpair[0]) + "\t" + str(segpair[1])
        #First bring down exactly one unbalanced call from the previous row, if any exist that match exactly.
        mixed_samples = samples.copy()
        shuffle(mixed_samples)
        chosen_sample = ""
        #First bring down exactly one unbalanced call from the previous row, if any exist that match exactly.
        for sample in mixed_samples:
            (intA, intB, __, __, __, whichmatch) = row[sample]
            if (r>0 and intA != intB and intA == output[current+r-1][sample][0][1] and intB == output[current+r-1][sample][1][1]):
                thisN = output[current+r-1][sample][0][0]
                thisS = output[current+r-1][sample][1][0]
                output[current+r][sample] = ((thisN, intA), (thisS, intB))
                chosen_sample = sample
                #Only do this for exactly one sample: others with match or not inside this row
                break
        #Now deal with all the other samples, haplotyping them according to match/antimatches.
        for sample in samples:
            if sample==chosen_sample:
                continue
            (intA, intB, __, __, __, whichmatch) = row[sample]
            if (intA == intB):
                output[current+r][sample] = ((N, intA), (S, intB))
            else:
                (thisN, thisS) = getMatchFrom(output[current+r], whichmatch, labels, sample, N, S, intA, intB)
                output[current+r][sample] = ((thisN, intA), (thisS, intB))
    return n+1



validated_segments = {}
validated_labels = {}
for v in vallist:
    if v.find("_") == -1:
        continue
    (patient, __) = v.split("_")
    if somepatientsonly and patient not in somepatients:
        continue
    vfile = open(validation_input + v, "r")
    print("Reading", v)
    for line in vfile:
        if line.find("chr") != -1:
            validated_labels[patient] = line.split("\t")[8:]
            continue
        (vpatient, chr, start, end, nBAFs, matches, antimatches, fails) = line.split()[0:8]
        if somechromsonly and chr in somechroms:
            continue
        assert(vpatient == patient)
        whichmatch = line.rstrip().split()[8:]
        if not(patient in validated_segments):
            validated_segments[patient] = {}
        if not(chr in validated_segments[patient]):
            validated_segments[patient][chr] = {}
        segpair = (start, end)
        validated_segments[patient][chr][segpair] = (nBAFs, matches, antimatches, fails, whichmatch)

CN_segments = {}
for c in CNlist:
    if "nonint_CNs" not in c:
        continue
    (patient, sample, gamma, ploidy) = c.split("_")[0:4]
    if ploidy != lsl.getBestPloidyFor(patient, sample):
        continue
    if somepatientsonly and patient not in somepatients:
        continue
    cfile = open(CN_input + c, "r")
    print("Reading", c)
    for line in cfile:
        if line.find("chr") != -1:
            continue
        (__, __, chr, start, end, __, __, intA, intB) = line.rstrip().split()
        if not patient in CN_segments:
            CN_segments[patient] = {}
        if not sample in CN_segments[patient]:
            CN_segments[patient][sample] = {}
        if not chr in CN_segments[patient][sample]:
            CN_segments[patient][sample][chr] = {}
        segpair = (start, end)
        CN_segments[patient][sample][chr][segpair] = (intA, intB)



for patient in CN_segments:
    if somepatientsonly and patient not in somepatients:
        continue
    print("Processing data for patient", patient)
    all_segments = {}
    for sample in CN_segments[patient]:
#        ps_out = open(tree_output+patient+"_"+sample+"_processed.txt", "w")
#        ps_out.write("chr\tstart\tend\tintA\tintB\tmatches\tantimatches\tfails\n")
        for chr in CN_segments[patient][sample]:
            for segpair in CN_segments[patient][sample][chr]:
                (intA, intB) = CN_segments[patient][sample][chr][segpair]
                try:
                    (vnBAF_SNPs, matches, antimatches, fails, whichmatch) = validated_segments[patient][chr][segpair]
                except:
                    (vnBAF_SNPs, matches, antimatches, fails, whichmatch) = ("0", "0", "0", "0", [])
                #if (vnBAF_SNPs != "0" and vnBAF_SNPs != nBAF_SNPs):
                #    print("BAF/validation difference in nBAF_SNPs:", patient, sample, chr, segpair, nBAF_SNPs, vnBAF_SNPs)
#                ps_out.write(chr + "\t" + segpair[0] + "\t" + segpair[1] + "\t" + intA + "\t" + intB + "\t" + matches + "\t" + antimatches + "\t" + fails + "\n")
                intchr = int(chr)
                if intchr not in all_segments:
                    all_segments[intchr] = {}
                intsegpair = (int(segpair[0]), int(segpair[1]))
                if intsegpair not in all_segments[intchr]:
                    all_segments[intchr][intsegpair] = {}
                if (intA == "nan" or intA == "NA"):
                    all_segments[intchr][intsegpair][sample] = (intA, intB, int(matches), int(antimatches), int(fails), whichmatch)
                else:
                    if int(intB)<int(intA):
                        all_segments[intchr][intsegpair][sample] = (intB, intA, int(matches), int(antimatches), int(fails), whichmatch)
                    else:
                        all_segments[intchr][intsegpair][sample] = (intA, intB, int(matches), int(antimatches), int(fails), whichmatch)

#        ps_out.close()
    joint_out = open(tree_output + patient + "_characters.txt", "w")
    samples = list(CN_segments[patient].keys())
    samples.sort()
    joint_out.write("chr\tstart\tend")
    for sample in samples:
        joint_out.write("\t" + sample)
    joint_out.write("\n")
    n = 1
    output = []
    for chr in range(1,23):
        collecting = False
        if chr not in all_segments:
            continue
        segpairs = list(all_segments[chr].keys())
        segpairs.sort()
        oneset = []
        for segpair in segpairs:
#            if isAllWT(all_segments[chr][segpair]):
#                print("All wt")
#                if collecting==True:
#                    n = writeOneSet(chr, oneset, n, samples, output, validated_labels[patient])
#                collecting = False
            if isAllEven(all_segments[chr][segpair]):
                #print("All even"))
                if collecting==True:
                    n = writeOneSet(chr, oneset, n, samples, output, validated_labels[patient])
                oneset = []
                oneset.append((segpair, all_segments[chr][segpair]))
                writeOneSet(chr, oneset, 0, samples, output, validated_labels[patient])
                oneset = []
                collecting = False
            else:
                collecting = True
                oneset.append((segpair, all_segments[chr][segpair]))
        if collecting == True:
            #got to the end of a chromosome while collecting data
            n = writeOneSet(chr, oneset, n, samples, output, validated_labels[patient])
        #n = n+1
#   writeJointOutput(joint_out, chr, segpairs, output, samples)
    for r in range(len(output)):
        joint_out.write(output[r]["label"])
        for sample in samples:
            NS = output[r][sample]
            joint_out.write("\t" + NS[0][0] + ":" + str(NS[0][1]) + "::" + NS[1][0] + ":" + str(NS[1][1]))
        joint_out.write("\n")
    joint_out.close()
