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

import lucianSNPLibrary as lsl

tag = "_g250_diploid/"

CN_input = "CN_joint_log2rs" + tag
BAF_input = "BAF_joint_vals" + tag
validation_input = "segmentation_validation" + tag

tree_output = "joint_processed_segmentation" + tag

if not(path.isdir(tree_output + "/")):
    mkdir(tree_output + "/")

CNlist = []
for (_, _, f) in walk(CN_input):
    CNlist += f

BAFlist = []
for (_, _, f) in walk(BAF_input):
    BAFlist += f

vallist = []
for (_, _, f) in walk(validation_input):
    vallist += f

validated_segments = {}
validated_labels = {}
for v in vallist:
    if v.find("_") == -1:
        continue
    (patient, __) = v.split("_")
#    if patient != "862":
#        continue
    vfile = open(validation_input + v, "r")
    print("Reading", v)
    for line in vfile:
        if line.find("chr") != -1:
            validated_labels[patient] = line.split("\t")[8:]
            continue
        (patient, chr, start, end, nBAFs, matches, antimatches, fails) = line.split()[0:8]
        whichmatch = line.split()[8:]
        if not(patient in validated_segments):
            validated_segments[patient] = {}
        if not(chr in validated_segments[patient]):
            validated_segments[patient][chr] = {}
        segpair = (start, end)
        validated_segments[patient][chr][segpair] = (nBAFs, matches, antimatches, fails, whichmatch)

BAF_segments = {}
for b in BAFlist:
    if b.find("_")== -1:
        continue
    (patient, sample, __) = b.split("_")
#    if patient != "862":
#        continue
    bfile = open(BAF_input + b, "r")
    print("Reading", b)
    for line in bfile:
        if line.find("chr") != -1:
            continue
        (chr, start, end, rawA, rawB, intA, intB, avgBAF, nBAF_SNPs) = line.split()
        if not patient in BAF_segments:
            BAF_segments[patient] = {}
        if not sample in BAF_segments[patient]:
            BAF_segments[patient][sample] = {}
        if not chr in BAF_segments[patient][sample]:
            BAF_segments[patient][sample][chr] = {}
        segpair = (start, end)
        BAF_segments[patient][sample][chr][segpair] = (rawA, rawB, intA, intB, avgBAF, nBAF_SNPs)

CN_segments = {}
for c in CNlist:
    if c.find("_")== -1:
        continue
    (patient, sample, __) = c.split("_")
#    if patient != "862":
#        continue
    cfile = open(CN_input + c, "r")
    print("Reading", c)
    for line in cfile:
        if line.find("chr") != -1:
            continue
        (chr, start, end, rawA, rawB, intA, intB, avgCN, nCN_SNPs) = line.split()
        if not patient in CN_segments:
            CN_segments[patient] = {}
        if not sample in CN_segments[patient]:
            CN_segments[patient][sample] = {}
        if not chr in CN_segments[patient][sample]:
            CN_segments[patient][sample][chr] = {}
        segpair = (start, end)
        CN_segments[patient][sample][chr][segpair] = (rawA, rawB, intA, intB, avgCN, nCN_SNPs)



for patient in CN_segments:
#    if patient != "862":
#        continue
    print("Processing data for patient", patient)
    all_segments = {}
    seg_maxBAFs = {}
    seg_nCNSNPs = {}
    for sample in CN_segments[patient]:
        ps_out = open(tree_output+patient+"_"+sample+"_processed.txt", "w")
        ps_out.write("chr\tstart\tend\trawA\trawB\tintA\tintB\tavgCN\tnCN_SNPs\tavgBAF\tnBAF_SNPs\tmatches\tantimatches\tfails\n")
        for chr in CN_segments[patient][sample]:
            for segpair in CN_segments[patient][sample][chr]:
                (rawA, rawB, intA, intB, avgCN, nCN_SNPs) = CN_segments[patient][sample][chr][segpair]
                try:
                    (brawA, brawB, bintA, bintB, avgBAF, nBAF_SNPs) = BAF_segments[patient][sample][chr][segpair]
                except:
                    (brawA, brawB, bintA, bintB, avgBAF, nBAF_SNPs) = (rawA, rawB, intA, intB, "NA", 0)
#                if rawA != brawA:
#                    print("CN/BAF difference in rawA", patient, sample, chr, segpair, rawA, brawA)
#                if rawB != brawB:
#                    print("CN/BAF difference in rawB", patient, sample, chr, segpair, rawB, brawB)
#                if intA != bintA:
#                    print("CN/BAF difference in intA", patient, sample, chr, segpair, intA, bintA)
#                if intB != bintB:
#                    print("CN/BAF difference in intB", patient, sample, chr, segpair, intB, bintB)
                try:
                    (vnBAF_SNPs, matches, antimatches, fails, whichmatch) = validated_segments[patient][chr][segpair]
                except:
                    (vnBAF_SNPs, matches, antimatches, fails, whichmatch) = ("0", "0", "0", "0", [])
                #if (vnBAF_SNPs != "0" and vnBAF_SNPs != nBAF_SNPs):
                #    print("BAF/validation difference in nBAF_SNPs:", patient, sample, chr, segpair, nBAF_SNPs, vnBAF_SNPs)
                ps_out.write(chr + "\t" + segpair[0] + "\t" + segpair[1] + "\t" + rawA + "\t" + rawB + "\t" + intA + "\t" + intB + "\t" + avgCN + "\t" + nCN_SNPs + "\t" + avgBAF + "\t" + str(nBAF_SNPs) + "\t" + matches + "\t" + antimatches + "\t" + fails + "\n")
                intchr = int(chr)
                if intchr not in all_segments:
                    all_segments[intchr] = {}
                    seg_maxBAFs[intchr] = {}
                    seg_nCNSNPs[intchr] = {}
                intsegpair = (int(segpair[0]), int(segpair[1]))
                if intsegpair not in all_segments[intchr]:
                    all_segments[intchr][intsegpair] = {}
                    seg_maxBAFs[intchr][intsegpair] = nBAF_SNPs
                    seg_nCNSNPs[intchr][intsegpair] = nCN_SNPs
                if (intA == "nan"):
                    all_segments[intchr][intsegpair][sample] = (intA, intB, avgCN, int(nCN_SNPs), avgBAF, int(nBAF_SNPs), int(matches), int(antimatches), int(fails), whichmatch)
                else:
                    all_segments[intchr][intsegpair][sample] = (int(intA), int(intB), float(avgCN), int(nCN_SNPs), avgBAF, int(nBAF_SNPs), int(matches), int(antimatches), int(fails), whichmatch)

        ps_out.close()
    joint_out = open(tree_output + patient + "_characters.txt", "w")
    samples = list(CN_segments[patient].keys())
    samples.sort()
    joint_out.write("chr\tstart\tend\tmaxnBAFs\tnCN_SNPs")
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
#            if lsl.isAllWT(all_segments[chr][segpair]):
#                print("All wt")
#                if collecting==True:
#                    n = lsl.writeOneSet(chr, oneset, n, joint_out, samples, output, validated_labels[patient])
#                collecting = False
            if lsl.isAllEven(all_segments[chr][segpair]):
                #print("All even"))
                if collecting==True:
                    n = lsl.writeOneSet(chr, oneset, n, joint_out, samples, output, validated_labels[patient], seg_maxBAFs[chr], seg_nCNSNPs[chr])
                oneset = []
                oneset.append((segpair, all_segments[chr][segpair]))
                lsl.writeOneSet(chr, oneset, 0, joint_out, samples, output, validated_labels[patient], seg_maxBAFs[chr], seg_nCNSNPs[chr])
                oneset = []
                collecting = False
            else:
                collecting = True
                oneset.append((segpair, all_segments[chr][segpair]))
        if collecting == True:
            #got to the end of a chromosome while collecting data
            n = lsl.writeOneSet(chr, oneset, n, joint_out, samples, output, validated_labels[patient], seg_maxBAFs[chr], seg_nCNSNPs[chr])
        #n = n+1
#   writeJointOutput(joint_out, chr, segpairs, output, samples)
    for r in range(len(output)):
        joint_out.write(output[r]["label"])
        for sample in samples:
            NS = output[r][sample]
            joint_out.write("\t" + NS[0][0] + ":" + str(NS[0][1]) + "::" + NS[1][0] + ":" + str(NS[1][1]))
        joint_out.write("\n")
    joint_out.close()
