#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:50:09 2017

@author: lpsmith
"""

from __future__ import division
from os import path
from os import mkdir

import lucianSNPLibrary as lsl
import numpy

joint_seg_files = ["joint_seg_v6_diploid.txt", "joint_seg_v6_tetraploid.txt"]
comparison_dir = "joint_seg_comparisons/"

onePatientOnly = False
onePatient = "817"

if not(path.isdir(comparison_dir)):
    mkdir(comparison_dir)

#avlabels, avrev_labels = lsl.getSNPLabelsAveraged(False)
labels, rev_labels = lsl.getSNPLabelsAll(False)


def readsegfile(infilename):
    rfile = open(infilename, "r")
    data = {}
    for line in rfile:
        if line.find("patient") != -1:
            continue
        (patient, sample, chrom, segstart, segend, rawA, rawB, intA, intB) = line.rstrip().split()
        if patient not in data:
            data[patient] = {}
            print "Reading data from patient", patient, "in file", infilename
        if sample not in data[patient]:
            data[patient][sample] = {}
        if chrom not in data[patient][sample]:
            data[patient][sample][chrom] = []
        data[patient][sample][chrom].append((int(segstart), int(segend), int(intA), int(intB)))
    return data

def sortLabels(keys):
    data = {}
    for key in keys:
        (chrom, pos) = key
        if chrom not in data:
            data[chrom] = []
        data[chrom].append(int(pos))
    for chrom in data:
        data[chrom].sort()
    return data

#avlist = sortLabels(avrev_labels.keys())
alllist = sortLabels(rev_labels.keys())

def getSeg(segment):
    (segstart, segend, segA, segB) = segment
    if (segB > segA):
        temp = segB
        segB = segA
        segA = temp
    return (segstart, segend, segA, segB)

def compareChroms(chrom, segs1, segs2, comparison, seg1lengths, seg2lengths):
    comparison["nsegs1"] += len(segs1)
    comparison["nsegs2"] += len(segs2)
    seg1index = 0
    seg2index = 0
    (seg1start, seg1end, seg1A, seg1B) = getSeg(segs1[0])
    (seg2start, seg2end, seg2A, seg2B) = getSeg(segs2[0])
    seg1start = 1
    seg2start = 1
    prevseg1 = (seg1A, seg1B)
    prevseg2 = (seg2A, seg2B)
    prevseg1start = 1
    prevseg2start = 1

    for pos in alllist[chrom]:
        while pos > seg1end:
            seg1index += 1
            if seg1index >= len(segs1):
                (seg1start, seg1end, seg1A, seg1B) = (seg1start, pos+1, seg1A, seg1B)
            else:
                (seg1start, seg1end, seg1A, seg1B) = getSeg(segs1[seg1index])
            thisseg1 = (seg1A, seg1B)
            if (thisseg1 == prevseg1):
                seg1start = prevseg1start
            else:
                seg1lengths.append(seg1end-seg1start)
            prevseg1 = thisseg1
            prevseg1start = seg1start
        while pos > seg2end:
            seg2index += 1
            if seg2index >= len(segs2):
                (seg2start, seg2end, seg2A, seg2B) = (seg2start, pos+1, seg2A, seg2B)
            else:
                (seg2start, seg2end, seg2A, seg2B) = getSeg(segs2[seg2index])
            thisseg2 = (seg2A, seg2B)
            if (thisseg2 == prevseg2):
                seg2start = prevseg2start
            else:
                seg2lengths.append(seg2end-seg2start)
            prevseg2 = thisseg2
            prevseg2start = seg2start
        if pos < seg1start or pos < seg2start:
            continue
        if (seg1A, seg1B) == (seg2A, seg2B):
            comparison["match"] += 1
            if (seg1A, seg1B) not in comparison:
                comparison[(seg1A, seg1B)] = 0
            comparison[(seg1A, seg1B)] += 1
        else:
            comparison["mismatch"] += 1
            key = (seg1A, seg1B, seg2A, seg2B)
            if key not in comparison:
                comparison[key] = 0
            comparison[key] += 1
    seg1lengths.append(seg1end-seg1start)
    seg2lengths.append(seg2end-seg2start)
    #print seg1lengths

for f1 in range(0, len(joint_seg_files)-1):
    for f2 in range(f1+1, len(joint_seg_files)):
        file1 = joint_seg_files[f1]
        file2 = joint_seg_files[f2]
        dset1 = readsegfile(file1)
        dset2 = readsegfile(file2)
        allkeys = set()
        allsamples = []
        allcomparisons = {}
        seg1lengths = []
        seg2lengths = []
        for patient in dset1:
            if onePatientOnly and patient != onePatient:
                continue
            if patient not in dset2:
                print "Skipping patient", patient, ": not present in", file2
                continue
            print "Processing patient", patient
            for sample in dset1[patient]:
                if sample not in dset2[patient]:
                    print "Skipping patient/sample", patient, sample, ": sample not present in", file2
                    continue
                comparison = {}
                comparison["match"] = 0
                comparison["mismatch"] = 0
                comparison["nsegs1"] = 0
                comparison["nsegs2"] = 0
                for chrom in dset1[patient][sample]:
                    if chrom not in dset2[patient][sample]:
                        print "Error: an entire chromosome is missing from", file2, patient, sample
                        foo()
                    compareChroms(chrom, dset1[patient][sample][chrom], dset2[patient][sample][chrom], comparison, seg1lengths, seg2lengths)
                for key in comparison:
                    allkeys.add(key)
                allcomparisons[(patient, sample)] = comparison
        outfilename = "compare_" + file1 + "__" + file2
        if file1 > file2:
            outfilename = "compare_" + file2 + "__" + file1
        outfile = open(comparison_dir + outfilename, "w")
        outfile.write("patient\tsample")
        allkeys = list(allkeys)
        allkeys.sort()
        matchkeys = []
        mismatchkeys = []
        for key in allkeys:
            if len(key)==4:
                mismatchkeys.append(key)
            else:
                matchkeys.append(key)
        for key in matchkeys:
            outfile.write("\t" + str(key))
        for key in mismatchkeys:
            outfile.write("\t" + str(key))
        outfile.write("\n")
        for (patient, sample) in allcomparisons:
            comparison = allcomparisons[(patient, sample)]
            outfile.write(patient + "\t" + sample)
            for key in matchkeys:
                if key in comparison:
                    outfile.write("\t" + str(comparison[key]))
                else:
                    outfile.write("\t0")
            for key in mismatchkeys:
                if key in comparison:
                    outfile.write("\t" + str(comparison[key]))
                else:
                    outfile.write("\t0")
            outfile.write("\n")
        outfile.close()

#        lengthvec = [20, 50, 100, 1000, 10000, 100000, 1000000]
#        seg1binnedlengths = binLengths(lengthvec, seg1lengths)
#        seg2binnedlengths = binLengths(lengthvec, seg2lengths)
        lsl.createPrintAndSaveHistogram(numpy.log10(seg1lengths), comparison_dir + file1 + "_seglengths.txt", 0.01, xdata="Segment lengths", axis=(), show=True)
        lsl.createPrintAndSaveHistogram(numpy.log10(seg2lengths), comparison_dir + file2 + "_seglengths.txt", 0.01, xdata="Segment lengths", axis=(), show=True)

