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
from os import readlink
from os.path import isfile
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

BAF_links = False
if BAF_links:
    BAF_dir = "gamma_template/"
else:
    BAF_dir = "pASCAT_input_combined_everything/"

subdirs = ["diploid", "tetraploid"]
#subdirs = ["diploid"]
subdirdict = {}
for subdir in subdirs:
    subdirdict[subdir] = []
gamma_outdir = "gamma_test_output/"
#gamma_outdir = "gamma_test_output/"
outdir = "gamma_test_output/summaries/"
#gamma_list = ["100", "150", "200", "250", "300", "350", "400", "450", "500", "600"]
#gamma_list = ["50"]
#gamma_list = ["Test"]
gamma_list = ["100", "500", "1000", "3000"]
#gamma_list = ["0", "50", "100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000"]
#gamma_list = []
bafrawdata = {}
patient_samples = {}

showgraphs = False

onlyonepatient = False
onepatient = "575"

twopatientcompare = False
twopatients = [("360", "672")]

compareEverythingToGamma100 = True

compareXiaohongSegments = True

#twopatients = [("521", "252")]
#samples from patient 360:
onepatientsamples = ["360_21620", "360_18370", "360_18354", "360_18358", "360_21608", "360_21606", "360_21604", "360_21602", "360_18356", "360_21600", "360_24412", "360_18350", "360_18444", "360_18446", "360_18352", "360_18366", "360_21618", "360_21598", "360_18362", "360_21596", "360_21614", "360_21594", "360_21616", "360_21592", "360_21610", "360_21612", "360_24409", "360_24810", "360_24813", "360_18372"]
ignore_cnvis = True
mirror_percentages = False

#median for patient 521:  .5154885
#Reasonable cutoffs for 521: 0.48 and 0.563 (about 750k left each)
#A more agressive but balanced cutoff:  0.45 and 0.6  (about 215k left each)
#The best cutoff is no cutoff at all:  0.5 and 0.5
bafErrorBarLow = 0.5
bafErrorBarHigh = 0.5

print("Using BAF filtering for comparison from ", str(bafErrorBarLow), "to", str(bafErrorBarHigh))

bafWtLow = 0.4
bafWtHigh = 0.65

#original cutoff: 0.35 to 0.65, or 0.15 out from 0.5.  The above values account for dye bias.

def readBafNormal(patient):
    bafrawdata = {}
    bafwt = {}
#    allcnvis = []
    if BAF_links:
        bafnormal = readlink(BAF_dir + patient + "_Normal_BAF.txt")
    else:
        bafnormal = BAF_dir + patient + "_Normal_BAF.txt"
    if not(isfile(bafnormal)):
        print("ERROR:  no Normal BAF file found for patient", patient)
        return ({}, {})
    bafnormal = open(bafnormal, "r")
    print("Reading BAF normal data for patient", patient)
    for line in bafnormal:
        lvec = line.split()
        if line.find("Chr") != -1:
            continue
        if ignore_cnvis and line.find("cnvi") != -1:
            continue
        try:
            value = float(lvec[3])
        except:
            continue
#        if (line.find("cnvi") != -1):
#            allcnvis.append(value)
        if (value < bafWtLow or value > bafWtHigh):
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if chr not in bafrawdata:
            bafrawdata[chr] = {}
        if chr not in bafwt:
            bafwt[chr] = {}
        bafrawdata[chr][pos] = {}
        bafwt[chr][pos] = value
    bafnormal.close()
#    print("All normal CNVI bafs:")
#    lsl.createPrintAndSaveHistogram(allcnvis, "", .01)
    return bafrawdata, bafwt

def readBafSamples(baffile, bafrawdata):
    labels = []
    allsamples = []
    if BAF_links:
        baffile = readlink(BAF_dir + patient + "_BAF.txt")
    else:
        baffile = BAF_dir + patient + "_BAF.txt"
    if not(isfile(baffile)):
        print("ERROR:  no BAF file found for patient", patient)
        bafrawdata = {}
        return
    print("Reading BAF sample data for patient", patient)
    baffile = open(baffile, "r")
    allbafs = {}
    allbafs["cnvi"] = []
    allbafs["normal"] = []
    for line in baffile:
        lvec = line.split()
        if line.find("Chr") != -1:
            labels = lvec
            for l in range(2,len(labels)):
                allsamples.append(labels[l].split('"')[1])
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if chr not in bafrawdata:
            continue
        if pos not in bafrawdata[chr]:
            continue
        for p in range(3,len(lvec)):
            sample = labels[p-1].split('"')[1]
            try:
                bafrawdata[chr][pos][sample] = float(lvec[p])
                if line.find("cnvi") != -1:
                    allbafs["cnvi"].append(float(lvec[p]))
                else:
                    bafrawdata[chr][pos][sample] = float(lvec[p])
                    allbafs["normal"].append(float(lvec[p]))
            except:
                continue
    baffile.close()
    if showgraphs:
        print("Sample CNVI bafs for normal CNVIs that were 0.5 in wt:")
        lsl.createPrintAndSaveHistogram(allbafs["cnvi"], "", .01)
#    print("all other bafs:")
#    lsl.createPrintAndSaveHistogram(allbafs["normal"], "", .01)
    return allbafs, allsamples

def compareNormalBafs(bafwt, patient1, patient2):
    scatterx = []
    scattery = []
    for chr in bafwt[patient1]:
        for pos in bafwt[patient1][chr]:
            if pos in bafwt[patient2][chr]:
                scatterx.append(bafwt[patient1][chr][pos])
                scattery.append(bafwt[patient2][chr][pos])
    plt.scatter(scatterx, scattery, marker=".")
    plt.gcf().set_size_inches(5,5)
    plt.show()
    plt.close()
    print("Same thing as a heat map:")
    heatmap, xedges, yedges = numpy.histogram2d(scatterx, scattery, bins=200)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower')
    plt.show()
    plt.close()

def combineTwoBafs(patient1, patient2):
    bafrawdata = {}
    bafwt = {}
    brd = {}
    allwtbafs = []
    for patient in (patient1, patient2):
        brd[patient] = {}
        bafwt[patient] = {}
        if BAF_links:
            bafnormal = readlink(BAF_dir + patient + "_Normal_BAF.txt")
        else:
            bafnormal = BAF_dir + patient + "_Normal_BAF.txt"
        if not(isfile(bafnormal)):
            print("ERROR:  no Normal BAF file found for patient", patient)
            return ({}, {})
        bafnormal = open(bafnormal, "r")
        print("Reading BAF normal data for patient", patient)
        for line in bafnormal:
    #        if (line.find("cnvi") != -1):
    #            continue
            lvec = line.split()
            if line.find("Chr") != -1:
                continue
            try:
                value = float(lvec[3])
            except:
                continue
            if (value < 0.1 or value > 0.9):
                continue
            allwtbafs.append(value)
            if (value < bafWtLow or value > bafWtHigh):
                continue
            chr = lvec[1].split('"')[1]
            pos = int(lvec[2])
            if chr not in brd[patient]:
                brd[patient][chr] = {}
            if chr not in bafwt[patient]:
                bafwt[patient][chr] = {}
                #print("Adding", chr, "to patient", patient)
            brd[patient][chr][pos] = {}
            bafwt[patient][chr][pos] = value
        bafnormal.close()
#        print("Number of 0.5 BAFs for patient", patient, ":")
#        for chr in brd[patient]:
#            print(len(brd[patient][chr]))
    if showgraphs:
        compareNormalBafs(bafwt, patient1, patient2)
        lsl.createPrintAndSaveHistogram(allwtbafs, "", 0.01)

    todelete = []
    for chr in brd[patient1]:
        for pos in brd[patient1][chr]:
            if pos not in brd[patient2][chr]:
                todelete.append((chr, pos))
    for (chr, pos) in todelete:
        del brd[patient1][chr][pos]
    bafrawdata = brd[patient1]
#    print("Number of 0.5 BAFs for both combined:")
#    for chr in bafrawdata:
#        print(len(bafrawdata[chr]))

    allsamples = []
    for patient in [patient1, patient2]:
        labels = []
        if BAF_links:
            baffile = readlink(BAF_dir + patient + "_BAF.txt")
        else:
            baffile = BAF_dir + patient + "_BAF.txt"
        if not(isfile(baffile)):
            print("ERROR:  no BAF file found for patient", patient)
            bafrawdata = {}
            return
        print("Reading BAF sample data for patient", patient)
        baffile = open(baffile, "r")
        allbafs = {}
        allbafs["cnvi"] = []
        allbafs["normal"] = []
        for line in baffile:
            lvec = line.split()
            if line.find("Chr") != -1:
                labels = lvec
                for lv in range(2,len(lvec)):
                    allsamples.append(lvec[lv].split('"')[1])
                continue
            chr = lvec[1].split('"')[1]
            pos = int(lvec[2])
            if chr not in bafrawdata:
                continue
            if pos not in bafrawdata[chr]:
                continue
            for p in range(3,len(lvec)):
                sample = labels[p-1].split('"')[1]
                try:
                    bafrawdata[chr][pos][sample] = float(lvec[p])
                    if line.find("cnvi") != -1:
                        allbafs["cnvi"].append(float(lvec[p]))
                    else:
                        bafrawdata[chr][pos][sample] = float(lvec[p])
                        allbafs["normal"].append(float(lvec[p]))
                except:
                    continue
        baffile.close()
    outfile = open(outdir + "two_patients_input.txt", "w")
    outfile.write("Chr\tpos")
    for sample in allsamples:
        outfile.write("\t" + sample)
    outfile.write("\n")
    for chr in bafrawdata:
        for pos in bafrawdata[chr]:
            outfile.write(chr + "\t" + str(pos))
            for sample in allsamples:
                outfile.write("\t")
                if sample in bafrawdata[chr][pos]:
                    outfile.write(str(bafrawdata[chr][pos][sample]))
                else:
                    outfile.write("NA")
            outfile.write("\n")
    bafwt = bafwt[patient1]
    return bafrawdata, bafwt, allbafs, allwtbafs, allsamples

#cnfile = open(outdir + "copynumber_summaries.txt", "w")
#cnfile.write("patient\tg0\tg0\tg50\tg50\tg100\tg100\tg150\tg150\n")
def readCopynumberFiles(patient, g):
    isegs = {}
    total_segments = {}
    total_seglengths = {}
#    cnfile.write(patient)
    glist = ("0", "50", "100")
    if len(g)>0:
        glist = [g]
    for gamma in glist:
        isegs[gamma] = {}
        total_segments[gamma] = 0
        total_seglengths[gamma] = 0
        root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"
        isegfilename = root_dir + patient + "_copynumber_segments.txt"
        if not(isfile(isegfilename)):
#            cnfile.write("\t\t")
            continue
        isegfile = open(isegfilename, "r")
        for line in isegfile:
            if line.find("Chr") != -1:
                continue
            (chr, start, end, nlogr, nbaf) = line.split()
            nbaf = int(nbaf)
            if nbaf <10:
                continue
            start = int(start)
            end = int(end)
            total_segments[gamma] += 1
            total_seglengths[gamma] += end-start
            if chr not in isegs[gamma]:
                isegs[gamma][chr] = []
            isegs[gamma][chr].append([start, end, {}])
#        cnfile.write("\t" + str(total_segments[gamma]))
#        cnfile.write("\t" + str(total_seglengths[gamma]))
#    cnfile.write("\n")
#    print(total_segments)
#    print(total_seglengths)
    if len(g)>0:
        return isegs[g]
    return isegs

def writeSummaryFileHeader(summary_out, ASC_not_X):
    if ASC_not_X:
        summary_out.write("gamma\t")
    summary_out.write("patient\t")
    summary_out.write("sample\t")
    if ASC_not_X:
        summary_out.write("subdir\t")
        summary_out.write("ploidy\t")
        summary_out.write("purity\t")
    summary_out.write("good matches\t")
    summary_out.write("mismatches\t")
    summary_out.write("missed matches\t")
    summary_out.write("all unbalanced SCA, diploid\t")
    summary_out.write("all unbalanced SCA, tetraploid\t")
    summary_out.write("always-valid SCA, diploid\t")
    summary_out.write("always-valid SCA, tetraploid\t")
    summary_out.write("mixed SCA, diploid\t")
    summary_out.write("mixed SCA, tetraploid\t")
    summary_out.write("always-mismatched SCA, diploid\t")
    summary_out.write("always-mismatched SCA, tetraploid\t")
    summary_out.write("missed SCA, diploid\t")
    summary_out.write("missed SCA, tetraploid\t")
    summary_out.write("\n")

def readPloidyFile(ploidyfile):
    ploidies = {}
    pf = open(ploidyfile)
    for line in pf:
        if line.find('x') != -1:
            continue
        (id, val) = line.split()
        sample = id.split('"')[1].split('_')[1]
        ploidies[sample] = val
    pf.close()
    return ploidies

def readPurityFile(purityfile):
    purities = {}
    pf = open(purityfile)
    for line in pf:
        if line.find('x') != -1:
            continue
        (id, val) = line.split()
        sample = id.split('"')[1].split('_')[1]
        purities[sample] = val
    pf.close()
    return purities

def readSegmentationFile(ascsegfile, allsamples):
    totsca = {}
    totsca["overall"] = set()
    for sample in allsamples:
        totsca[sample] = set()
    osegs = {}
    for line in open(ascsegfile):
        if line.find("Chr") != -1:
            continue
        (__, sample, chr, start, end, nprobes, mbaf, logr, nA, nB) = line.split()
        sample = sample.split('"')[1]
        start = int(start)
        end = int(end)
        if (nA == nB):
            continue
        if sample not in osegs:
            osegs[sample] = {}
        if chr not in osegs[sample]:
            osegs[sample][chr] = []
        osegs[sample][chr].append((start, end))
        totsca[sample].add((chr, start, end))
        totsca["overall"].add((chr, start, end))
    return (osegs, totsca)

def addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca):
    if chr == "":
        return
    if chr == "23":
        return
    if chr == "24":
        return
    if chr=="0":
        #print("There's a chromosome zero segment, weirdly:", full_sample, chr, start, end)
        return
    svec = full_sample.split("-")
    if len(svec) < 4:
        svec = full_sample.split("_")
    patient = svec[0]
    sample = svec[0] + "_" + svec[1]
    if patient not in Xiaohong_segments:
        Xiaohong_segments[patient] = {}
        totsca[patient] = {}
    if sample not in Xiaohong_segments[patient]:
        Xiaohong_segments[patient][sample] = {}
        totsca[patient][sample] = set()
    if chr not in Xiaohong_segments[patient][sample]:
        Xiaohong_segments[patient][sample][chr] = []
    start = int(start)
    end = int(end)
    Xiaohong_segments[patient][sample][chr].append((start, end))
    totsca[patient][sample].add((chr, start, end))
    if "overall" not in totsca[patient]:
        totsca[patient]["overall"] = set()
    totsca[patient]["overall"].add((chr, start, end))
    #print("Adding", chr, start, end)

def readXiaohongWGSLOHFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (__, full_sample, __, chr, start, end) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca)
    xfile.close()

def readXiaohong1MLOHFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (full_sample, chr, start, end, __, __) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca)
    xfile.close()

def readXiaohongCopynumFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    lastseg = ["", 0, 0, 0]
    for line in xfile:
        lvec = line.split()
        if len(lvec) == 7:
            (full_sample, chr, start, end, __, __, code) = line.split()
        elif len(lvec) == 10:
            (__, full_sample, __, chr, start, end, __, __, __, code) = line.split()
        else:
            print("Incorrect line length:", line)
            continue
        start =int(start)
        end = int(end)
        if code == "Double_d" or code=="Balanced_gain" or code == "34" or code=="41":
            #balanced gain or loss: continue
            addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)
            lastseg = ["", 0, 0, 0]
            continue
        if chr == lastseg[0] and code == lastseg[3] and start-lastseg[2] < 5000:
            #print("Combining", chr, str(lastseg[1]), str(lastseg[2]), str(lastseg[3]), "with", start, end, code)
            lastseg[2] = end
        else:
            addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)
            lastseg = [chr, start, end, code]
    xfile.close()
    addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)

def readAllXiaohongSegmentation():
    print("Reading all Xiaohong segmentation.")
    Xiaohong_segments = {}
    totsca = {}
    files = []
    Xdir_WGS = "Xiaohong_WGS_segmentation/"
    Xdir_1M = "CN_Xiaohong_segmentation/"
    for (__, __, f) in walk(Xdir_WGS):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("test") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohongWGSLOHFile(Xdir_WGS + f, Xiaohong_segments, totsca)
        else:
            readXiaohongCopynumFile(Xdir_WGS + f, Xiaohong_segments, totsca)
    files = []
    for (__, __, f) in walk(Xdir_1M):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohong1MLOHFile(Xdir_1M + f, Xiaohong_segments, totsca)
        else:
            readXiaohongCopynumFile(Xdir_1M + f, Xiaohong_segments, totsca)
    return Xiaohong_segments, totsca

def markInputSegmentsWithUnbalancedSamples(isegs, osegs, subdir):
    for sample in osegs:
        for chr in osegs[sample]:
            for oseg in osegs[sample][chr]:
                for iseg in isegs[chr]:
                    if (iseg[0] >= oseg[0] and iseg[1] <= oseg[1]):
                        if subdir not in iseg[2]:
                            iseg[2][subdir] = set()
                        iseg[2][subdir].add(sample)

def collectMatchInfo(isegs, bafrawdata, subdir):
    evidence = {}
    balanced_evidence = {}
    scatter_x = {}
    scatter_y = {}
    key1s = ["in-patient"]
    if twopatientcompare:
        key1s = ("in-patient", "cross-patient")
    key2s = ("balanced-balanced", "balanced-unbalanced", "unbalanced-unbalanced")
    key3s = ("no filter", "medium filter", "high filter")
    for key1 in key1s:
        scatter_x[key1] = {}
        scatter_y[key1] = {}
        for key2 in key2s:
            scatter_x[key1][key2] = {}
            scatter_y[key1][key2] = {}
            for key3 in key3s:
                scatter_x[key1][key2][key3] = []
                scatter_y[key1][key2][key3] = []
    for chr in isegs:
        evidence[chr] = {}
        balanced_evidence[chr] = {}
        for iseg in isegs[chr]:
            isegrange = (iseg[0], iseg[1])
            unbal_samples = set()
            if subdir in iseg[2]:
                unbal_samples = iseg[2][subdir]
            evidence[chr][isegrange] = {}
            balanced_evidence[chr][isegrange] = {}
            for pos in bafrawdata[chr]:
                if pos >= iseg[0] and pos <= iseg[1]:
                    rdsamples = list(bafrawdata[chr][pos].keys())
                    rdsamples.sort()
                    for sample1 in range(0, len(rdsamples)-1):
                        for sample2 in range(sample1+1, len(rdsamples)):
                            s1 = rdsamples[sample1]
                            s2 = rdsamples[sample2]
                            try:
                                val1 = bafrawdata[chr][pos][s1]
                                val2 = bafrawdata[chr][pos][s2]
                                key1 = "cross-patient"
                                if (s1 in onepatientsamples and s2 in onepatientsamples) or (s1 not in onepatientsamples and s2 not in onepatientsamples):
                                    key1 = "in-patient"
                                if s1 not in unbal_samples and s2 not in unbal_samples:
                                    key2 = "balanced-balanced"
                                elif s1 in unbal_samples and s2 in unbal_samples:
                                    key2 = "unbalanced-unbalanced"
                                else:
                                    key2 = "balanced-unbalanced"
                                for key3 in key3s:
                                    if key3 == "no filter":
                                        scatter_x[key1][key2][key3].append(val2)
                                        scatter_y[key1][key2][key3].append(val1)
                                    if key3 == "medium filter" and (val1<0.48 or val1>0.563) and (val2<0.48 or val2>0.563):
                                        scatter_x[key1][key2][key3].append(val2)
                                        scatter_y[key1][key2][key3].append(val1)
                                    if key3 == "high filter" and (val1<0.45 or val1>0.6) and (val2<0.45 or val2>0.6):
                                        scatter_x[key1][key2][key3].append(val2)
                                        scatter_y[key1][key2][key3].append(val1)
                            except:
                                continue
                            if val1 > bafErrorBarLow and val1 < bafErrorBarHigh:
                                continue
                            if val2 > bafErrorBarLow and val2 < bafErrorBarHigh:
                                continue
                            segpair = (s1, s2)
                            if s1 in unbal_samples and s2 in unbal_samples:
                                #both are unbalanced!  Put this in 'evidence'
                                if segpair not in evidence[chr][isegrange]:
                                    evidence[chr][isegrange][segpair] = [0, 0]
                                if (val1 > 0.5 and val2 > 0.5) or (val1 < 0.5 and val2 < 0.5):
                                    evidence[chr][isegrange][segpair][0] += 1
                                    #print("match", val1, val2)
                                else:
                                    evidence[chr][isegrange][segpair][1] += 1
                                    #print("antimatch", val1, val2)
                            else:
                                # One or both are balanced.  Put in 'balanced_evidence'.
                                if segpair not in balanced_evidence[chr][isegrange]:
                                    balanced_evidence[chr][isegrange][segpair] = [0, 0]
                                if (val1 > 0.5 and val2 > 0.5) or (val1 < 0.5 and val2 < 0.5):
                                    balanced_evidence[chr][isegrange][segpair][0] += 1
                                    #print("match", val1, val2)
                                else:
                                    balanced_evidence[chr][isegrange][segpair][1] += 1
                                    #print("antimatch", val1, val2)
    if showgraphs:
        for key1 in key1s:
            for key2 in key2s:
                for key3 in key3s:
                    print("Value comparison for", key1, ",", key2, ",", key3, "with", str(len(scatter_x[key1][key2][key3])), "data points.")

                    if key3 == "no filter":
                        if len(scatter_x[key1][key2][key3]) > 10000000:
                            print("Too many data points: reducing to 10%")
                            sx2 = []
                            sy2 = []
                            for xval in range(0, int(numpy.floor(len(scatter_x[key1][key2][key3])/10))):
                                sx2.append(scatter_x[key1][key2][key3][xval*10])
                                sy2.append(scatter_y[key1][key2][key3][xval*10])
                            scatter_x[key1][key2][key3] = sx2
                            scatter_y[key1][key2][key3] = sy2

                        plt.scatter(scatter_x[key1][key2][key3], scatter_y[key1][key2][key3], marker=".")
                        plt.gcf().set_size_inches(5,5)
                        plt.show()
                        plt.close()

                    heatmap, xedges, yedges = numpy.histogram2d(scatter_x[key1][key2][key3], scatter_y[key1][key2][key3], bins=50)
                    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    plt.clf()
                    plt.imshow(heatmap.T, extent=extent, origin='lower')
                    plt.gcf().set_size_inches(5,5)
                    plt.show()
                    plt.close()
    return (evidence, balanced_evidence)


def processEvidence(evidence, balanced_evidence, osegs, allsamples):
    good_samples = {}
    bad_samples = {}
    missed_samples = {}
    good_sca = {}
    bad_sca = {}
    missed_sca = {}

#    ev_ratios = []
#    balev_ratios = []

    #For information:
    balanced_percs = []
    unbalanced_percs = []
    crosscheck_percs = []
    allbal_percs = []
    balpercs_crosspatient = []
    balpercs_inpatient = []
    unbalpercs_crosspatient = []
    unbalpercs_inpatient = []
    for sample in allsamples:
        good_samples[sample] = 0
        bad_samples[sample] = 0
        missed_samples[sample] = 0
        good_sca[sample] = set()
        bad_sca[sample] = set()
        missed_sca[sample] = set()
    good_samples["overall"] = 0
    bad_samples["overall"] = 0
    missed_samples["overall"] = 0
    good_sca["overall"] = set()
    bad_sca["overall"] = set()
    missed_sca["overall"] = set()
    for chr in evidence:
        for isegrange in evidence[chr]:
            segpercs = []
            minnbaf = 100000000
            chrsegrange = (chr, isegrange[0], isegrange[1])
#            nbases = isegrange[1] - isegrange[0]
            for segpair in evidence[chr][isegrange]:
                [match, antimatch] = evidence[chr][isegrange][segpair]
                tot = match+antimatch
#                ev_ratios.append(math.log(nbases/tot))
                if tot < 10:
                    continue
                minnbaf = min(minnbaf, tot)
                perc = match/tot
                if mirror_percentages and antimatch>match:
                    perc = antimatch/tot
                segpercs.append(perc)
                unbalanced_percs.append(perc)
                if (segpair[0] in onepatientsamples and segpair[1] in onepatientsamples) or (segpair[0] not in onepatientsamples and segpair[1] not in onepatientsamples):
                    unbalpercs_inpatient.append(perc)
                else:
                    unbalpercs_crosspatient.append(perc)
                if perc<0.95 or perc<0.05:
                    #print("bad:", chrsegrange)
                    bad_samples[segpair[0]] += 1
                    bad_samples[segpair[1]] += 1
                    bad_samples["overall"] += 1
                    bad_sca[segpair[0]].add(chrsegrange)
                    bad_sca[segpair[1]].add(chrsegrange)
                    bad_sca["overall"].add(chrsegrange)
                else:
                    #print("good:", chrsegrange)
                    good_samples[segpair[0]] += 1
                    good_samples[segpair[1]] += 1
                    good_samples["overall"] += 1
                    good_sca[segpair[0]].add(chrsegrange)
                    good_sca[segpair[1]].add(chrsegrange)
                    good_sca["overall"].add(chrsegrange)
            for segpair in balanced_evidence[chr][isegrange]:
                [match, antimatch] = balanced_evidence[chr][isegrange][segpair]
                tot = match+antimatch
#                balev_ratios.append(math.log(nbases/tot))
                if tot < 20:
                    continue
                minnbaf = min(minnbaf, tot)
                perc = match/tot
                if mirror_percentages and antimatch>match:
                    perc = antimatch/tot
                #print(match, antimatch, tot, perc)
                balanced_percs.append(perc)
                if (segpair[0] in onepatientsamples and segpair[1] in onepatientsamples) or (segpair[0] not in onepatientsamples and segpair[1] not in onepatientsamples):
                    balpercs_inpatient.append(perc)
                else:
                    balpercs_crosspatient.append(perc)
                unbal_samples = []
                for iseg in isegs[chr]:
                    if iseg[0] == isegrange[0] and subdir in iseg[2]:
                        unbal_samples = iseg[2][subdir]
                if segpair[0] in unbal_samples or segpair[1] in unbal_samples:
                    crosscheck_percs.append(perc)
                else:
                    allbal_percs.append(perc)
                if perc>=0.95 or perc <= 0.05:
                    #print("missed_:", chrsegrange)
                    missed_samples[segpair[0]] += 1
                    missed_samples[segpair[1]] += 1
                    missed_samples["overall"] += 1
                    missed_sca[segpair[0]].add(chrsegrange)
                    missed_sca[segpair[1]].add(chrsegrange)
                    missed_sca["overall"].add(chrsegrange)
    if showgraphs:
        print("All percent matches for all balanced-to-anything segments:")
        lsl.createPrintAndSaveHistogram(balanced_percs, "", .001)
        if twopatientcompare:
            print("Only in-patient percent matches for all balanced-to-anything segments:")
            lsl.createPrintAndSaveHistogram(balpercs_inpatient, "", .001)
            print("Only cross-patient percent matches for all balanced-to-anything segments:")
            lsl.createPrintAndSaveHistogram(balpercs_crosspatient, "", .001)

    #    print("Number of bases/useable SNPs for unbalanced segments:")
    #    lsl.createPrintAndSaveHistogram(ev_ratios, "", 0.1)
    #    print("Mean of unbalanced segment ratios:", numpy.mean(ev_ratios))
    #    print("Median of unbalanced segment ratios:", numpy.median(ev_ratios))
    #    print("Number of bases/useable SNPs for balanced segments:")
    #    lsl.createPrintAndSaveHistogram(balev_ratios, "", 0.1)
    #    print("Mean of balanced segment ratios:", numpy.mean(balev_ratios))
    #    print("Median of balanced segment ratios:", numpy.median(balev_ratios))

        print("Percent matches for all balanced-to-unbalanced checks:")
        lsl.createPrintAndSaveHistogram(crosscheck_percs, "", .001)
        print("Percent matches for all balanced-to-balanced checks:")
        lsl.createPrintAndSaveHistogram(allbal_percs, "", .001)

        print("All percent matches for unbalanced-to-unbalanced segments:")
        lsl.createPrintAndSaveHistogram(unbalanced_percs, "", .001)
        if twopatientcompare:
            print("Only in-patient percent matches for unbalanced-to-unbalanced segments:")
            lsl.createPrintAndSaveHistogram(unbalpercs_inpatient, "", .001)
            print("Only cross-patient percent matches for unbalanced-to-unbalanced segments:")
            lsl.createPrintAndSaveHistogram(unbalpercs_crosspatient, "", .001)

    return (good_samples, bad_samples, missed_samples, good_sca, bad_sca, missed_sca)

def printInfoAbout(missedsca, bafrawdata):
    for seg in missedsca["overall"]:
        (chr, start, end) = seg
        print("segment", seg)
        if (end==4965152):
            for pos in bafrawdata[chr]:
                if pos >= start and pos <= end:
                    print(bafrawdata[chr][pos])

def writeSummaries(summary_out, good_samples, bad_samples, missed_samples, ASC_not_X, totsca, gamma, patient, subdir):
    missed_lengths = {}
    for samp in good_samples:
        missed_lengths[samp] = []
        validatedsca = 0
        mixedsca = 0
        invalidsca = 0
        balancedsca = 0
        sumsca = 0
        if samp in totsca:
            for chrsegpair in totsca[samp]:
                length = chrsegpair[2] - chrsegpair[1]
                #print(length)
                sumsca += length
        if samp in good_sca:
            for chrsegpair in good_sca[samp]:
                length = chrsegpair[2] - chrsegpair[1]
                if chrsegpair in bad_sca[samp]:
                    mixedsca += length
                else:
                    validatedsca += length
        if samp in bad_sca:
            for chrsegpair in bad_sca[samp]:
                length = chrsegpair[2] - chrsegpair[1]
                if chrsegpair not in good_sca[samp]:
                    invalidsca += length
        if samp in missed_sca:
            for chrsegpair in missed_sca[samp]:
                length = chrsegpair[2] - chrsegpair[1]
                balancedsca += length
                missed_lengths[samp].append(length)
        samponly = samp
        if samponly.find("_") != -1:
            samponly = samponly.split("_")[1]
        if ASC_not_X:
            summary_out.write(gamma + "\t")
        summary_out.write(patient + "\t")
        summary_out.write(samponly + "\t")
        if ASC_not_X:
            summary_out.write(subdir + "\t")
            if samp in ploidies:
                summary_out.write(ploidies[samp] + "\t")
                summary_out.write(purities[samp] + "\t")
            else:
                summary_out.write("---\t")
                summary_out.write("---\t")
        summary_out.write(str(good_samples[samp]) + "\t")
        summary_out.write(str(bad_samples[samp]) + "\t")
        summary_out.write(str(missed_samples[samp]) + "\t")
        if subdir == "tetraploid":
            summary_out.write("\t")
        summary_out.write(str(sumsca) + "\t\t")
        summary_out.write(str(validatedsca) + "\t\t")
        summary_out.write(str(mixedsca) + "\t\t")
        summary_out.write(str(invalidsca) + "\t\t")
        summary_out.write(str(balancedsca) + "\t\t")
        summary_out.write("\n")
            #lsl.createPrintAndSaveHistogram(allpercs, outdir + patient + "_" + str(gamma) + ".txt", .001, xdata="Percent match or anti-match")
#            print("Bad samples:")
#            print(bad_samples)
#            print("Good samples:")
#            print(good_samples)
#                for sampIgnore in allpercsminus:
#                    if bad_samples[sampIgnore] < 10 or good_samples[sampIgnore] / bad_samples[sampIgnore] > 0.8:
#                        continue
#                    print("Histogram without any data from sample", sampIgnore, ", which had more bad samples than good:")
#                    lsl.createPrintAndSaveHistogram(allpercsminus[sampIgnore], "", 0.001, xdata="Percent match or anti-match")
#                    summary_out.write(gamma + "\t")
#                    summary_out.write(patient + "\twithout_" + sampIgnore + "\t")
#                    numbad = sum(perc < 0.8 for perc in allpercsminus[sampIgnore])
#                    summary_out.write(str(len(allpercsminus[sampIgnore])-numbad) + "\t")
#                    summary_out.write(str(numbad) + "\t")
#                    summary_out.write("\n")
#    patient_file.close()

if compareXiaohongSegments:
    (Xiaohong_segments, X_totsca) = readAllXiaohongSegmentation()

files = []
for (__, __, f) in walk(BAF_dir):
    files += f
for f in files:
    if f.find("_Normal_BAF.txt") == -1:
        continue
    patient = f.split("_")[0]
    if (onlyonepatient and patient != onepatient):
        continue

#Two-patient comparison code:
#for twopatients in twopatients:
#    patient1, patient2 = twopatients
#    patient = patient1 + patient2
#    bafrawdata, bafwt, allbafs, allwtbafs, allsamples = combineTwoBafs(patient1, patient2)

#Normal one-patient code:
    bafrawdata, bafwt = readBafNormal(patient)
    if (len(bafrawdata)==0):
        continue
    allbafs, allsamples = readBafSamples(patient, bafrawdata)

    if not onlyonepatient and not twopatientcompare and isfile(outdir + patient + "_gamma_summary.txt"):
        print("Skipping patient", patient, ": analysis already exists.")
        continue
    summary_out = open(outdir + patient + "_gamma_summary.txt", "w")
    writeSummaryFileHeader(summary_out, True)
#    patient_file = open(outdir + "gamma_full_" + patient + ".txt", "w")
#    patient_file.write("gamma\tsubdir\tchr\tstart\tend\tminNBafs\tavg\tstdev\tlist\n")

    if compareEverythingToGamma100:
        all_isegs = readCopynumberFiles(patient, "")

    for gamma in gamma_list:
        print("Processing results from a gamma of", gamma)
        root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"
        isegs = {}
        if compareEverythingToGamma100:
            if gamma in all_isegs:
                isegs = deepcopy(all_isegs[gamma])
            elif "100" in all_isegs:
                isegs = deepcopy(all_isegs["100"])
            else:
                isegs = deepcopy(all_isegs["50"])
        else:
            isegs = readCopynumberFiles(patient, gamma)

        for subdir in subdirs:
            ploidyfile = root_dir + subdir + "/" + patient + "_fcn_ascat_ploidy.txt"
            if not(isfile(ploidyfile)):
                print("No ploidy for", patient, "gamma", gamma, "for the", subdir, "constrained run: ASCAT failure for this patient.")
                continue
            ploidies = readPloidyFile(ploidyfile)
            purityfile = root_dir + subdir + "/" + patient + "_fcn_ascat_cont.txt"
            if not(isfile(purityfile)):
                print("No purity for", patient, "gamma", gamma, "for the", subdir, "constrained run: ASCAT failure for this patient.")
                continue
            purities = readPurityFile(purityfile)
            ascsegfile = root_dir + subdir + "/" + patient + "_fcn_ascat_segments.txt"
            if not(isfile(ascsegfile)):
                print("ASCAT seems to have failed for", root_dir, ",",subdir,",",patient,".")
                continue

            (osegs, totsca) = readSegmentationFile(ascsegfile, allsamples)
            markInputSegmentsWithUnbalancedSamples(isegs, osegs, subdir)
            (evidence, balanced_evidence) = collectMatchInfo(isegs, bafrawdata, subdir)

            (good_samples, bad_samples, missed_samples, good_sca, bad_sca, missed_sca) = processEvidence(evidence, balanced_evidence, osegs, allsamples)
            #printInfoAbout(missed_sca, bafrawdata)

            #now count up everything for each sample and output to a file:
            writeSummaries(summary_out, good_samples, bad_samples, missed_samples, True, totsca, gamma, patient, subdir)

    summary_out.close()

    if compareXiaohongSegments:
        print("Running SCA analysis for Xiaohong data for patient", patient)
        xsummary_out = open(outdir + patient + "_xiaohong_summary.txt", "w")
        writeSummaryFileHeader(xsummary_out, False)
        if "100" in all_isegs:
            isegs = deepcopy(all_isegs["100"])
        else:
            isegs = deepcopy(all_isegs["50"])
        subdir = "Xiaohong" #just for 
        markInputSegmentsWithUnbalancedSamples(isegs, Xiaohong_segments[patient], "Xiaohong")
        (evidence, balanced_evidence) = collectMatchInfo(isegs, bafrawdata, subdir)
        (good_samples, bad_samples, missed_samples, good_sca, bad_sca, missed_sca) = processEvidence(evidence, balanced_evidence, Xiaohong_segments[patient], allsamples)
        writeSummaries(xsummary_out, good_samples, bad_samples, missed_samples, False, X_totsca[patient], "", patient, subdir)
        xsummary_out.close()



