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

import numpy

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

BAF_dir = "gamma_template/"
subdirs = ["diploid", "tetraploid"]
subdirdict = {}
for subdir in subdirs:
    subdirdict[subdir] = []
gamma_outdir = "gamma_test_output/"
#gamma_outdir = "gamma_test_output/"
outdir = "gamma_test_output/summaries/"
#gamma_list = ["100", "150", "200", "250", "300", "350", "400", "450", "500", "600"]
#gamma_list = ["0"]
gamma_list = ["0", "50", "100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000", "2500", "3000"]
#gamma_list = ["800"]
bafrawdata = {}
patient_samples = {}

onlyonepatient = True
onepatient = "521"

files = []
for (__, __, f) in walk(BAF_dir):
    files += f
for f in files:
    if f.find("_Normal_BAF.txt") == -1:
        continue
    patient = f.split("_")[0]
    if (onlyonepatient and patient != onepatient):
        continue
    bafrawdata = {}
    bafnormal = readlink(BAF_dir + patient + "_Normal_BAF.txt")
    if not(isfile(bafnormal)):
        print "ERROR:  no Normal BAF file found for patient", patient
        continue
    bafnormal = open(bafnormal, "r")
    baffile = readlink(BAF_dir + patient + "_BAF.txt")
    if not(isfile(baffile)):
        print "ERROR:  no BAF file found for patient", patient
        continue
    labels = []
    print "Reading BAF normal data for patient", patient
    for line in bafnormal:
        lvec = line.split()
        if line.find("Chr") != -1:
            continue
        try:
            value = float(lvec[3])
        except:
            continue
        if (value < .35 or value > .65):
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if chr not in bafrawdata:
            bafrawdata[chr] = {}
        bafrawdata[chr][pos] = {}
    bafnormal.close()

    print "Reading BAF sample data for patient", patient
    baffile = open(baffile, "r")
    for line in baffile:
        lvec = line.split()
        if line.find("Chr") != -1:
            labels = lvec
            continue
        chr = lvec[1].split('"')[1]
        pos = int(lvec[2])
        if chr not in bafrawdata:
            continue
        if pos not in bafrawdata[chr]:
            continue
        for p in range(3,len(lvec)):
            sample = labels[p-1].split('"')[1].split('_')[1]
            try:
                bafrawdata[chr][pos][sample] = float(lvec[p])
            except:
                continue
    baffile.close()

    summary_out = open(outdir + patient + "_gamma_summary.txt", "w")
    summary_out.write("gamma\t")
    summary_out.write("patient\t")
    summary_out.write("sample\t")
    summary_out.write("subdir\t")
    summary_out.write("ploidy\t")
    summary_out.write("purity\t")
    summary_out.write("good matches\t")
    summary_out.write("mismatches\t")
    summary_out.write("all unbalanced SCA, diploid\t")
    summary_out.write("all unbalanced SCA, tetraploid\t")
    summary_out.write("always-valid SCA, diploid\t")
    summary_out.write("always-valid SCA, tetraploid\t")
    summary_out.write("mixed SCA, diploid\t")
    summary_out.write("mixed SCA, tetraploid\t")
    summary_out.write("always-mismatched SCA, diploid\t")
    summary_out.write("always-mismatched SCA, tetraploid\t")
    summary_out.write("\n")
#    patient_file = open(outdir + "gamma_full_" + patient + ".txt", "w")
#    patient_file.write("gamma\tsubdir\tchr\tstart\tend\tminNBafs\tavg\tstdev\tlist\n")

    for gamma in gamma_list:
        print "Processing results from a gamma of", gamma
        root_dir = gamma_outdir + "pASCAT_input_g" + gamma + "/"

        isegs = {}
        isegfilename = root_dir + patient + "_copynumber_segments.txt"
        if not(isfile(isegfilename)):
            print "copynumber not run for patient", patient, "gamma", gamma
            continue
        isegfile = open(isegfilename, "r")
        for line in isegfile:
            if line.find("Chr") != -1:
                continue
            (chr, start, end, nlogr, nbaf) = line.split()
            if nbaf <10:
                continue
            if chr not in isegs:
                isegs[chr] = []
            isegs[chr].append([int(start), int(end), {}, {}])
        for subdir in subdirs:
            totsca = {}
            totsca["overall"] = set()
            osegs = {}
            balsegs = {}
            ploidies = {}
            ploidyfile = root_dir + subdir + "/" + patient + "_fcn_ascat_ploidy.txt"
            if not(isfile(ploidyfile)):
                print "No ploidy for", patient, ", gamma", gamma, ": ASCAT failure for this patient."
                continue
            pf = open(ploidyfile)
            for line in pf:
                if line.find('x') != -1:
                    continue
                (id, val) = line.split()
                sample = id.split('"')[1].split('_')[1]
                ploidies[sample] = val
            pf.close()
            purities = {}
            purityfile = root_dir + subdir + "/" + patient + "_fcn_ascat_cont.txt"
            if not(isfile(purityfile)):
                print "No purity for", patient, ", gamma", gamma, ": ASCAT failure for this patient."
                continue
            pf = open(purityfile)
            for line in pf:
                if line.find('x') != -1:
                    continue
                (id, val) = line.split()
                sample = id.split('"')[1].split('_')[1]
                purities[sample] = val
            pf.close()
            ascsegfile = root_dir + subdir + "/" + patient + "_fcn_ascat_segments.txt"
            if not(isfile(ascsegfile)):
                print "ASCAT seems to have failed for", root_dir, ",",subdir,",",patient,"."
                continue
            for line in open(ascsegfile):
                if line.find("Chr") != -1:
                    continue
                (__, sample, chr, start, end, nprobes, mbaf, logr, nA, nB) = line.split()
                sample = sample.split('"')[1].split("_")[1]
                if sample not in osegs:
                    osegs[sample] = {}
                    totsca[sample] = set()
                start = int(start)
                end = int(end)
                if (nA == nB):
                    if chr not in balsegs[sample]:
                        balsegs[sample][chr]= []
                    balsegs[sample][chr].append((start, end))
                else:
                    if chr not in osegs[sample]:
                        osegs[sample][chr] = []
                    osegs[sample][chr].append((start, end))
                    totsca[sample].add((chr, start, end))
                    totsca["overall"].add((chr, start, end))
            for sample in osegs:
                for chr in osegs[sample]:
                    for oseg in osegs[sample][chr]:
#                        print oseg
                        for iseg in isegs[chr]:
                            if (iseg[0] >= oseg[0] and iseg[1] <= oseg[1]):
                                if subdir not in iseg[2]:
                                    iseg[2][subdir] = []
                                if sample not in iseg[2][subdir]:
                                    iseg[2][subdir].append(sample)
            for sample in balsegs:
                for chr in balsegs[sample]:
                    for balseg in balsegs[sample][chr]:
#                        print balseg
                        for iseg in isegs[chr]:
                            if (iseg[0] >= balseg[0] and iseg[1] <= balseg[1]):
                                if subdir not in iseg[3]:
                                    iseg[3][subdir] = []
                                if sample not in iseg[3][subdir]:
                                    iseg[3][subdir].append(sample)
#                                    print "Just added sample", sample, "to", iseg
#                                    print "balseg is", balseg
            evidence = {}
            for chr in isegs:
                evidence[chr] = {}
                for iseg in isegs[chr]:
                    isegrange = (iseg[0], iseg[1])
                    if subdir not in iseg[2]:
                        continue
                    samples = iseg[2][subdir]
                    samples.sort()
                    evidence[chr][isegrange] = {}
                    for pos in bafrawdata[chr]:
                        if pos >= iseg[0] and pos <= iseg[1]:
                            for sample1 in range(0, len(samples)-1):
                                for sample2 in range(sample1+1, len(samples)):
                                    s1 = samples[sample1]
                                    s2 = samples[sample2]
                                    try:
                                        val1 = bafrawdata[chr][pos][s1]
                                        val2 = bafrawdata[chr][pos][s2]
                                    except:
                                        continue
                                    if val1 > 0.48 and val1 < 0.52:
                                        continue
                                    if val2 > 0.48 and val2 < 0.52:
                                        continue
                                    #This section is redundant with the above bit that appends the sample to iseg[2]:
#                                    found = False
#                                    if s1 not in osegs or s2 not in osegs:
#                                        foo()
#                                        continue
#                                    if chr not in osegs[s1] or chr not in osegs[s2]:
#                                       bar()
#                                       continue
#                                    for s1seg in osegs[s1][chr]:
#                                        if pos >= s1seg[0] and pos <= s1seg[1]:
#                                            found = True
#                                    if found==False:
#                                        foo()
#                                        continue
#                                    found = False
#                                    for s2seg in osegs[s2][chr]:
#                                        if pos >= s2seg[0] and pos <= s2seg[1]:
#                                            found = True
#                                    if found==False:
#                                        continue
                                    segpair = (samples[sample1], samples[sample2])
                                    if segpair not in evidence[chr][isegrange]:
                                        evidence[chr][isegrange][segpair] = [0, 0]
                                    if (val1 > 0.5 and val2 > 0.5) or (val1 < 0.5 and val2 < 0.5):
                                        evidence[chr][isegrange][segpair][0] += 1
                                        #print "match", val1, val2
                                    else:
                                        evidence[chr][isegrange][segpair][1] += 1
                                        #print "antimatch", val1, val2
            allpercs = []
            allpercsminus = {}
            goodsamples = {}
            badsamples = {}
            goodsca = {}
            badsca = {}
            for sample in osegs:
                allpercsminus[sample] = []
                goodsamples[sample] = 0
                badsamples[sample] = 0
                goodsca[sample] = set()
                badsca[sample] = set()
            goodsamples["overall"] = 0
            badsamples["overall"] = 0
            goodsca["overall"] = set()
            badsca["overall"] = set()
            for chr in evidence:
                for isegrange in evidence[chr]:
                    segpercs = []
#                    patient_file.write(gamma + "\t")
#                    patient_file.write(subdir + "\t")
#                    patient_file.write(chr + "\t")
#                    patient_file.write(str(isegrange[0]) + "\t")
#                    patient_file.write(str(isegrange[1]) + "\t")
                    minnbaf = 100000000
                    sca = isegrange[1] - isegrange[0]
                    chrsegrange = (chr, isegrange[0], isegrange[1])
                    for segpair in evidence[chr][isegrange]:
                        [match, antimatch] = evidence[chr][isegrange][segpair]
                        tot = match+antimatch
                        if tot < 10:
                            continue
                        minnbaf = min(minnbaf, tot)
                        perc = match/tot
                        if (antimatch>match):
                            perc = antimatch/tot
                        allpercs.append(perc)
                        segpercs.append(perc)
                        for sampIgnore in osegs:
                            if sampIgnore not in segpair:
                                allpercsminus[sampIgnore].append(perc)
                        if perc<0.8:
                            #print "bad:", chrsegrange
                            badsamples[segpair[0]] += 1
                            badsamples[segpair[1]] += 1
                            badsamples["overall"] += 1
                            badsca[segpair[0]].add(chrsegrange)
                            badsca[segpair[1]].add(chrsegrange)
                            badsca["overall"].add(chrsegrange)
                        else:
                            #print "good:", chrsegrange
                            goodsamples[segpair[0]] += 1
                            goodsamples[segpair[1]] += 1
                            goodsamples["overall"] += 1
                            goodsca[segpair[0]].add(chrsegrange)
                            goodsca[segpair[1]].add(chrsegrange)
                            goodsca["overall"].add(chrsegrange)
#                    if len(segpercs) < 1:
#                        patient_file.write("[None]\t[None]")
#                    else:
#                        patient_file.write(str(minnbaf) + "\t")
#                        patient_file.write(str(numpy.average(segpercs)) + "\t")
#                        patient_file.write(str(numpy.std(segpercs)))
#                        for perc in segpercs:
#                            patient_file.write("\t" + str(perc))
#                    patient_file.write("\n")
            for samp in goodsamples:
                validatedsca = 0
                mixedsca = 0
                invalidsca = 0
                sumsca = 0
                for chrsegpair in totsca[samp]:
                    length = chrsegpair[2] - chrsegpair[1]
                    #print length
                    sumsca += length
                for chrsegpair in goodsca[samp]:
                    length = chrsegpair[2] - chrsegpair[1]
                    if chrsegpair in badsca[samp]:
                        mixedsca += length
                    else:
                        validatedsca += length
                for chrsegpair in badsca[samp]:
                    length = chrsegpair[2] - chrsegpair[1]
                    if chrsegpair not in goodsca[samp]:
                        invalidsca += length
                summary_out.write(gamma + "\t")
                summary_out.write(patient + "\t")
                summary_out.write(samp + "\t")
                summary_out.write(subdir + "\t")
                if samp in ploidies:
                    summary_out.write(ploidies[samp] + "\t")
                    summary_out.write(purities[samp] + "\t")
                else:
                    summary_out.write("---\t")
                    summary_out.write("---\t")
                summary_out.write(str(goodsamples[samp]) + "\t")
                summary_out.write(str(badsamples[samp]) + "\t")
                if subdir == "tetraploid":
                    summary_out.write("\t")
                summary_out.write(str(sumsca) + "\t\t")
                summary_out.write(str(validatedsca) + "\t\t")
                summary_out.write(str(mixedsca) + "\t\t")
                summary_out.write(str(invalidsca))
                summary_out.write("\n")
            #lsl.createPrintAndSaveHistogram(allpercs, outdir + patient + "_" + str(gamma) + ".txt", .001, xdata="Percent match or anti-match")
#            print "Bad samples:"
#            print badsamples
#            print "Good samples:"
#            print goodsamples
#                for sampIgnore in allpercsminus:
#                    if badsamples[sampIgnore] < 10 or goodsamples[sampIgnore] / badsamples[sampIgnore] > 0.8:
#                        continue
#                    print "Histogram without any data from sample", sampIgnore, ", which had more bad samples than good:"
#                    lsl.createPrintAndSaveHistogram(allpercsminus[sampIgnore], "", 0.001, xdata="Percent match or anti-match")
#                    summary_out.write(gamma + "\t")
#                    summary_out.write(patient + "\twithout_" + sampIgnore + "\t")
#                    numbad = sum(perc < 0.8 for perc in allpercsminus[sampIgnore])
#                    summary_out.write(str(len(allpercsminus[sampIgnore])-numbad) + "\t")
#                    summary_out.write(str(numbad) + "\t")
#                    summary_out.write("\n")
#    patient_file.close()

summary_out.close()


