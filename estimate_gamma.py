# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:50:09 2017

@author: lpsmith
"""

from __future__ import division
from os import walk

import numpy

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

nygcsamples = ["23599","23602","23605","24660","24409","24412","24810","24813","24139","24145","24148","24900",]

outfile = "joint_summary_gamma_test.txt"
out = open(outfile, "w")
out.write("patient\tgamma\ttotlogr\tploidysum\tploidysum1m\tploidysum25m\tnbreaks\tnusedbreaks\tnusedbreaks1m\tnusedbreaks25m\tfrequencies\n")
for tag in ["50", "100", "150", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900", "1000", "1200", "1400", "1600", "2000"]:

    #print tag
    root_dir = "gamma_test_output/pASCAT_input_g" + tag + "/"

    segs_init = {}
    segs_final = {}


    files = []
    for (__, __, f) in walk(root_dir + "unconstrained/"):
        files = f
    for f in files:
        if f.find("fcn_ascat_segments") != -1:
            patient = f.split("_")[0]
            #print "Patient", patient, "using a gamma of", tag
            if patient not in segs_init:
                segs_init[patient] = {}
                segs_final[patient] = {}
            initialsegs = open(root_dir + patient + "_copynumber_segments.txt")
            totlogr = 0
            for line in initialsegs:
                if line.find("Chr") != -1:
                    continue
                (chr, start, end, nlogr, nbaf) = line.split()
                if chr not in segs_init[patient]:
                    segs_init[patient][chr] = []
                segs_init[patient][chr].append((start, end))
                totlogr += int(nlogr)

            finalsegs = open(root_dir + "unconstrained/" + f, "r")
            prevsample = 0
            prevchr = 0
            prevnA = 0
            prevnB = 0
            for line in finalsegs:
                if line.find("Chr") != -1:
                    continue
                lvec = line.split('"')
                (patient, sample) = lvec[3].split("_")
                (__, __, chr, start, end, nprobes, mbaf, logr, nA, nB) = line.split()
                if sample not in segs_final[patient]:
                    segs_final[patient][sample] = {}
                if chr not in segs_final[patient][sample]:
                    segs_final[patient][sample][chr] = set()
                segs_final[patient][sample][chr].add(start)
                segs_final[patient][sample][chr].add(end)
                if (prevsample == sample and prevchr == chr and prevnA == nA and prevnB == nB):
                    print "Error: two sequential segments with same nA and nB:"
                    print sample, chr, start, end, nA, nB

            ploidyfile = open(root_dir + "unconstrained/" + patient + "_fcn_ascat_ploidy.txt")
            ploidysum = 0
            ploidysum1m = 0
            ploidysum25m = 0
            for line in ploidyfile:
                if line.find("x") != -1:
                    continue
                (__, sample, ploidy) = line.split('"')
                ploidy = float(ploidy)
                ploidysum += ploidy
                sample = sample.split("_")[1]
                if (sample in nygcsamples):
                    ploidysum25m += ploidy
                else:
                    ploidysum1m += ploidy

            #Now analyze the overlap:
            nbreaks = 0
            nusedbreaks = 0
            n1musedtot = 0
            n25musedtot = 0
            usedfreqs = {}
            for chr in segs_init[patient]:
                for n in range(1,len(segs_init[patient][chr])):
                    break1 = segs_init[patient][chr][n][0]
                    break2 = segs_init[patient][chr][n-1][1]
                    nbreaks += 1
                    nused = 0
                    n25mused = 0
                    n1mused = 0
                    for sample in segs_final[patient]:
                        segset = segs_final[patient][sample][chr]
                        if break1 in segset or break2 in segset:
                            nused += 1
                            if sample in nygcsamples:
                                n25mused += 1
                            else:
                                n1mused += 1
                    if nused > 0:
                        nusedbreaks += 1
                    if n1mused > 0:
                        n1musedtot += 1
                    if n25mused > 0:
                        n25musedtot += 1
                    if nused not in usedfreqs:
                        usedfreqs[nused] = 0
                    usedfreqs[nused] += 1
            out.write(patient + "\t")
            out.write(tag + "\t")
            out.write(str(totlogr) + "\t")
            out.write(str(ploidysum) + "\t")
            out.write(str(ploidysum1m) + "\t")
            out.write(str(ploidysum25m) + "\t")
            out.write(str(nbreaks) + "\t")
            out.write(str(nusedbreaks) + "\t")
            out.write(str(n1musedtot) + "\t")
            out.write(str(n25musedtot) + "\t")
            out.write(str(usedfreqs) + "\n")

            #print "nbreaks:", nbreaks
            #print "nusedbreaks:", nusedbreaks
            #print "usedfreqs", usedfreqs
out.close()