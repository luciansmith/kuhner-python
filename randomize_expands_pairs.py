#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:02:52 2016

@author: lpsmith
"""

#Create parallel BAF segmentation files to match the CN segmentation files

from __future__ import division
from os import walk
from os.path import isfile

import lucianSNPLibrary as lsl

#Use if you want to analyze the 'rejoined' data:
#expands_directory = "expands_rejoined_input/results_rejoined/"
#outDirectory = "expands_pairs_rejoined_shuffled/"

#Use if you want to analyze the original Xiaohong-segmented data:
expands_directory = "expands_results/"
outDirectory = "expands_pairs_shuffled/"

flist = []
filesets = {}
#SNPfiles.append(["1034", "20008"])
for (_, _, f) in walk(expands_directory):
    flist += f
    
for f in flist:
    if f.find(".sps") != -1 and f.find(".cbs") == -1 and f.find(".spstats") == -1:
        treename = f.replace(".sps", ".tree")
        if not(isfile(expands_directory + treename)):
            continue
        (patient, sample, tag) = f.split("_")
        if not(patient in filesets):
            filesets[patient] = list()
        filesets[patient].append(f)

patientdata = {}
for patient in filesets:
    fileset = filesets[patient]
    if len(fileset) == 1:
        continue
    for outfile in fileset:
        (patient, sample, __) = outfile.split("_")
        efile = open(expands_directory + outfile, "r")
        numpops = 0
        segments = []
        for line in efile:
            if line.find("expands version") != -1:
                continue
            if line.find("chr") != -1:
                numpops = len(line.split()) - 16
                continue
            rowvec = line.rstrip().split()
            scenario = rowvec[14]
            if (scenario == "4"):
                #Don't count the ones that separate the CN and the BAF segments into different populations.
                continue
            chr = int(rowvec[0])
            start = int(float(rowvec[1]))
            end = int(float(rowvec[2]))
            pops = rowvec[15:15+numpops]
            for p in range(0,len(pops)):
                pops[p] = int(pops[p])
            segment = [chr, start, end, pops]
            segments.append(segment)
        if not(patient in patientdata):
            patientdata[patient] = []
        ss = (sample, segments)
        patientdata[patient].append(ss)

#Now we need to unify the segments, first by simply finding where they are:
nrandom = 1000
for patient in patientdata:
    data = patientdata[patient]
    points = lsl.getPoints(data)
    pairs = lsl.calculateOverlap(data, points)
    randompairs = {}
    print "Shuffling data for patient", patient
    for n in range(0,nrandom):
        lsl.randomizePopAssignment(data)
        shuffled_pairs = lsl.calculateOverlap(data, points)
        lsl.collectRandomPairs(randompairs, shuffled_pairs, nrandom)
    outfile = open(outDirectory + patient + "_pairs_shuffled.txt", "w")
    outfile.write("p1_chr\tp1_pos\tp2_chr\tp2_pos\tnsamples\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n")
    print "Writing data for patient", patient
    for pair in pairs:
        samples = pairs[pair]
        outfile.write(str(pair[0]) + "\t")
        outfile.write(str(pair[1]) + "\t")
        outfile.write(str(pair[2]) + "\t")
        outfile.write(str(pair[3]) + "\t")
        outfile.write(str(len(samples)))
        finalvec = lsl.tenvec(nrandom)
        if pair in randompairs:
            finalvec = randompairs[pair]
        for val in finalvec:
            outfile.write(str("\t" + str(val)))
        outfile.write("\n")
    for pair in randompairs:
        if pair in randompairs:
            continue
        outfile.write(str(pair[0]) + "\t")
        outfile.write(str(pair[1]) + "\t")
        outfile.write(str(pair[2]) + "\t")
        outfile.write(str(pair[3]) + "\t")
        outfile.write(str("N/A" + "\t"))
        for val in randompairs[pair]:
            outfile.write("\t" + str(val))
        outfile.write("\n")
        
    outfile.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        