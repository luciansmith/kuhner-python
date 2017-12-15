#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:02:52 2016

@author: lpsmith
"""

#Create parallel BAF segmentation files to match the CN segmentation files

from __future__ import division
from os import walk
import numpy

import lucianSNPLibrary as lsl

#Use if you want to analyze the 'rejoined' data:
#expands_directory = "expands_rejoined_input/results_rejoined/"
#outDirectory = "expands_pairs_rejoined/"

#Use if you want to analyze the original Xiaohong-segmented data:
expands_directory = "expands_results/"
outDirectory = "expands_pairs/"

def addPopsToList(samepops, point, sample, pops):
    if not(sample in samepops):
        samepops[sample] = {}
    if not (tuple(pops) in samepops[sample]):
        samepops[sample][tuple(pops)] = []
    samepops[sample][tuple(pops)].append(point)

flist = []
filesets = {}
#SNPfiles.append(["1034", "20008"])
for (_, _, f) in walk(expands_directory):
    flist += f
    
for f in flist:
    if f.find(".sps") != -1 and f.find(".cbs") == -1 and f.find(".spstats") == -1:
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
            segment = (chr, start, end, pops)
            segments.append(segment)
        if not(patient in patientdata):
            patientdata[patient] = []
        ss = (sample, segments)
        patientdata[patient].append(ss)

#Now we need to unify the segments, first by simply finding where they are:
for patient in patientdata:
    data = patientdata[patient]
    allsegments = set()
    for segments in data:
        for segment in segments[1]:
            start = (segment[0], segment[1])
            end   = (segment[0], segment[2])
            allsegments.add(start)
            allsegments.add(end)
    points = sorted(list(allsegments))
    #Now go through the segments again, noting the 
    samepops = {}
    for segments in data:
        sample = segments[0]
        seg = segments[1]
        pnum = 0
        snum = 0
        while (pnum < len(points) and snum < len(seg)):
            point = points[pnum]
            chr   = seg[snum][0]
            start = seg[snum][1]
            end   = seg[snum][2]
            pops  = seg[snum][3]
            if chr < point[0]:
                snum = snum+1
            elif chr > point[0]:
                pnum = pnum+1
            elif point[1] == 0 and start == 0:
                #include
                addPopsToList(samepops, point, sample, pops)
                pnum = pnum+1
            elif point[1] < start:
                pnum = pnum+1
            elif point[1] >= end:
                snum = snum+1
            else:
                addPopsToList(samepops, point, sample, pops)
                pnum = pnum+1
    pairs = {}
    for sample in samepops:
        for pops in samepops[sample]:
            points = samepops[sample][pops]
            for p in range(0, len(points)):
                for q in range(p+1, len(points)):
                    pair = (points[p][0], points[p][1], points[q][0], points[q][1])
                    if not(pair in pairs):
                        pairs[pair] = []
                    pairs[pair].append(sample)
    outfile = open(outDirectory + patient + "_pairs.txt", "w")
    outfile.write("p1_chr\tp1_pos\tp2_chr\tp2_pos\tnsamples\tsamples\n")
    for pair in pairs:
        samples = pairs[pair]
        outfile.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" + str(pair[2]) + "\t" + str(pair[3]) + "\t" + str(len(samples)))
        for sample in samples:
            outfile.write("\t" + sample)
        outfile.write("\n")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        