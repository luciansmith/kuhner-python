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

def combinePopsInList(combinedpops, point, sample, pops):
    if not(point in combinedpops):
        combinedpops[point] = {}
    combinedpops[point][sample] = pops


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
for patient in patientdata:
    data = patientdata[patient]
    points = lsl.getPoints(data)
    combinedpops = {}
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
                combinePopsInList(combinedpops, point, sample, pops)
                pnum = pnum+1
            elif point[1] < start:
                pnum = pnum+1
            elif point[1] >= end:
                snum = snum+1
            else:
                combinePopsInList(combinedpops, point, sample, pops)
                pnum = pnum+1
    foundfour = 0
    notfour = 0
    for n in range(0, len(points)-1):
        for s in range(n+1, len(points)):
            gametes1 = []
            gametes2 = []
            if not(points[n] in combinedpops):
                continue
            if not(points[s] in combinedpops):
                continue
            for sample in combinedpops[points[n]]:
                if not(sample in combinedpops[points[s]]):
                    continue
                gametes1 += combinedpops[points[n]][sample]
                gametes2 += combinedpops[points[s]][sample]
            if (len(gametes1) != len(gametes2)):
                print "Error!"
            g4test = set()
            g4test.add("zerozero") #Because we know that the ancestral allele was WT.
            for g in range(0,len(gametes1)):
                g1 = gametes1[g]
                g2 = gametes2[g]
                if (g1==1):
                    if (g2==1):
                        g4test.add("oneone")
                    else:
                        g4test.add("onezero")
                else:
                    if (g2==1):
                        g4test.add("zeroone")
                    else:
                        g4test.add("zerozero")
            if len(g4test) > 3:
                foundfour += 1
            else:
                notfour += 1
    print "Patient", patient, "4-gamete total =", foundfour, "not four = ", notfour


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        