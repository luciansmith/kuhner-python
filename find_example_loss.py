# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:27:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile
import numpy

# read the processed data from the file Xiaohong sent me, which has segmentation and non-2N log2R values

#filenames = ["1032_20834_78Y_34.txt"]
filenames = []
for (_, _, f) in walk("CN_float_Xiaohong/"):
    filenames += f
    break        

outfilename = "short_segments/all_short.txt"
outfile = open(outfilename, "w")
outline = "patient\tsample\tchr\tstart\tend\tXiao_log2r\tcall\tnum_log2rs\tavg(log2r)\tboundary log2rs\n"
outfile.write(outline)

for file in filenames:
    split= file.split("_")
    if (len(split) < 4):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    segmented_file = open("CN_float_Xiaohong/" + file, "r")
    log2r_filename = "CN_filtered_data/" + patient + "_" + sample + "_copynumber_all.txt"
    if (not(isfile(log2r_filename))):
        continue
    log2r_file = open(log2r_filename, "r")
        
    #Go through the segmentation file and build up a database of segments.
    lastchr = 1
    startpos = 0
    inf = float('infinity')
    segments = []
    
    for line in segmented_file:
        (chr, start, end, log2r, call) = line.rstrip().split()
        if(chr == "chr"):
            continue
        chr = int(chr)
        if (chr >= 23):
            #Don't do stuff with the sex chromosomes
            continue
        start = int(start)
        end = int(end)
        log2r = float(log2r)
        while (chr > lastchr):
            segments.append((lastchr, startpos, inf, 0, "wt"))
            lastchr = lastchr+1
            startpos = 0
        if (start > startpos):
            segments.append((chr, startpos, start, 0, "wt"))
            
        segments.append((chr, start, end, log2r, call))
        startpos = end
        
    while (lastchr <= 22):
        segments.append((lastchr, startpos, inf, 0, "wt"))
        startpos = 0
        lastchr = lastchr+1
            
    #print segments
            
    seglog2rs = {}
    for segment in segments:
        seglog2rs[segment] = []
    
    for line in log2r_file:
        (id, chr, pos, log2r) = line.rstrip().split()
        if (id=="SNPid"):
            continue
        if (log2r == "?"):
            continue
        chr = int(chr)
        if (chr==0):
            continue
        if (chr >= 23):
            continue
        pos = int(pos)
        if (pos==0):
            continue
        log2r = float(log2r)
        found = False
        for segment in segments:
            if (chr == segment[0] and
                pos >  segment[1] and
                pos <= segment[2]):
                    oldarray = seglog2rs[segment]
                    oldarray.append(log2r)
                    seglog2rs[segment] = oldarray
                    found = True
        if (found == False):
            print "Couldn't find segment for chr ", chr, ", pos ", pos, "in file ", log2r_filename
        
    #Write out the segment data:
    
    for s in range(1,len(segments)-1):
        segment = segments[s]
        log2rs = seglog2rs[segment]
        numlogs = len(log2rs)
        avg = numpy.mean(log2rs)
        if (numlogs <= 20 and avg < -.8 and avg > -1.3):
            if segment[0] != segments[s-1][0]:
                continue
            if segment[0] != segments[s+1][0]:
                continue
            #The segment in question has same-chromosome neighbors!
            log2rs = seglog2rs[segments[s-1]]
            outline = patient + "\t" + sample + "\t" + str(segments[s-1][0]) + "\t" + str(segments[s-1][1]) + "\t" + str(segments[s-1][2]) + "\t" + str(segments[s-1][3]) + "\t" + str(segments[s-1][4]) + "\t" + str(len(log2rs)) + "\t" + str(numpy.mean(log2rs))
            lsize = len(log2rs)
            for l in range(max(0,lsize-10),lsize):
                outline += "\t" + str(log2rs[l])
            outline += "\n"
            outfile.write(outline)

            log2rs = seglog2rs[segments[s]]
            outline = patient + "\t" + sample + "\t" + str(segments[s][0]) + "\t" + str(segments[s][1]) + "\t" + str(segments[s][2]) + "\t" + str(segments[s][3]) + "\t" + str(segments[s][4]) + "\t" + str(len(log2rs)) + "\t" + str(numpy.mean(log2rs))
            lsize = len(log2rs)
            for log2r in log2rs:
                outline += "\t" + str(log2r)
            outline += "\n"
            outfile.write(outline)

            log2rs = seglog2rs[segments[s+1]]
            outline = patient + "\t" + sample + "\t" + str(segments[s+1][0]) + "\t" + str(segments[s+1][1]) + "\t" + str(segments[s+1][2]) + "\t" + str(segments[s+1][3]) + "\t" + str(segments[s+1][4]) + "\t" + str(len(log2rs)) + "\t" + str(numpy.mean(log2rs))
            lsize = len(log2rs)
            for l in range(0,min(10, lsize)):
                outline += "\t" + str(log2rs[l])
            outline += "\n"
            outfile.write(outline)

            outfile.write("\n")

outfile.close()
