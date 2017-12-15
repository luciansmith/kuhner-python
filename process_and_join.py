# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:27:59 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile
import numpy

def shouldMerge(seg1, seg2, seglog2rs):
    if (len(seg2) == 0):
        return False
    nstds = 1
    if seg1[0] != seg2[0]:
        return False
    if seg1[4] != seg2[4]:
        return False
    seg1l2r = seglog2rs[tuple(seg1)]
    seg2l2r = seglog2rs[seg2]
    mean1 = numpy.mean(seg1l2r)
    mean2 = numpy.mean(seg2l2r)
    std1 = numpy.std(seg1l2r)
    std2 = numpy.std(seg2l2r)
    if (abs(mean1-mean2) > nstds*(std1+std2)):
        return False
    print "Joining segments", seg1, "and", seg2
    return True
        

# read the processed data from the file Xiaohong sent me, which has segmentation and non-2N log2R values

#filenames = ["1006_21104_190V_38.txt"]
filenames = []
for (_, _, f) in walk("CN_float_Xiaohong/"):
    filenames += f
    break        

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
        
    while (lastchr <= 24):
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
        pos = int(pos)
        if (pos==0):
            continue
        log2r = float(log2r)
        found = False
        for segment in segments:
            if (chr == segment[0] and
                pos >  segment[1] and
                pos <= segment[2]):
                    seglog2rs[segment].append(log2r)
                    found = True
                    break
        if (found == False):
            print "Couldn't find segment for chr ", chr, ", pos ", pos, "in file ", log2r_filename
        
    #Write out the segment data:
    outfilename = "CN_calc_log2rs_rejoined/" + patient + "_" + sample + "_avglog2rs.txt"
    outfile = open(outfilename, "w")
    outline = "chr\tstart\tend\tXiao_log2r\tcall\tnum_log2rs\tavg(log2r)\tstdev\n"
    outfile.write(outline)
    
    s = 0
    for segment in segments:
        segment = list(segment)
        log2rs = seglog2rs[segments[s]]
        numlogs = len(log2rs)
        meanval = "?"
        stdev = "?"
        if (numlogs > 0):
            meanval = str(numpy.mean(log2rs))
            stdev = str(numpy.std(log2rs))
            if (s < len(segments)-1):
                nextseg = segments[s+1]
                nextlog2rs = seglog2rs[nextseg]
                while shouldMerge(segment, nextseg, seglog2rs):
                    log2rs += nextlog2rs
                    log2rs = seglog2rs[tuple(segment)]
                    segment[2] = nextseg[2]
#                    segment[4] = "[merged]"
                    seglog2rs[tuple(segment)] = log2rs
                    del segments[s+1]
                    nextseg = []
                    nextlog2rs = []
                    if (s < len(segments)-1):
                        nextseg = segments[s+1]
                        nextlog2rs = seglog2rs[nextseg]
                meanval = str(numpy.mean(log2rs))
                stdev = str(numpy.std(log2rs))
        s += 1
        
        outline = str(segment[0]) + "\t" + str(segment[1]) + "\t" + str(segment[2]) + "\t" + str(segment[3]) + "\t" + str(segment[4]) + "\t" + str(len(log2rs)) + "\t" + meanval + "\t" + stdev + "\n"
        outfile.write(outline)
    outfile.close()
