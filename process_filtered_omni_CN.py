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

discrepancies = open("discrepancies.txt", "w")
outline = "chr\tstart\tend\tXiao_log2r\tcall\tnum_log2rs\tavg(log2r)\tstdev\n"
discrepancies.write(outline)

for file in filenames:
    fvec= file.split("_")
    if (len(fvec) < 4):
        continue
    patient = fvec[0]
    sample = fvec[1]
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
    outfilename = "CN_calc_log2rs/" + patient + "_" + sample + "_avglog2rs.txt"
    print "Writing data for", patient, sample
    outfile = open(outfilename, "w")
    outline = "chr\tstart\tend\tXiao_log2r\tcall\tnum_log2rs\tavg(log2r)\tstdev\n"
    outfile.write(outline)
    
    for segment in segments:
        log2rs = seglog2rs[segment]
        invlog2rs = []
        for log2r in log2rs:
            invlog2rs.append(numpy.power(2,log2r))
        numlogs = len(log2rs)
        meanval = "?"
        discrepancy = False
        if (numlogs > 0):
            meanval = str(numpy.mean(log2rs))
            stdev = str(numpy.std(log2rs))
            try:
                if abs(segment[3] - meanval > 0.000001):
                    discrepancy = True
            except:
                discrepancy = False
        outline = str(segment[0]) + "\t" + str(segment[1]) + "\t" + str(segment[2]) + "\t" + str(segment[3]) + "\t" + str(segment[4]) + "\t" + str(numlogs) + "\t" + meanval + "\t" + stdev + "\n"
        outfile.write(outline)
        if (discrepancy):
            discrepancies.write(outline)
    outfile.close()
discrepancies.close()