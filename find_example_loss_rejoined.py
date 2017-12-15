# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:27:59 2016

@author: lpsmith
""" 

from __future__ import division
from os import walk
from os.path import isfile
import numpy

rejoin = ""
#rejoin = "_rejoined"

# read the processed data from the file Xiaohong sent me, which has segmentation and non-2N log2R values

#filenames = ["660_19272_avglog2rs.txt"]
filenames = []
for (_, _, f) in walk("CN_calc_log2rs" + rejoin + "/"):
    filenames += f
    break        

outfilename = "short_segments/all_short" + rejoin + ".txt"
outfile = open(outfilename, "w")
outline = "patient\tsample\tchr\tstart\tend\tXiao_log2r\tcall\tnum_log2rs\tavg(log2r)\tstdev\n"
outfile.write(outline)

for file in filenames:
    split= file.split("_")
    if (len(split) < 3):
        continue
    patient = split[0]
    sample = split[1]
    if (patient[0] < '0' or patient[0] > '9'):
        continue
    segmented_file = open("CN_calc_log2rs" + rejoin + "/" + file, "r")
        
    two_prev = ""
    prev_line = ""
    write_next = False
    for line in segmented_file:
        if (line.find("chr") != -1):
            continue
        this_line = patient + "\t" + sample + "\t" + line
        if write_next:
            outfile.write(two_prev)
            outfile.write(prev_line)
            outfile.write(this_line)
            outfile.write("\n")
            write_next = False
        (chr, start, end, xLog2r, call, nlog2r, log2r, stdev) = line.rstrip().split()
        chr = int(chr)
        if (chr >= 23):
            #Don't do stuff with the sex chromosomes
            continue
        nlog2r = int(nlog2r)
        if (nlog2r <= 12 and nlog2r > 5 and start != "0" and end != "inf"):
            write_next = True
        else:
            write_next = False
        two_prev = prev_line
        prev_line = this_line

outfile.close()
