#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 13:27:15 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile

import lucianSNPLibrary as lsl

expands_directory = "expands_full_input_alldupes/"
expands_output = "expands_full_input_alldupes_noLOH/"

LOHfilename = "CN_Xiaohong_segmentation/LOH_for_Lucian.txt"

LOHfile = open(LOHfilename, "r")
LOHsegs = set()
for line in LOHfile:
    (id, chr, start, end, __, __) = line.split()
    idvec = id.split("_")
    patient = idvec[0]
    sample = idvec[1]
    LOHsegs.add((patient, sample, chr, start, end))
    LOHsegs.add((patient, sample, chr, "0", end))

flist = []
for (_, _, f) in walk(expands_directory):
    flist += f
    
for f in flist:
    if f.find("BAF") == -1 and f.find("CN") == -1:
        continue
    filevec = f.split("_")
    patient = filevec[0]
    sample = filevec[1]
    infile = open(expands_directory + f, "r")
    outfile = open(expands_output + f, "w")
    for line in infile:
        invec = line.split()
        chr = invec[0]
        start = invec[1]
        end = invec[2]
        if (patient, sample, chr, start, end) in LOHsegs:
            print "Not writing out", patient, sample, chr, start, end
            continue
        #Write it out
        outfile.write(line)