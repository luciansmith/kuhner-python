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

expands_directory = "CN_calc_log2rs/"
flist = [] #["999_20716_avglog2rs.txt"]
for (_, _, f) in walk(expands_directory):
    flist += f
    
for f in flist:
    print "Reading", f
    if f.find("avglog2rs") == -1:
        continue
    filevec = f.split("_")
    patient = filevec[0]
    sample = filevec[1]
    infile = open(expands_directory + f, "r")
    prevline = []
    for line in infile:
        thisline = line.split()
        if thisline[0:2] == prevline[0:2]:
            print "Duplicate segments in", patient, sample, ":"
            print prevline
            print thisline
        prevline = thisline
