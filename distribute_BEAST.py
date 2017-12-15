#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:44:27 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

import numpy

import lucianSNPLibrary as lsl

tag = "_25M_v2_joined_best/"

BEAST_indir = "BEAST"+ tag
BEAST_outdir = "BEAST" + "_allcombos" + tag

if not(path.isdir(BEAST_outdir + "/")):
    mkdir(BEAST_outdir + "/")

def writeAllOutfiles(oneset, pairs, index, segments, data, f):
    if len(pairs) == 0:
        #print oneset, pairs, index
        fname = BEAST_outdir + f
        fname = fname.replace("_phased", "_" + str(index) + "_phased")
        outfile = open(fname, "w")
        outfile.write("Chr\tStart\tEnd\tnLogR\tNBAF")
        for sample in oneset:
            outfile.write("\t" + sample)
        outfile.write("\n")
        for n in range(len(segments)):
            segment = segments[n]
            for segbit in segment:
                outfile.write(segbit + "\t")
            for sample in oneset:
                outfile.write(data[sample][n] + "\t")
            outfile.write("\n")
        outfile.close()
        return index+1
    else:
        newpair = pairs[0]
        pairs.remove(newpair)
        index = writeAllOutfiles(oneset, pairs[:], index, segments, data, f)
        orig = newpair.replace("b","")
        oneset = [newpair if x==orig else x for x in oneset]
        index = writeAllOutfiles(oneset, pairs[:], index, segments, data, f)
    return index

def writeEdgeOutfiles(oneset, pairs, index, segments, data, f):
    if len(pairs) == 0:
        #print oneset, pairs, index
        fname = BEAST_outdir + f
        fname = fname.replace("_phased", "_" + str(index) + "_phased")
        outfile = open(fname, "w")
        outfile.write("Chr\tStart\tEnd\tnLogR\tNBAF")
        for sample in oneset:
            outfile.write("\t" + sample)
        outfile.write("\n")
        for n in range(len(segments)):
            segment = segments[n]
            for segbit in segment:
                outfile.write(segbit + "\t")
            for sample in oneset:
                outfile.write(data[sample][n] + "\t")
            outfile.write("\n")
        outfile.close()
    else:
        writeEdgeOutfiles(oneset, [], "diploidier", segments, data, f)
        for newpair in pairs:
            orig = newpair.replace("b","")
            oneset = [newpair if x==orig else x for x in oneset]
        writeEdgeOutfiles(oneset, [], "tetraploidier", segments, data, f)
    return index

BEASTlist = []
for (_, _, f) in walk(BEAST_indir):
    BEASTlist += f

for f in BEASTlist:
    if f.find("phased") == -1:
        continue
    infile = open(BEAST_indir + f, "r")
    headers = infile.readline().split()
    samples = headers[5:]
    data = {}
    for sample in samples:
        data[sample] = []
    segments = []
    for line in infile:
        linevec = line.split()
        segments.append(linevec[0:5])
        for n in range(len(samples)):
            data[samples[n]].append(linevec[5+n])
    pairs = []
    originals = []
    for sample in samples:
        if sample.find("b") != -1:
            pairs.append(sample)
        else:
            originals.append(sample)
    
    index = 0
    writeEdgeOutfiles(originals, pairs, "allknown", segments, data, f)
    