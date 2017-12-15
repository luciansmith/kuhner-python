# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 13:57:27 2016

@author: lpsmith
"""

from sets import Set
import lucianSNPLibrary as lsl

LOHfilename = "CN_Xiaohong_segmentation/LOH_for_Lucian.txt"
LOHfile = open(LOHfilename, "r")
segs = {}
for line in LOHfile:
    (id, chr, start, end, __, __) = line.split()
    if not(id in segs):
        segs[id] = []
    segs[id].append([int(chr), int(start), int(end), "0.0", "LOH"])
    

#Then read the file Xiaohong sent with the log2R data in it
CNfilename = "CN_Xiaohong_segmentation/copy_number_log2_Lucian.txt"
CNfile = open(CNfilename,"r")

ids = Set()
for line in CNfile:
    (id, chr, start, end, __, ln2r, call) = line.rstrip().split()
    if not(id in segs):
        segs[id] = []
    segs[id].append([int(chr), int(start), int(end), ln2r, call])

for id in segs:
    outvec = segs[id]
    outvec = lsl.resegmentLOHes(outvec)
    outfile = open("CN_float_Xiaohong/" + id + ".txt", "w")
    outfile.write("chr\tstartpos\tendpos\tCN_Estimate\tCall\n")
    for seg in outvec:
        outfile.write(str(seg[0]) + "\t" + str(seg[1]) + "\t" + str(seg[2]) + "\t" + seg[3] + "\t" + seg[4] + "\n")
    
#print ids
#print "total number of IDs = " + str(len(ids))
