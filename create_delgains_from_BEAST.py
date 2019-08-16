#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 13:49:04 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

import lucianSNPLibrary as lsl

tag = "_g500_better_ploidy/"

beast_input = "BEAST" + tag
character_output = "joint_delsandgains" + tag

somepatientsonly = False
somepatients = ["1001"]

if not(path.isdir(character_output + "/")):
    mkdir(character_output + "/")

def SaveSegment(currseg, changeslist, sample):
    if currseg[3] == 1:
        #Don't save wt segments
        return
    if currseg[3] == "Unknown":
        #Don't save unknown segments
        return
    chr = currseg[0]
    segid = tuple(currseg[1:4])
    if chr not in changeslist:
        changeslist[chr] = {}
    if segid not in changeslist[chr]:
        changeslist[chr][segid] = {}
    changeslist[chr][segid][sample] = True
    

def SetSampleValueIn(changeslist, newsegid, sample):
    if sample in changeslist[newsegid]:
        #already set
        return
    for segid in changeslist:
        if sample in changeslist[segid]:
            if newsegid[2] == "Unknown" or segid[2] == "Unknown":
                changeslist[newsegid][sample] = "Unknown"
                return
            if newsegid[0] == segid[0] and newsegid[1] == segid[1]:
                if newsegid[2] == segid[2]:
                    print("What, exactly, was different?", newsegid, segid)
                if newsegid[2] < segid[2] and newsegid[2] > 1:
                    #We have a '3' compared with a '2': the '3' must have a '2' in its past.
                    changeslist[newsegid][sample] = True
                    return
            elif newsegid[0] >= segid[0] and newsegid[1] <= segid[1] and changeslist[segid][sample]==True:
                #The new segment ID is encompassed by an old one:
                #print "new:", newsegid, "old:", segid
                if newsegid[2] == segid[2]:
                    changeslist[newsegid][sample] = "Unknown"
                    return
                if newsegid[2] > 1 and newsegid[2] < segid[2]:
                    changeslist[newsegid][sample] = "Unknown"
                    return
    changeslist[newsegid][sample] = False

beast_files = []
for (_, _, f) in walk(beast_input):
    beast_files += f

for f in beast_files:
    if f.find("phased") == -1:
        continue
    patient = f.split("_")[0]
    if somepatientsonly and patient not in somepatients:
        continue
#    if (patient != "1042"):
#        continue
    print("Processing patient", patient)
    beastfile = open(beast_input + f, "r")
    samples = beastfile.readline().split()[3:]
    chrs = []
    starts = []
    ends = []
    gainsnlosses = {}
    for sample in samples:
        gainsnlosses[sample] = []
    for line in beastfile:
        linevec = line.split()
        (chr, start, end) = linevec[0:3]
        chrs.append(int(chr))
        starts.append(int(start))
        ends.append(int(end))
        for n in range(len(samples)):
            try:
                gainsnlosses[samples[n]].append(int(linevec[3+n]))
            except:
                gainsnlosses[samples[n]].append("Unknown")
    changeslist = {}
    for sample in samples:
        currseg = [chrs[0], starts[0], ends[0], gainsnlosses[sample][0]]
        for n in range(1, len(chrs)):
            if (chrs[n] != currseg[0] or gainsnlosses[sample][n] != currseg[3]):
                #save the current segment, move on to the next
                SaveSegment(currseg, changeslist, sample)
                currseg = [chrs[n], starts[n], ends[n], gainsnlosses[sample][n]]
            else:
                #Append to the current segment:
                currseg[2] = ends[n]
        SaveSegment(currseg, changeslist, sample)
    for chr in changeslist:
        for segid in changeslist[chr]:
            for sample in samples:
                if sample not in changeslist[chr][segid]:
                    SetSampleValueIn(changeslist[chr], segid, sample)
    delgainOut = open(character_output + f, "w")
    delgainOut.write("chr\tstart\tend\tcopynumber")
    for sample in samples:
        delgainOut.write("\t" + sample)
    delgainOut.write("\n")
    chrs = list(changeslist.keys())
    chrs.sort()
    for chr in chrs:
        segids = list(changeslist[chr].keys())
        segids.sort()
        for segid in segids:
            delgainOut.write(str(chr))
            for val in segid:
                delgainOut.write("\t" + str(val))
            for sample in samples:
                delgainOut.write("\t" + str(changeslist[chr][segid][sample]))
            delgainOut.write("\n")
    delgainOut.close()
    