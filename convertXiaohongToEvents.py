#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:43:25 2017

@author: lpsmith
"""

#Take *all* BAF and CN data and create expands input.

from __future__ import division
from os import walk
from os import path
from os import readlink
from os.path import isfile
from os import mkdir
from shutil import copy2 as copy
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt

import lucianSNPLibrary as lsl

#Use this value to set up whether to use the 'rejoined' segments or not

outdir = "xiaohong_events/"

onlysomepatients = False
#somepatients = ["55", "59", "74", "591", "595", "597"]
#somepatients = ["403", "512", "568", "852", "572"]
somepatients = ["772"]

#chrs = [12]
chrs = range(1,23)

def addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, call):
    if chr == "":
        return
    if chr == "23":
        return
    if chr == "24":
        return
    if chr=="0":
        #print("There's a chromosome zero segment, weirdly:", full_sample, chr, start, end)
        return
    chr = int(chr)
    if chr not in chrs:
        return
    start = int(start)
    end = int(end)
    added = False
    svec = full_sample.split("-")
    if len(svec) < 4:
        svec = full_sample.split("_")
    patient = svec[0]
    sample = svec[1]
    if onlysomepatients and patient not in somepatients:
        return
    if patient not in Xiaohong_segments:
        Xiaohong_segments[patient] = {}
    if sample not in Xiaohong_segments[patient]:
        Xiaohong_segments[patient][sample] = {}
    if chr not in Xiaohong_segments[patient][sample]:
        Xiaohong_segments[patient][sample][chr] = []
    else:
        lastseg = Xiaohong_segments[patient][sample][chr][-1]
        if lastseg[2] == call and start - lastseg[1] < 5000:
            lastseg[1] = end
            added = True
    if not added:
        Xiaohong_segments[patient][sample][chr].append([start, end, call])

def readXiaohongWGSLOHFile(f, Xiaohong_segments):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (__, full_sample, __, chr, start, end) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, "CNLOH")
    xfile.close()

def readXiaohong1MLOHFile(f, Xiaohong_segments):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (full_sample, chr, start, end, __, __) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, "CNLOH")
    xfile.close()

def readXiaohongCopynumFile(f, Xiaohong_segments):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        lvec = line.split()
        if len(lvec) == 7:
            (full_sample, chr, start, end, __, __, call) = line.split()
        elif len(lvec) == 10:
            (__, full_sample, __, chr, start, end, __, __, __, call) = line.split()
        else:
            print("Incorrect line length:", line)
            continue
        if call=="41":
            call = "Double_d"
        elif call=="34":
            call = "Balanced_gain"
        elif call=="22" or call=="23":
            call = "Loss"
        elif call=="32" or call=="33":
            call = "Gain"
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, call)
    xfile.close()

def addCentromereToXiaohong(xsegs):
    for patient in xsegs:
        for sample in xsegs[patient]:
            for chr in xsegs[patient][sample]:
                (start, end) = lsl.getChromosomeCentromere(int(chr))
                xsegs[patient][sample][chr].append([start, end, "Centromere"])

def readAllXiaohongSegmentation():
    print("Reading all Xiaohong segmentation.")
    Xiaohong_segments = {}
    files = []
    Xdir_WGS = "Xiaohong_WGS_segmentation/"
#    Xdir_1M = "CN_Xiaohong_segmentation/"
    for (__, __, f) in walk(Xdir_WGS):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("test") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohongWGSLOHFile(Xdir_WGS + f, Xiaohong_segments)
        else:
            readXiaohongCopynumFile(Xdir_WGS + f, Xiaohong_segments)
#    files = []
#    for (__, __, f) in walk(Xdir_1M):
#        files += f
#    for f in files:
#        if f.find("read") != -1:
#            continue
#        if f.find("LOH") != -1:
#            readXiaohong1MLOHFile(Xdir_1M + f, Xiaohong_segments)
#        else:
#            readXiaohongCopynumFile(Xdir_1M + f, Xiaohong_segments)
    addCentromereToXiaohong(Xiaohong_segments)
    return Xiaohong_segments

def regularizeXiaohongSegmentation(xsegs):
    for patient in xsegs:
        for sample in xsegs[patient]:
            for chr in xsegs[patient][sample]:
                xlist = xsegs[patient][sample][chr]
                xlist.sort()
                newlist = []
                prevx = xlist[0]
                for n in range(1,len(xlist)):
                    xnew = xlist[n]
                    if xnew[0] == prevx[1]:
                        xnew[0] = xnew[0]+1
                    if xnew[0] > prevx[1]:
                        newlist.append(prevx)
                        prevx = xnew
                    else:
                        if prevx[2] == "Centromere":
                            continue
                        if xnew[2] == "Centromere":
                            prevx = xnew
                            continue
                        if prevx[2] == "CNLOH":
                            continue
                        if xnew[2] == "CNLOH":
                            prevx = xnew
                            continue
                        print("Overlapping internals???")
                        assert(False)
                newlist.append(prevx)
                xsegs[patient][sample][chr] = newlist

def addWTToSegs(xsegs):
    for patient in xsegs:
        for sample in xsegs[patient]:
            for chr in chrs:
                if chr not in xsegs[patient][sample]:
                    xsegs[patient][sample][chr] = [[0, lsl.getChromosomeMax(chr), "wt"]]
                    continue
                xlist = xsegs[patient][sample][chr]
                xlist.sort()
                newlist = []
                prevx = xlist[0]
                if (prevx[0] > 0):
                    newlist.append([0, prevx[0]-1, "wt"])
                    prevx = newlist[0]
                for n in range(0,len(xlist)):
                    xnew = xlist[n]
                    gap = xnew[0] - prevx[1]
                    if gap>10000:
                        newlist.append((prevx[1]+1, xnew[0]-1, "wt"))
                    newlist.append(xnew)
                    prevx = xnew
                end = lsl.getChromosomeMax(int(chr))
                gap = end - prevx[1]
                if gap>10000:
                    newlist.append([prevx[1]+1, end, "wt"])
                xsegs[patient][sample][chr] = newlist

def printSegs(xsegs):
    for patient in xsegs:
        for sample in xsegs[patient]:
            outfile = open(outdir + patient + "_" + sample + "_xcall.txt", "w")
            outfile.write("patient\tsample\tchr\tstartpos\tendpos\tintA\tintB\n")
            for chr in chrs:
                for seg in xsegs[patient][sample][chr]:
                    outfile.write(patient)
                    outfile.write("\t" + sample)
                    outfile.write("\t" + str(chr))
                    outfile.write("\t" + str(seg[0]))
                    outfile.write("\t" + str(seg[1]))
                    call = seg[2]
                    if (call=="wt") or call=="Centromere":
                        outfile.write("\t1\t1\n")
                    elif (call=="CNLOH"):
                        outfile.write("\t2\t0\n")
                    elif (call=="Double_d"):
                        outfile.write("\t0\t0\n")
                    elif (call=="Balanced_gain"):
                        outfile.write("\t2\t2\n")
                    elif (call=="Loss"):
                        outfile.write("\t1\t0\n")
                    elif (call=="Gain"):
                        outfile.write("\t2\t1\n")
        

Xiaohong_segments = readAllXiaohongSegmentation()
regularizeXiaohongSegmentation(Xiaohong_segments)
addWTToSegs(Xiaohong_segments)
printSegs(Xiaohong_segments)

