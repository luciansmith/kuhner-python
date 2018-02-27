# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 13:36:23 2018

@author: Lucian
"""
from __future__ import division
from os import walk
from os import path
from os import readlink
from os.path import isfile
from copy import deepcopy

import numpy
import math
import matplotlib.pyplot as plt

import lucianSNPLibrary as lsl

def addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca):
    if chr == "":
        return
    if chr == "23":
        return
    if chr == "24":
        return
    if chr=="0":
        #print("There's a chromosome zero segment, weirdly:", full_sample, chr, start, end)
        return
    svec = full_sample.split("-")
    if len(svec) < 4:
        svec = full_sample.split("_")
    patient = svec[0]
    sample = svec[0] + "_" + svec[1]
    if patient not in Xiaohong_segments:
        Xiaohong_segments[patient] = {}
        totsca[patient] = {}
    if sample not in Xiaohong_segments[patient]:
        Xiaohong_segments[patient][sample] = {}
        totsca[patient][sample] = set()
    if chr not in Xiaohong_segments[patient][sample]:
        Xiaohong_segments[patient][sample][chr] = []
    start = int(start)
    end = int(end)
    Xiaohong_segments[patient][sample][chr].append((start, end))
    totsca[patient][sample].add((chr, start, end))
    if "overall" not in totsca[patient]:
        totsca[patient]["overall"] = set()
    totsca[patient]["overall"].add((chr, start, end))
    #print("Adding", chr, start, end)

def readXiaohongWGSLOHFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (__, full_sample, __, chr, start, end) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca)
    xfile.close()

def readXiaohong1MLOHFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    for line in xfile:
        (full_sample, chr, start, end, __, __) = line.split()
        addXiaohongSegment(Xiaohong_segments, full_sample, chr, start, end, totsca)
    xfile.close()

def readXiaohongCopynumFile(f, Xiaohong_segments, totsca):
    print("reading", f)
    xfile = open(f, "r")
    lastseg = ["", 0, 0, 0]
    for line in xfile:
        lvec = line.split()
        if len(lvec) == 7:
            (full_sample, chr, start, end, __, __, code) = line.split()
        elif len(lvec) == 10:
            (__, full_sample, __, chr, start, end, __, __, __, code) = line.split()
        else:
            print("Incorrect line length:", line)
            continue
        start =int(start)
        end = int(end)
        if code == "Double_d" or code=="Balanced_gain" or code == "34" or code=="41":
            #balanced gain or loss: continue
            addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)
            lastseg = ["", 0, 0, 0]
            continue
        if chr == lastseg[0] and code == lastseg[3] and start-lastseg[2] < 5000:
            #print("Combining", chr, str(lastseg[1]), str(lastseg[2]), str(lastseg[3]), "with", start, end, code)
            lastseg[2] = end
        else:
            addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)
            lastseg = [chr, start, end, code]
    xfile.close()
    addXiaohongSegment(Xiaohong_segments, full_sample, lastseg[0], lastseg[1], lastseg[2], totsca)

def readAllXiaohongSegmentation():
    print("Reading all Xiaohong segmentation.")
    Xiaohong_segments = {}
    totsca = {}
    files = []
    Xdir_WGS = "Xiaohong_WGS_segmentation/"
    Xdir_1M = "CN_Xiaohong_segmentation/"
    for (__, __, f) in walk(Xdir_WGS):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("test") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohongWGSLOHFile(Xdir_WGS + f, Xiaohong_segments, totsca)
        else:
            readXiaohongCopynumFile(Xdir_WGS + f, Xiaohong_segments, totsca)
    files = []
    for (__, __, f) in walk(Xdir_1M):
        files += f
    for f in files:
        if f.find("read") != -1:
            continue
        if f.find("LOH") != -1:
            readXiaohong1MLOHFile(Xdir_1M + f, Xiaohong_segments, totsca)
        else:
            readXiaohongCopynumFile(Xdir_1M + f, Xiaohong_segments, totsca)
    return Xiaohong_segments, totsca

Xsegs, totsca = readAllXiaohongSegmentation()
alldiffs = []
for patient in Xsegs:
    for sample in Xsegs[patient]:
        for chr in Xsegs[patient][sample]:
            Xsegs[patient][sample][chr].sort()
            for i in range(1,len(Xsegs[patient][sample][chr])):
                endlast = Xsegs[patient][sample][chr][i-1][1]
                startnext = Xsegs[patient][sample][chr][i][0]
                if endlast < startnext:
                    alldiffs.append(numpy.log10(startnext-endlast))

lsl.createPrintAndSaveHistogram(alldiffs, "", 0.01)