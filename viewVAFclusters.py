#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 12:45:13 2019

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy
from ete3 import Tree
from cycler import cycler

import numpy
import math
import matplotlib.pyplot as plt
import csv

import lucianSNPLibrary as lsl

onlysomepatients = True
somepatients = ["391"]

VAFdir = "VAFclusters/"
outdir = "VAFclusters_pngs/"

if not path.isdir(outdir):
    mkdir(outdir)

VAFfiles = []
for __, _, files in walk(VAFdir):
    VAFfiles += files

def sortLabels(labels):
    newlist = []
    sublist = []
    for label in labels:
        if len(label.split(", '"))==1:
            sublist.append(label)
    sublist.sort()
    newlist.extend(sublist)
    for n in range(9,1,-1):
        for dels in range(0, 9):
            sublist = []
            for label in labels:
                delcount = label.count("-")
                if len(label.split(", '"))==n and delcount == dels:
                    sublist.append(label)
            sublist.sort()
            newlist.extend(sublist)
    return newlist

for file in VAFfiles:
    if "_VAFs" not in file:
        continue
    (patient, sample) = file.split("_")[0:2]
    if onlysomepatients and patient not in somepatients:
        continue
    index = 0
    data = {}
    chrlines = []
    prevchr = '1'
    for line in open(VAFdir + file, "r"):
        lvec = line.rstrip().split("\t")
        if "Patient" in line:
            labels = lvec[5:]
            for label in labels:
                data[label] = [[], []]
            continue
        for n in range(5, len(lvec)):
            if lvec[n] != "":
                data[labels[n-5]][0].append(index)
                data[labels[n-5]][1].append(float(lvec[n]))
        if lvec[2] != prevchr:
            chrlines.append(index)
            prevchr = lvec[2]
        index += 1
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
    #colors = ['#000000', '#1f77b4', '#ff7f0e', '#d62728', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#8c564b', '#0000ff', '#00ff00', '#ff00ff', '#00ffff', '#ff00ff', '#ffff00', '#5500ff', '#ff5500', '#00ff55', '#88ff00', '#ff0088', '#0088ff']
    kelly_colors = ['#000000', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26', '#0000ff', '#00ff00', '#ff0000']
    labels = sortLabels(labels)
    n=0
    for label in labels:
        if len(data[label][0]) < 25:
            continue
        axL = ax.scatter(data[label][0], data[label][1], label=label, s=4, c=kelly_colors[n])
        n += 1
    for index in chrlines:
        axBreak = ax.plot([index-1, index], [0, 1], marker="", color="black")
    plt.ylim(bottom=0.0, top=1.0)
    plt.xlabel(patient + ", " + sample)
    ax.legend()
    plt.savefig(outdir + patient + "_" + sample + "_VAFclusters.png")
    plt.show()


