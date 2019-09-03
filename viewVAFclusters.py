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

onlysomepatients = False
somepatients = ["160"]

VAFdir = "VAFclusters/"
outdir = "VAFclusters_pngs/"

if not path.isdir(outdir):
    mkdir(outdir)

VAFfiles = []
for __, _, files in walk(VAFdir):
    VAFfiles += files

for file in VAFfiles:
    if "_VAFs" not in file:
        continue
    (patient, sample) = file.split("_")[0:2]
    if onlysomepatients and patient not in somepatients:
        continue
    index = 0
    data = {}
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
        index += 1
    fig = plt.figure(figsize=(30, 10))
    ax = fig.add_subplot(111)
    colors = ['#1f77b4', '#ff7f0e', '#d62728', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#8c564b']
    for n, label in enumerate(labels):
        axL = ax.scatter(data[label][0], data[label][1], label=label, s=1, c=colors[n])
    plt.ylim(bottom=0.0, top=1.0)
    plt.xlabel(patient + ", " + sample)
    ax.legend()
    plt.savefig(outdir + patient + "_" + sample + "_VAFclusters.png")
    plt.show()


