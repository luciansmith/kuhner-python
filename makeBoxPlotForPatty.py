#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:32:32 2019

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
import seaborn as sns

import numpy
import math
import matplotlib.pyplot as plt
import csv

import lucianSNPLibrary as lsl

onlysomepatients = False
somepatients = ["160"]

drivers_file = "20190903_NumEAdriversPerPatient.txt"

NP = []
prog = []

for line in open(drivers_file, 'r'):
    if "ID" in line:
        continue
    lvec = line.rstrip().split("\t")
    if lvec[0]=="":
        continue
    if not lvec[1] =="":
        NP.append(int(lvec[1]))
    elif not lvec[2] =="":
        prog.append(int(lvec[2]))

all_data = [NP, prog]
subplot = plt.subplot()
plot1 = subplot.boxplot(all_data, patch_artist="true", widths=(0.5, 0.5))
colors = ["blue", "orangered"]
for patch, color in zip(plot1['boxes'], colors):
    patch.set_facecolor(color)
    
for line in plot1['medians']:
    line.set_color("black")

for flier in plot1['fliers']:
    flier.set_marker("*")
    flier.set_markersize(3)

subplot.set_xticklabels(("NP", "Prog"))

plt.ylim(top=20)
plt.savefig("NP_Prog_EAdrivers_boxplot.png")
plt.show()