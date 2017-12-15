#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:22:45 2016

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile

import lucianSNPLibrary as lsl
import random

expands_directory = "expands_full_input/"
alldupes_directory = "expands_full_input_alldupes/"

flist = []
filesets = {}
for (_, _, f) in walk(expands_directory):
    flist += f
    
for f in flist:
    if f.find("BAF") == -1:
        continue
#    print "Workin on", f
    baffile = open(expands_directory + f, "r")
    bafdata = []
    for line in baffile:
        if (line.find("chr") != -1):
            continue
        bafdata.append(line.split())
    
    CNname = f.replace("BAF", "CN")
    CNfile = open(expands_directory + CNname, "r")
    CNdata = []
    for line in CNfile:
        if (line.find("chr") != -1):
            continue
        CNdata.append(line.split())
    
    bafout = open(alldupes_directory + f, "w")
    bafout.write("chr\tstartpos\tendpos\tAF_Tumor\tPN_B\n")
    CNout = open(alldupes_directory + CNname, "w")
    CNout.write("chr\tstartpos\tendpos\tCN_Estimate\n")
    bl = 0
    cl = 0
    while (bl < len(bafdata) and cl < len(CNdata)):
        if bafdata[bl][0] == CNdata[cl][0] and bafdata[bl][1] == CNdata[cl][1] and bafdata[bl][2] == CNdata[cl][2]:
            bafout.write(bafdata[bl][0] + "\t" + bafdata[bl][1] + "\t" + bafdata[bl][2] + "\t" + bafdata[bl][3] + "\t" + bafdata[bl][4] + "\n")
            CNout.write(CNdata[cl][0] + "\t" + CNdata[cl][1] + "\t" + CNdata[cl][2] + "\t" + CNdata[cl][3] + "\n")
            bl += 1
            cl += 1
            continue
        n = 1
        found = False
        while (bl+n < len(bafdata)):
            if bafdata[bl+n][0] == CNdata[cl][0] and bafdata[bl+n][1] == CNdata[cl][1] and bafdata[bl+n][2] == CNdata[cl][2]:
                bl += n
                found = True
                continue
            n += 1
        n = 1
        while (found==False and cl+n < len(CNdata)):
            
            if cl+n < len(CNdata) and bafdata[bl][0] == CNdata[cl+n][0] and bafdata[bl][1] == CNdata[cl+n][1] and bafdata[bl][2] == CNdata[cl+n][2]:
                cl += 1
                found = True
                continue
            n += 1
        if not(found):
            bl += 1
            cl += 1
    
    
    
    
    
    
    
    
    
    
    
    
