# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:22:25 2018

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

oldprobes = {}
oldprobesetfile = open("probe_set_build37_forpartek.txt", "r")
for line in oldprobesetfile:
    if line.find("Chr") != -1:
        continue
    (name, chr, pos) = line.split()
    oldprobes[name] = (chr, pos)
    
newprobes = {}
newprobesetfile = open("HumanOmni1-Quad_v1-0_H_SNPlist.txt", "r")
for line in newprobesetfile:
    if line.find("Chr") != -1:
        continue
    (name, __, chr, pos) = line.split()
    newprobes[name] = (chr, pos)

same = 0
missing = 0
missing_from_placed = 0
missing_from_0 = 0
changed = 0
changed_cnvi = 0
changed_0_to_0 = 0
changed_0_to_placed = 0
changed_placed_to_0 = 0
changed_placed_to_placed = 0

new_chromosome = 0
new_pos_1_1000 = 0
new_pos_1000_10000 = 0
new_pos_10000_100000 = 0
new_pos_100000_1000000 = 0
new_pos_1000000_plus = 0


added = 0
for name in oldprobes:
    if name in newprobes:
        if newprobes[name] == oldprobes[name]:
            same += 1
        else:
            changed += 1
            if name.find("cnvi") != -1:
                changed_cnvi += 1
            else:
                if oldprobes[name][0] == "0" or oldprobes[name][1] == "0":
                    if newprobes[name][0] == "0" or newprobes[name][1] == "0":
                        changed_0_to_0 += 1
                    else:
                        changed_0_to_placed += 1
                        print(name, "was unknown at", oldprobes[name], "; is now placed at", newprobes[name])
                else:
                    if newprobes[name][0] == "0" or newprobes[name][1] == "0":
                        changed_placed_to_0 += 1
                    else:
                        changed_placed_to_placed += 1
                        if newprobes[name][0] == oldprobes[name][0]:
                            oldpos = int(oldprobes[name][1])
                            newpos = int(newprobes[name][1])
                            diff = abs(oldpos-newpos)
                            if diff<1000:
                                new_pos_1_1000 += 1
                            elif diff<10000:
                                new_pos_1000_10000 += 1
                            elif diff<100000:
                                new_pos_10000_100000 += 1
                            elif diff<1000000:
                                new_pos_100000_1000000 += 1
                            else:
                                new_pos_1000000_plus += 1
                        else:
                            new_chromosome += 1
                
            #print(name, "changed from", oldprobes[name], "to", newprobes[name])
    else:
        missing += 1
        if oldprobes[name][0] == "0" or oldprobes[name][1] == "0":
            missing_from_0 += 1
        else:
            missing_from_placed += 1
        #print(name, "changed from", oldprobes[name], "to not being present in the new set.")

for name in newprobes:
    if name in oldprobes:
        continue
    added += 1
    print(name, "wasn't present before, and is now mapped to", newprobes[name])

print("Same:", same)
print("Missing, overall:", missing)
print("Missing; was placed:", missing_from_placed)
print("Missing; was 0:", missing_from_0)
print("Changed, overall:", changed)
print("Changed from 0 to 0", changed_0_to_0)
print("Changed from 0 to placed", changed_0_to_placed)
print("Changed from placed to 0", changed_placed_to_0)
print("Changed from placed to placed", changed_placed_to_placed)
print("Added:", added)

print("Moved to a new chromosome:", new_chromosome)
print("Moved within 1000 bases:", new_pos_1_1000)
print("Moved from 1000 to 10000 bases:", new_pos_1000_10000)
print("Moved from 10000 to 100000 bases:", new_pos_10000_100000)
print("Moved from 100000 to 1000000 bases:", new_pos_100000_1000000)
print("Moved greater than 1000000 bases:", new_pos_1000000_plus)

print("Changed, overall, cnvi", changed_cnvi)
