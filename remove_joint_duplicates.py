#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:14:16 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

tag = "_g500_better_ploidy/"

tree_input = "joint_processed_segmentation" + tag
joint_output = "joint_segmentation_nodupes" + tag

if not(path.isdir(joint_output + "/")):
    mkdir(joint_output + "/")

segment_files = []
for (_, _, f) in walk(tree_input):
    segment_files += f

for f in segment_files:
    if f.find("characters") == -1:
        continue
    patient = f.split("_")[0]
    print("Processing", patient)
    segfile = open(tree_input + f, "r")
    joint_out = open(joint_output + f, "w")
    prev = segfile.readline().split()
    for line in segfile:
        linevec = line.split()
        if (linevec[3:] == prev[3:]) and (linevec[0] == prev[0]):
            prev[2] = linevec[2]
        else:
            joint_out.write(prev[0])
            for n in range(1,len(prev)):
                joint_out.write("\t" + str(prev[n]))
            joint_out.write("\n")
            prev = linevec
    joint_out.write(prev[0])
    for n in range(1,len(prev)):
        joint_out.write("\t" + str(prev[n]))
    joint_out.write("\n")
    joint_out.close()

