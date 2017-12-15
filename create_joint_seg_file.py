#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 14:23:07 2017

@author: lpsmith
"""
from __future__ import division
import lucianSNPLibrary as lsl
from os import walk
#from os import mkdir
#from os import path
#import string

is1M = False
isAveraged = True

onepatientonly = False
onepatient = "951"
    

# read the probeset file, which correlates name to position.
if (isAveraged):
    labels, rev_labels = lsl.getSNPLabelsAveraged(False)
else:
    labels, rev_labels = lsl.getSNPLabels(is1M, False)

#for tag in ["_combined_avSNPs", "_combined_avSNPs_only25M", "_combined_avSNPs_only1M"]:
for tag in ["_25M_v2"]:

    pASCAT_output_dir = "pASCAT_output" + tag + "/"
    joint_out_name = "joint_seg" + tag
    
    for constraint in ["diploid", "tetraploid", "unconstrained"]:
        directory = pASCAT_output_dir + constraint + "/"
        joint_name = joint_out_name + "_" + constraint + ".txt"
    
        outfile = open(joint_name,"w")
        outfile.write("patient\tbiopsy\tchrom\tsegstart\tsegend\trawA\trawB\tintA\tintB\n")
        filenames = []
        for (_, _, f) in walk(directory):
            filenames += f
            break
        for f in filenames:
            if (onepatientonly and f.find(onepatient)==-1):
                continue
            if (f.find("raw_segments") == -1):
                continue
            
            print "Analyzing", f
            lsl.collatepASCATOutput(directory + f, outfile, labels)
        outfile.close()
