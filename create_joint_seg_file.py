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

onepatientonly = True
onepatient = ["74", "422", "450", "512", "728", "995", "997"]
gamma_dir = "partial_gamma_test_output/gamma_test/pASCAT_input_"

# read the probeset file, which correlates name to position.
labels, rev_labels = lsl.getSNPLabelsAll(False)

#for tag in ["_combined_avSNPs", "_combined_avSNPs_only25M", "_combined_avSNPs_only1M"]:
for tag in ["g1000"]:

    pASCAT_output_dir = gamma_dir + tag + "/"
    joint_out_name = "joint_seg_" + tag

    for constraint in ["unconstrained"]: #["diploid", "tetraploid", "unconstrained"]:
        directory = pASCAT_output_dir + constraint + "/"
        joint_name = joint_out_name + "_" + constraint + ".txt"

        outfile = open(joint_name,"w")
        outfile.write("patient\tbiopsy\tchrom\tsegstart\tsegend\trawA\trawB\tintA\tintB\n")
        filenames = []
        for (_, _, f) in walk(directory):
            filenames += f
            break
        for f in filenames:
            if onepatientonly:
                found = False
                for patient in onepatient:
                    if f.find(patient + "_")==0:
                        found = True
                        continue
                if not found:
                    continue
            if (f.find("raw_segments") == -1):
                continue

            print "Analyzing", f
            lsl.collatepASCATOutput(directory, f, outfile, labels)
        outfile.close()
