#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:50:09 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
from os.path import isfile

import lucianSNPLibrary as lsl

pASCAT_root = "gamma_test_output/pASCAT_input_g"
bestdir = "best_analyses/"
outdir = "nonintegerCNs/"

use500 = True

onlysomepatients = False
somepatients = ["512"]

firstpatients = ["17", "42", "55", "59", "74", "43", "184", "163", "396", "1047"]

if not path.isdir(outdir):
    mkdir(outdir)

missing = open(outdir + "missing_samples.tsv", "w")

(labels, rev_labels) = lsl.getSNPLabelsAll(False)

files = []
for (__, __, f) in walk(bestdir):
    files += f
for f in files:
    if "catch" in f:
        continue
    patient = f.split("_")[0]
    if onlysomepatients and patient not in somepatients:
        continue
#    if patient in firstpatients:
#        #These need to be re-run, since the SNPs were off.
#        continue
    if use500:
        canonical_ascats = lsl.getSingleGammaCallsFor(patient, "500")
    else:
        canonical_ascats = lsl.getCanonicalAscatCallsFor(patient)
    for canon in canonical_ascats:
        (sample, gamma, constraint) = canon
        if gamma=="None":
            continue
        if constraint=="eight":
            print("Found different 'eight' solution for", patient, sample, gamma)
        rawsegs_dir = pASCAT_root + gamma + "/" + constraint + "/" 
        rawsegs_file = patient + "_" + sample + "_raw_segments.txt"
        outname = outdir + patient + "_" + sample + "_g" + gamma + "_" + constraint + "_nonint_CNs.txt"
        
        outfile = open(outname,"w")
        outfile.write("patient\tbiopsy\tchrom\tsegstart\tsegend\trawA\trawB\tintA\tintB\n")

        #print("Analyzing", patient, sample, constraint)
        if not (lsl.collatepASCATOutput(rawsegs_dir, rawsegs_file, outfile, labels)):
            missing.write(patient)
            missing.write("\t" + sample)
            missing.write("\t" + gamma)
            missing.write("\t" + constraint)
            missing.write("\n")
        
        outfile.close()

missing.close()