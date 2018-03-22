#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:09:06 2018

@author: lpsmith
"""

import os


VAFdir = "VAF_pngs/"

VAFsort = []

for root, dirs, files in os.walk(VAFdir):
    if len(dirs)>0:
        continue
    root = root.replace(VAFdir, "")
    for file in files:
        if "png" in file:
            (patient, sample, __, __) = file.split("_")
            VAFsort.append((int(patient), sample, root))
            
VAFsort.sort()
VAFout = open(VAFdir + "VAF_sort_summary.tsv", "w")
VAFout.write("Patient")
VAFout.write("\tSample")
VAFout.write("\tCategory")
VAFout.write("\n")

for VAF in VAFsort:
    for VAFitem in VAF:
        VAFout.write(str(VAFitem) + "\t")
    VAFout.write("\n")

VAFout.close()