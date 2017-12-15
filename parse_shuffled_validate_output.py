#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 13:22:09 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile
import sys

import lucianSNPLibrary as lsl
import numpy

outdir = "shuffled_output/"
corefilename = "validate_expands_full"


sys.stdout.write("Filename\tCan't Fail\tCan't Pass\tPass Perfectly\tPass Nigh-Perfectly\tPass As Perfect As Possible\tPass Well\tPass OK\tFail\tFailures\n")
flist = []
for (_, _, f) in walk(outdir):
    flist += f
    
for f in flist:
    if f.find(corefilename) == -1:
        continue
    isPatient = True
    if f.find("persample") != -1:
        isPatient = False
    #(__, __, __, __, __, id, __) = f.split("_")
    freqs = [0, 0, 0, 0, 0, 0, 0, 0, []]
    vfile = open(outdir + f, "r")
    sys.stdout.write(f + "\t")
    for line in vfile:
        if line.find("Patient") != -1:
            continue
        if line.find("---") != -1:
            (patient, nsamples, nsegs, __, __, __) = line.split()
            continue
        (patient, sample, npops, nsegs, nfourg, nthreeg, rfreq, avgfourg, stdfourg) = line.split()
        #Order goes:
            # 0 Can't fail (nfourg==0; stdfourg==0)
            # 1 Can't pass (nfourg!=0; stdfourg==0)
            # 2 Pass perfectly (nfourg==0; rfreq == 0)
            # 3 Pass nigh-perfectly (nfourg!==0; rfreq == 0)
            # 4 Pass as well as possible (nfourg == 0; rfreq != 0)
            # 5 Pass well (0 < rfreq <= .1)
            # 6 Pass OK   (.1 < rfreq < .5)
            # 7 Do not pass (rfreq >= .5)
            # 8 List failures (patient_sample)
        nfourg = int(nfourg)
        nthreeg = int(nthreeg)
        rfreq = float(rfreq)
        avgfourg = float(avgfourg)
        stdfourg = float(stdfourg)
        if (stdfourg == 0):
            if nfourg == 0:
                freqs[0] += 1
            else:
                freqs[1] += 1
        elif rfreq == 0:
            if nfourg == 0:
                freqs[2] += 1
            else:
                freqs[3] += 1
        elif nfourg == 0:
            freqs[4] += 1
        elif rfreq <= 0.1:
            freqs[5] += 1
        elif rfreq < 0.5:
            freqs[6] += 1
        else:
            freqs[7] += 1
            freqs[8].append(patient + "_" + sample)
    for n in range(8):
        sys.stdout.write(str(freqs[n]) + "\t")
    for el in freqs[8]:
        sys.stdout.write(el + ";")
    sys.stdout.write("\n")
    sys.stdout.flush()
        