
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:30:20 2017

@author: lpsmith
"""

#This opens the evidence file and makes decisions accordingly.

from __future__ import division
#from os import walk

import shutil
import glob
from os import path
from os import mkdir

#import lucianSNPLibrary as lsl

evidence = "calling_evidence.txt"
oddsfile = "calling_evidence_odds.txt"

def calculateBayes(orig, d_h, d_notH):
    return (orig*d_h) / (d_h*orig + d_notH*(1-orig))
    # P(H|D) = [P(D|H)P(H)]/[P(D|H)P(H)+P(D|H')(1-P(H))]


def writeNewLine(f, line, odds):
    f.write(line.rstrip())
    for odd in odds:
        f.write("\t" + str(odd))
    f.write("\n")


f = open(evidence, "r")
outfile = open(oddsfile, "w")


orig = 0.623
print("Orig:", orig)
orig = calculateBayes(orig, 0.2, 0.2)
print("one step:", orig)
orig = calculateBayes(orig, 0.5, 0.5)
print("one step:", orig)
orig = calculateBayes(orig, 0.8, 0.8)
print("one step:", orig)


unused_VAFs = set()

for line in f:
    if "Patient" in line:
        writeNewLine(outfile, line, ["Prior from flow", "Odds"])
        continue
    (patient, sample, tom_call, VAF_type, flow_ratio, close_dip, close_tet, better_accuracy) = line.split('\t')
    (flowa, flowb) = flow_ratio.split('::')
    flowa = int(flowa)
    flowb = int(flowb)
    base_odds = flowa/(flowb)
    if (base_odds < 0.5):
        base_odds_skew = (flowa+3)/(flowb+3)
        if base_odds_skew > 0.5:
            base_odds = 0.5
        else:
            base_odds = base_odds_skew
    else:
        base_odds_skew = flowa/(flowb+3)
        if base_odds_skew < 0.5:
            base_odds = 0.5
        else:
            base_odds = base_odds_skew
    
    oddsvec = [1-base_odds, 1-base_odds]

    if tom_call == "yes" or tom_call=="gastric":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.95, 0.3) #P(tom_yes|dip, P(tom_yes|tet))
    elif tom_call =="no":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.05, 0.7) #P(tom_no|dip, P(tom_no|tet))
#    elif tom_call == "gastric":
#        oddsvec[1] = calculateBayes(oddsvec[1], 1.0, 0.0) #P(tom_gastric|dip, P(tom_no|tet))
    elif tom_call == "?":
        #Don't do anything
        assert(True)
    else:
        print("Unknown tom call:", tom_call)

    if VAF_type == "0.1 and 0.5 clean" or VAF_type=="0.1 and 0.4 anomaly":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.95, 0.3) #P(01/05|dip, P(01/05|tet))
    elif VAF_type == "Clear .25 peak":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.2, 0.8) #P(.25|dip, P(.25|tet))
    elif VAF_type == "Primarily 0.5":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.85, 0.4) #P(01/05|dip, P(01/05|tet))
    else:
        unused_VAFs.add(VAF_type)

    if close_dip == "True":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.9, 0.86) #P(close_dip_T|dip, P(close_dip_T|tet))
    elif close_dip == "False":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.1, 0.14) #P(close_dip_F|dip, P(close_dip_F|tet))
    elif close_dip == "Two":
        #This case may well not be informative at all
        oddsvec[1] = calculateBayes(oddsvec[1], 0.5, 0.5) #P(dip_is_two|dip, P(dip_is_two|tet))

    if close_tet == "True":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.5, 0.86) #P(close_tet_T|dip, P(close_tet_T|tet))
    elif close_tet == "False":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.5, 0.14) #P(close_tet_F|dip, P(close_tet_F|tet))


    if better_accuracy == "Diploid":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.7, 0.2) #P(dip_better|dip, P(dip_better|tet))
    elif close_tet == "Tetraploid":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.2, 0.7) #P(tet_better|dip, P(tet_better|tet))
    if better_accuracy == "Diploid only":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.09, 0.01) #P(dip_only|dip, P(dip_only|tet))
    elif close_tet == "Tetraploid only":
        oddsvec[1] = calculateBayes(oddsvec[1], 0.01, 0.09) #P(tet_only|dip, P(tet_only|tet))

#    if (has_tetraploid == "True"):
#        oddsvec[1] = calculateBayes(oddsvec[1], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
#        oddsvec[2] = calculateBayes(oddsvec[2], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
#        oddsvec[3] = calculateBayes(oddsvec[3], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
#        oddsvec[4] = calculateBayes(oddsvec[4], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)

    writeNewLine(outfile, line, oddsvec)

for unused in unused_VAFs:
    print("Unused VAF call", unused)

outfile.close()