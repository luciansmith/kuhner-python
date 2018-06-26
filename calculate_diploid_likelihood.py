
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

use_challenge = True

if use_challenge:
    evidence = "calling_evidence_challenge_inc.tsv"
    oddsfile = "calling_evidence_challenge_inc_odds.tsv"
else:
    evidence = "calling_evidence.tsv"
    oddsfile = "calling_evidence_odds.tsv"



evidence_keys = ["flow", "tom", "VAF", "close_diploid", "close_tetraploid", "better_acc"]
def calculateBayes(orig, d_h, d_notH):
    return (orig*d_h) / (d_h*orig + d_notH*(1-orig))
    # P(H|D) = [P(D|H)P(H)]/[P(D|H)P(H)+P(D|H')(1-P(H))]


def writeHeader(outfile):
    outfile.write("Patient")
    outfile.write("\tSample")
    outfile.write("\tFlow ratio")
    outfile.write("\tAdjusted Flow ratio")
    outfile.write("\tTom's Partek Column baseline 2n")
    outfile.write("\tTom's Partek: odds")
    outfile.write("\tVAF histogram category")
    outfile.write("\tVAF histogram odds")
    outfile.write("\tClose diploid flow?")
    outfile.write("\tClose diploid odds")
    outfile.write("\tClose tetraploid flow?")
    outfile.write("\tClose tetraploid odds")
    outfile.write("\tBetter accuracy?")
    outfile.write("\tBetter accuracy odds")
    outfile.write("\tPrior from flow")
    outfile.write("\tOdds of being diploid")
    outfile.write("\n")
    

def writeNewLine(outfile, patient, sample, evidence, oddsvec, odds_str):
    outfile.write(patient)
    outfile.write("\t" + sample)
    for key in evidence_keys:
        outfile.write("\t" + evidence[key])
        outfile.write("\t" + odds_str[key])
    for odds in oddsvec:
        outfile.write("\t" + str(odds))
    outfile.write("\n")


f = open(evidence, "r")
outfile = open(oddsfile, "w")
writeHeader(outfile)


unused_VAFs = set()

for line in f:
    evidence = {}
    odds_str = {}
    if "Patient" in line:
        continue
    (patient, sample, tom_call, VAF_type, flow_ratio, close_dip, close_tet, better_accuracy) = line.rstrip().split('\t')
    evidence["tom"] = tom_call
    evidence["VAF"] = VAF_type
    evidence["close_diploid"] = close_dip
    evidence["close_tetraploid"] = close_tet
    evidence["better_acc"] = better_accuracy
    (flowa, flowb) = flow_ratio.split('::')
    flowa = int(flowa)
    flowb = int(flowb)
    evidence["flow"] = str(flowb-flowa) + "::" + str(flowa)
    base_odds = 0.5
    if flowb>0:
        base_odds = flowa/(flowb)
    odds_str["flow"] = flow_ratio
    if (base_odds < 0.5):
        base_odds_skew = (flowa+3)/(flowb+3)
        odds_str["flow"] =  str(flowb-flowa) + "::" + str(flowa+3)
        if base_odds_skew > 0.5:
            base_odds = 0.5
            odds_str["flow"] = str(int(flowb/2)) + "::" + str(int(flowb/2))
        else:
            base_odds = base_odds_skew
    else:
        base_odds_skew = flowa/(flowb+3)
        odds_str["flow"] = str(flowb+3-flowa) + "::" + str(flowa)
        if base_odds_skew < 0.5:
            base_odds = 0.5
            odds_str["flow"] = str(int(flowb/2)) + "::" + str(int(flowb/2))
        else:
            base_odds = base_odds_skew
    
    oddsvec = [1-base_odds, 1-base_odds]

    ptom_dip = 0.50
    ptom_tet = 0.50
    if tom_call == "yes" or tom_call=="gastric":
        ptom_dip = 0.95
        ptom_tet = 0.30
    elif tom_call =="no":
        ptom_dip = 0.05
        ptom_tet = 0.70
    elif tom_call == "?":
        #Don't do anything
        assert(True)
    elif tom_call == "Unknown":
        #Don't do anything
        assert(True)
    else:
        print("Unknown tom call:", tom_call)
    oddsvec[1] = calculateBayes(oddsvec[1], ptom_dip, ptom_tet)
    odds_str["tom"] = str(int(100*ptom_dip)) + "::" + str(int(100*ptom_tet))

    pVAF_dip = 0.5
    pVAF_tet = 0.5
    if VAF_type == "0.1 and 0.5 clean" or VAF_type=="0.1 and 0.4 anomaly":
        pVAF_dip = 0.95
        pVAF_tet = 0.30
    elif VAF_type == "Clear .25 peak":
        pVAF_dip = 0.2
        pVAF_tet = 0.8
    elif VAF_type == "Primarily 0.5":
        pVAF_dip = 0.85
        pVAF_tet = 0.40
    else:
        unused_VAFs.add(VAF_type)
    oddsvec[1] = calculateBayes(oddsvec[1], pVAF_dip, pVAF_tet) #P(01/05|dip, P(01/05|tet))
    odds_str["VAF"] = str(int(100*pVAF_dip)) + "::" + str(int(100*pVAF_tet))

    pclose_d_dip = 0.5
    pclose_d_tet = 0.5
    if close_dip == "True":
        pclose_d_dip = 0.9
        pclose_d_tet = 0.86
    elif close_dip == "False":
        pclose_d_dip = 0.1
        pclose_d_tet = 0.14
    elif close_dip == "Two":
        #This case may well not be informative at all
        pclose_d_dip = 0.5
        pclose_d_tet = 0.5
    oddsvec[1] = calculateBayes(oddsvec[1], pclose_d_dip, pclose_d_tet) #P(close_dip_T|dip, P(close_dip_T|tet))
    odds_str["close_diploid"] = str(int(100*pclose_d_dip)) + "::" + str(int(100*pclose_d_tet))

    pclose_t_dip = 0.5
    pclose_t_tet = 0.5
    if close_tet == "True":
        pclose_t_dip = 0.5
        pclose_t_tet = 0.86
        oddsvec[1] = calculateBayes(oddsvec[1], 0.5, 0.86) #P(close_tet_T|dip, P(close_tet_T|tet))
    elif close_tet == "False":
        pclose_t_dip = 0.5
        pclose_t_tet = 0.14
    oddsvec[1] = calculateBayes(oddsvec[1], pclose_t_dip, pclose_t_tet) #P(close_tet_F|dip, P(close_tet_F|tet))
    odds_str["close_tetraploid"] = str(int(100*pclose_t_dip)) + "::" + str(int(100*pclose_t_tet))


    pbetter_acc_dip = 0.5
    pbetter_acc_tet = 0.5
    if better_accuracy == "Diploid":
        pbetter_acc_dip = 0.7
        pbetter_acc_tet = 0.2
    elif better_accuracy == "Tetraploid":
        pbetter_acc_dip = 0.2
        pbetter_acc_tet = 0.7
    if better_accuracy == "Diploid only":
        pbetter_acc_dip = 0.09
        pbetter_acc_tet = 0.01
    elif better_accuracy == "Tetraploid only":
        pbetter_acc_dip = 0.01
        pbetter_acc_tet = 0.09
    oddsvec[1] = calculateBayes(oddsvec[1], pbetter_acc_dip, pbetter_acc_tet) #P(dip_better|dip, P(dip_better|tet))
    odds_str["better_acc"] = str(int(100*pbetter_acc_dip)) + "::" + str(int(100*pbetter_acc_tet))

#    if (has_tetraploid == "True"):
#        oddsvec[1] = calculateBayes(oddsvec[1], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
#        oddsvec[2] = calculateBayes(oddsvec[2], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
#        oddsvec[3] = calculateBayes(oddsvec[3], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
#        oddsvec[4] = calculateBayes(oddsvec[4], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)

    writeNewLine(outfile, patient, sample, evidence, oddsvec, odds_str)

for unused in unused_VAFs:
    print("Unused VAF call", unused)

outfile.close()