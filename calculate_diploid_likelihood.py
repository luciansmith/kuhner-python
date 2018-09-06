
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
only_challenge_evidence = False

if use_challenge:
    evidence = "calling_evidence_challenge_inc.tsv"
    oddsfile = "calling_evidence_challenge_inc_odds.tsv"
elif only_challenge_evidence:
    evidence = "calling_evidence.tsv"
    oddsfile = "calling_evidence_only_challev_odds.tsv"
else:
    evidence = "calling_evidence.tsv"
    oddsfile = "calling_evidence_odds.tsv"



evidence_keys = ["flow", "tom", "VAF", "close_diploid", "close_tetraploid", "better_acc", "nygc", "goodness", "xcloser"]
if only_challenge_evidence:
    evidence_keys = ["flow", "close_diploid", "close_tetraploid", "better_acc", "goodness", "xcloser"]

def calculateBayes(orig, d_h, d_notH):
    return (orig*d_h) / (d_h*orig + d_notH*(1-orig))
    # P(H|D) = [P(D|H)P(H)]/[P(D|H)P(H)+P(D|H')(1-P(H))]


def writeHeader(outfile):
    outfile.write("Patient")
    outfile.write("\tSample")
    outfile.write("\tFlow ratio")
    outfile.write("\tAdjusted Flow ratio")
    if not only_challenge_evidence:
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
    if not only_challenge_evidence:
        outfile.write("\tNYGC closer?")
        outfile.write("\tNYGC closer odds")
    outfile.write("\tBetter goodness")
    outfile.write("\tBetter goodness odds")
    outfile.write("\tXiaohong closer")
    outfile.write("\tXiaohong closer odds")
    outfile.write("\tPrior from flow")
    outfile.write("\tOdds of being diploid")
    outfile.write("\tFinal call from Hutch")
    outfile.write("\n")
    

def writeNewLine(outfile, patient, sample, evidence, oddsvec, odds_str, final):
    outfile.write(patient)
    outfile.write("\t" + sample)
    for key in evidence_keys:
        outfile.write("\t" + evidence[key])
        outfile.write("\t" + odds_str[key])
    for odds in oddsvec:
        outfile.write("\t" + str(odds))
    outfile.write("\t" + final)
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
    (patient, sample, tom_call, VAF_type, flow_ratio, close_dip, close_tet, better_accuracy, __, nygc, goodness, xcloser, final) = line.rstrip().split('\t')
    evidence["tom"] = tom_call
    evidence["VAF"] = VAF_type
    evidence["close_diploid"] = close_dip
    evidence["close_tetraploid"] = close_tet
    evidence["better_acc"] = better_accuracy
    evidence["nygc"] = nygc
    evidence["goodness"] = goodness
    evidence["xcloser"] = xcloser
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

    if not only_challenge_evidence:
        ptom_dip = 0.50
        ptom_tet = 0.50
        if tom_call == "yes" or tom_call=="gastric":
            ptom_dip = 0.970
            ptom_tet = 0.054
        elif tom_call =="no":
            ptom_dip = 0.003
            ptom_tet = 0.629
        elif tom_call == "?":
            ptom_dip = 0.033
            ptom_tet = 0.314
        elif tom_call == "Unknown":
            #Don't do anything
            assert(True)
        else:
            print("Unknown tom call:", tom_call)
        oddsvec[1] = calculateBayes(oddsvec[1], ptom_dip, ptom_tet)
        odds_str["tom"] = str(int(100*ptom_dip)) + "::" + str(int(100*ptom_tet))

        pVAF_dip = 0.5
        pVAF_tet = 0.5
        if VAF_type == "0.1 and 0.4 anomaly":
            pVAF_dip = 0.125
            pVAF_tet = 0.054
        elif VAF_type == "0.1 and 0.5 clean":
            pVAF_dip = 0.310
            pVAF_tet = 0.054
        elif VAF_type == "Clear .25 peak":
            pVAF_dip = 0.043
            pVAF_tet = 0.200
        elif VAF_type == "Multi-hump":
            pVAF_dip = 0.243
            pVAF_tet = 0.429
        elif VAF_type == "Primarily 0.5":
            pVAF_dip = 0.079
            pVAF_tet = 0.027
        elif VAF_type == "Single peak anomaly":
            pVAF_dip = 0.118
            pVAF_tet = 0.200
        elif VAF_type == "Single sharp peak low VAF":
            pVAF_dip = 0.040
            pVAF_tet = 0.027
        else:
            unused_VAFs.add(VAF_type)
        oddsvec[1] = calculateBayes(oddsvec[1], pVAF_dip, pVAF_tet) #P(01/05|dip, P(01/05|tet))
        odds_str["VAF"] = str(int(100*pVAF_dip)) + "::" + str(int(100*pVAF_tet))

    pclose_d_dip = 0.5
    pclose_d_tet = 0.5
    #We are lumping 'true' and 'false' together here, since there were only 13 samples that weren't 'two', and both sets
    # indicated 'this is probably tetraploid'.
    if close_dip == "True":
        pclose_d_dip = 0.007
        pclose_d_tet = 0.091
    elif close_dip == "False":
        pclose_d_dip = 0.007
        pclose_d_tet = 0.091
    elif close_dip == "Two":
        pclose_d_dip = 0.997
        pclose_d_tet = 0.229
    oddsvec[1] = calculateBayes(oddsvec[1], pclose_d_dip, pclose_d_tet) #P(close_dip_T|dip, P(close_dip_T|tet))
    odds_str["close_diploid"] = str(int(100*pclose_d_dip)) + "::" + str(int(100*pclose_d_tet))

    pclose_t_dip = 0.5
    pclose_t_tet = 0.5
    if close_tet == "True":
        pclose_t_dip = 0.393
        pclose_t_tet = 0.629
        oddsvec[1] = calculateBayes(oddsvec[1], 0.5, 0.86) #P(close_tet_T|dip, P(close_tet_T|tet))
    elif close_tet == "False":
        pclose_t_dip = 0.601
        pclose_t_tet = 0.378
    oddsvec[1] = calculateBayes(oddsvec[1], pclose_t_dip, pclose_t_tet) #P(close_tet_F|dip, P(close_tet_F|tet))
    odds_str["close_tetraploid"] = str(int(100*pclose_t_dip)) + "::" + str(int(100*pclose_t_tet))


    pbetter_acc_dip = 0.5
    pbetter_acc_tet = 0.5
    if better_accuracy == "Diploid":
        assert(True)
        #This doesn't actually tell you anything (sigh)
    elif better_accuracy == "Tetraploid":
        pbetter_acc_dip = 0.144
        pbetter_acc_tet = 0.171
    if better_accuracy == "Diploid only":
#        if int(sample) < 23341:
#            pbetter_acc_dip = 1
#            pbetter_acc_tet = 0
        assert(True)
        #To few to assess
    elif better_accuracy == "Tetraploid only":
        pbetter_acc_dip = 0.003
        pbetter_acc_tet = 0.629
#        if int(sample) < 23341:
#            pbetter_acc_dip = 0
#            pbetter_acc_tet = 1
        
    oddsvec[1] = calculateBayes(oddsvec[1], pbetter_acc_dip, pbetter_acc_tet) #P(dip_better|dip, P(dip_better|tet))
    odds_str["better_acc"] = str(int(100*pbetter_acc_dip)) + "::" + str(int(100*pbetter_acc_tet))

    if not only_challenge_evidence:
        nygc_dip = 0.5
        nygc_tet = 0.5
        if nygc == "Diploid":
            nygc_dip = 0.986
            nygc_tet = 0.286
        elif nygc == "Tetraploid":
            nygc_dip = 0.017
            nygc_tet = 0.727
        oddsvec[1] = calculateBayes(oddsvec[1], nygc_dip, nygc_tet) #P(dip_better|dip, P(dip_better|tet))
        odds_str["nygc"] = str(int(100*nygc_dip)) + "::" + str(int(100*nygc_tet))

    goodness_dip = 0.5
    goodness_tet = 0.5
    if goodness == "Diploid":
        goodness_dip = 0.884
        goodness_tet = 0.216
    elif goodness == "Tetraploid":
        goodness_dip = 0.111
        goodness_tet = 0.143
    oddsvec[1] = calculateBayes(oddsvec[1], goodness_dip, goodness_tet) #P(dip_better|dip, P(dip_better|tet))
    odds_str["goodness"] = str(int(100*goodness_dip)) + "::" + str(int(100*goodness_tet))

    xcloser_dip = 0.5
    xcloser_tet = 0.5
    if xcloser == "Diploid":
        xcloser_dip = 0.970
        xcloser_tet = 0.054
    elif xcloser == "Tetraploid":
        xcloser_dip = 0.003
        xcloser_tet = 0.429
    elif xcloser == "Neither":
        xcloser_dip = 0.033
        xcloser_tet = 0.514
    oddsvec[1] = calculateBayes(oddsvec[1], xcloser_dip, xcloser_tet) #P(dip_better|dip, P(dip_better|tet))
    odds_str["xcloser"] = str(int(100*xcloser_dip)) + "::" + str(int(100*xcloser_tet))


    writeNewLine(outfile, patient, sample, evidence, oddsvec, odds_str, final)

for unused in unused_VAFs:
    print("Unused VAF call", unused)

outfile.close()