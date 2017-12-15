
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

tag = "25M_v2"

evidence = "joint_summary_" + tag + ".txt"
oddsfile = "joint_summary_" + tag + "_odds.txt"

CN_dir = "CN_joint_log2rs_" + tag + "_"
BAF_dir = "BAF_joint_vals_" + tag + "_"

un = "unconstrained"
diploid = "diploid"
tetraploid = "tetraploid"
out = "joined_best"

if not(path.isdir(CN_dir + out + "/")):
    mkdir(CN_dir + out + "/")

if not(path.isdir(BAF_dir + out + "/")):
    mkdir(BAF_dir + out + "/")

def calculateBayes(orig, d_h, d_notH):
    return (orig*d_h) / (d_h*orig + d_notH*(1-orig))
    # P(H|D) = [P(D|H)P(H)]/[P(D|H)P(H)+P(D|H')(1-P(H))]

def writeNewLine(f, line, odds):
    f.write(line.rstrip())
    for odd in odds:
        f.write("\t" + str(odd))
    f.write("\n")

def moveFiles(patient, sample, which, second):
    #return
    fnames = glob.glob(CN_dir + which + "/" + patient + "_" + sample + "*")
    if (len(fnames)==0):
        print "Can't find any file that matches", patient+"_"+sample, "in", CN_dir + which
        return
    elif(len(fnames) > 1):
        print "Too many files match: ", fnames
        return
    #print fnames
    shutil.copyfile(fnames[0], CN_dir+out+"/"+patient+"_"+sample+"_avglog2rs.txt")
    fnames = glob.glob(BAF_dir + which + "/" + patient + "_" + sample + "*")
    if (len(fnames)==0):
        print "Can't find any file that matches", patient+"_"+sample, "in", BAF_dir + which
        return
    elif(len(fnames) > 1):
        print "Too many files match: ", fnames
        return
    #print fnames
    shutil.copyfile(fnames[0], BAF_dir+out+"/"+patient+"_"+sample+"_avgbafvs.txt")

    if (second != ""):
        fnames = glob.glob(CN_dir + second + "/" + patient + "_" + sample + "*")
        if (len(fnames)==0):
            print "Can't find any file that matches", patient+"_"+sample, "in", CN_dir + second
            return
        elif(len(fnames) > 1):
            print "Too many files match: ", fnames
            return
        #print fnames
        shutil.copyfile(fnames[0], CN_dir+out+"/"+patient+"_"+sample+"b_avglog2rs.txt")

        fnames = glob.glob(BAF_dir + second + "/" + patient + "_" + sample + "*")
        if (len(fnames)==0):
            print "Can't find any file that matches", patient+"_"+sample, "in", BAF_dir + second
            return
        elif(len(fnames) > 1):
            print "Too many files match: ", fnames
            return
        #print fnames
        shutil.copyfile(fnames[0], BAF_dir+out+"/"+patient+"_"+sample+"b_avgbafvs.txt")

f = open(evidence, "r")
outfile = open(oddsfile, "w")

for line in f:
    if (line.find("patient") != -1):
        writeNewLine(outfile, line, ["Odds","odds no unmatch","odds no four","odds no best","odds has ploidy"])
        continue
    (patient, sample, un_ploidy, un_purity, un_psi, un_goodness, diploid_ploidy, diploid_purity, diploid_psi, diploid_goodness, tetraploid_ploidy, tetraploid_purity, tetraploid_psi, tetraploid_goodness, un_matches, exactly_four, has_diploid, has_tetraploid, dip_to_nondip_flows, goodness_diff) = line.translate(None, '"').split()

    if (diploid_purity == "NA"):
#        writeNewLine(outfile, line, ["auto-unconstrained"])
#        moveFiles(patient, sample, un, "")
        if (tetraploid_purity == "NA"):
            writeNewLine(outfile, line, ["auto-unconstrained"])
            moveFiles(patient, sample, un, "")
        else:
            writeNewLine(outfile, line, ["auto-tetraploid"])
            moveFiles(patient, sample, tetraploid, "")
        continue
    elif (tetraploid_purity == "NA"):
#        writeNewLine(outfile, line, ["auto-unconstrained"])
#        moveFiles(patient, sample, un, "")
        writeNewLine(outfile, line, ["auto-diploid"])
        moveFiles(patient, sample, diploid, "")
        continue

    odds = 0.5 # Should be the relative abundance of diploid vs. tetraploid entries (the prior)
    #Weigh the odds by the flow data overall, skewed towards 50% by 3 'missed' reads.
    (num_dip, num_nondip) = dip_to_nondip_flows.split("::")
    num_dip = int(num_dip)
    num_nondip = int(num_nondip)
    if (num_dip > num_nondip + 3):
        odds = num_dip / (num_dip + num_nondip + 3)
    elif (num_dip + 3 < num_nondip):
        odds = (num_dip + 3) / (num_dip + num_nondip)
    oddsvec = [odds, odds, odds, odds, odds]

    #Whether the unconstrained ASCAT gives the diploid or the tetraploid results:
    if un_matches=="diploid":
        oddsvec[0] = calculateBayes(oddsvec[0], 0.95, 0.5) #P(dip_match|dip, P(dip_match|tet))
        oddsvec[2] = calculateBayes(oddsvec[2], 0.95, 0.5) #P(dip_match|dip, P(dip_match|tet))
        oddsvec[3] = calculateBayes(oddsvec[3], 0.95, 0.5) #P(dip_match|dip, P(dip_match|tet))
        oddsvec[4] = calculateBayes(oddsvec[4], 0.95, 0.5) #P(dip_match|dip, P(dip_match|tet))
    elif un_matches=="tetraploid":
        oddsvec[0] = calculateBayes(oddsvec[0], 0.05, 0.5) #P(tet_match|dip, P(tet_match|tet))
        oddsvec[2] = calculateBayes(oddsvec[2], 0.05, 0.5) #P(tet_match|dip, P(tet_match|tet))
        oddsvec[3] = calculateBayes(oddsvec[3], 0.05, 0.5) #P(tet_match|dip, P(tet_match|tet))
        oddsvec[4] = calculateBayes(oddsvec[4], 0.05, 0.5) #P(tet_match|dip, P(tet_match|tet))

    #Whether the tetraploid results are within 0.1 of 4.0:
    if (exactly_four == "diploid"):
        oddsvec[0] = calculateBayes(oddsvec[0], 0.75, 0.2) #P(tet_is_four|dip), P(tet_is_four|tet)
        oddsvec[1] = calculateBayes(oddsvec[1], 0.75, 0.2) #P(tet_is_four|dip), P(tet_is_four|tet)
        oddsvec[3] = calculateBayes(oddsvec[3], 0.75, 0.2) #P(tet_is_four|dip), P(tet_is_four|tet)
        oddsvec[4] = calculateBayes(oddsvec[4], 0.75, 0.2) #P(tet_is_four|dip), P(tet_is_four|tet)
    elif (exactly_four == "tetraploid"):
        oddsvec[0] = calculateBayes(oddsvec[0], 0.25, 0.8) #P(tet_not_four|dip), p(tet_not_four|tet)
        oddsvec[1] = calculateBayes(oddsvec[1], 0.25, 0.8) #P(tet_not_four|dip), p(tet_not_four|tet)
        oddsvec[3] = calculateBayes(oddsvec[3], 0.25, 0.8) #P(tet_not_four|dip), p(tet_not_four|tet)
        oddsvec[4] = calculateBayes(oddsvec[4], 0.25, 0.8) #P(tet_not_four|dip), p(tet_not_four|tet)

    #Whether the better purity is diploid or not
#    if (best_purity == "diploid"):
#        oddsvec[0] = calculateBayes(oddsvec[0], 0.75, 0.75) #P(dip_best|dip), P(dip_best|tet)
#        oddsvec[1] = calculateBayes(oddsvec[1], 0.75, 0.75) #P(dip_best|dip), P(dip_best|tet)
#        oddsvec[2] = calculateBayes(oddsvec[2], 0.75, 0.75) #P(dip_best|dip), P(dip_best|tet)
#        oddsvec[4] = calculateBayes(oddsvec[4], 0.75, 0.75) #P(dip_best|dip), P(dip_best|tet)
#    elif (best_purity == "tetraploid"):
#        oddsvec[0] = calculateBayes(oddsvec[0], 0.25, 0.25) #P(tet_best|dip), P(tet_best|tet)
#        oddsvec[1] = calculateBayes(oddsvec[1], 0.25, 0.25) #P(tet_best|dip), P(tet_best|tet)
#        oddsvec[2] = calculateBayes(oddsvec[2], 0.25, 0.25) #P(tet_best|dip), P(tet_best|tet)
#        oddsvec[4] = calculateBayes(oddsvec[4], 0.25, 0.25) #P(tet_best|dip), P(tet_best|tet)

    #Whether the diploid ploidy is represented in the flow data
    if (has_diploid == "True"):
        oddsvec[0] = calculateBayes(oddsvec[0], 0.99, 0.8) #P(has_dip|dip), P(has_dip|tet)
        oddsvec[1] = calculateBayes(oddsvec[1], 0.99, 0.8) #P(has_dip|dip), P(has_dip|tet)
        oddsvec[2] = calculateBayes(oddsvec[2], 0.99, 0.8) #P(has_dip|dip), P(has_dip|tet)
        oddsvec[3] = calculateBayes(oddsvec[3], 0.99, 0.8) #P(has_dip|dip), P(has_dip|tet)

    #Whether the tetraploid ploidy is represented in the flow data
    if (has_tetraploid == "True"):
        oddsvec[0] = calculateBayes(oddsvec[0], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
        oddsvec[1] = calculateBayes(oddsvec[1], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
        oddsvec[2] = calculateBayes(oddsvec[2], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)
        oddsvec[3] = calculateBayes(oddsvec[3], 0.01, 0.2) #P(has_tet|dip), P(has_tet|tet)

    writeNewLine(outfile, line, oddsvec)
    second = ""
    if (odds >= 0.5):
        if (odds < 0.75):
            second = tetraploid
        moveFiles(patient, sample, diploid, second)
    else:
        if odds>0.25:
            second = diploid
        moveFiles(patient, sample, tetraploid, second)

outfile.close()