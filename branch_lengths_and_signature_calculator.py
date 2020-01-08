# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:43:41 2018

@author: Lucian
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy
from ete3 import Tree, TreeStyle, TextFace

import numpy
import math
import matplotlib.pyplot as plt
import csv

import lucianSNPLibrary as lps
#import imp
#lps = imp.load_source("lps","/home/mkkuhner/Papers/phylo/lucianSNPLibrary.py")

showTrees = True
needDoubleSplitInfo = False

onlysomepatients = True
somepatients = ["126"]
#somepatients = ["387", "483", "609", "611", "626", "631", "635", "652", "852"] #double split patients

#input
mutation_file = "lucian_from_kanika.csv"
treedir = "deconvtrees/"
optimdir = "optim_inputs/"

#output
lengthdir = "branch_lengths/"
sigdir = "branch_signatures/"


if not path.isdir(lengthdir):
    mkdir(lengthdir)

if not path.isdir(sigdir):
    mkdir(sigdir)

def getCNVCall(patient, sample, chrom, pos, CNVs):
    if patient not in CNVs:
        assert(False)
        return (-1, -1)
    if sample not in CNVs[patient]:
        assert(False)
        return (-1, -1)
    if chrom not in CNVs[patient][sample]:
        print("No chromosome", str(chrom), "found.")
        assert(False)
        return (-1, -1)
    for (start, end, call) in CNVs[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return call
    return (-1, -1)

def isDeleted(patient, sample, chrom, pos, deletions):
    if patient not in deletions:
        return False
    if sample not in deletions[patient]:
        return False
    if chrom not in deletions[patient][sample]:
        return False
    for (start, end) in deletions[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return True
    return False

def getSampleStatuses():
    statuses = {}
    with open("P01CA91955-WGS80-Full-Pilot-Samples.csv", "r") as csvfile:
        for lvec in csv.reader(csvfile):
            if "RandomID" in lvec:
                continue
            sample = lvec[6]
            statuses[sample] = {}
            statuses[sample]["patient"] = lvec[0]
            statuses[sample]["prog"] = lvec[1]
            statuses[sample]["gender"] = lvec[2]
            time = lvec[8]
            if time=="LongFollowUp":
                time = "T3"
            if time=="Index":
                time = "T2"
            if sample == "23521":
                time = "T2"
            statuses[sample]["time"] = time
    return statuses

def readMutations():
    mutations = {}
    with open(mutation_file, 'r') as csvfile:
        for lvec in csv.reader(csvfile):
            if "DNANum" in lvec[0]:
                continue
            (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
            if (is_snv=="f"):
                continue
            if (is_2p=="f"):
                continue
            if ("N" in sample):
                continue
            if sample not in samplePatientMap:
                continue
            refcnt = int(lvec[-2])
            bafcnt = int(lvec[-1])
            VAF = bafcnt/(refcnt+bafcnt)
            patient = samplePatientMap[sample][0]
    #        if patient not in patientSampleMap:
    #            patientSampleMap[patient] = set()
    #        patientSampleMap[patient].add(sample)
            if onlysomepatients and patient not in somepatients:
                continue
            if patient not in mutations:
                mutations[patient] = {}
                print("Reading mutations from patient", patient)
            if chr not in mutations[patient]:
                mutations[patient][chr] = {}
            pos = int(pos)
            if pos not in mutations[patient][chr]:
                mutations[patient][chr][pos] = {}
            if alt not in mutations[patient][chr][pos]:
                mutations[patient][chr][pos][alt] = {}
            mutations[patient][chr][pos][alt][sample] = VAF
    return mutations

def getCodeToSampleMap():
    lists = {}
    codes = {}
    alpha = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    with open("P01CA91955-WGS80-Full-Pilot-Samples.csv", "r") as csvfile:
        for lvec in csv.reader(csvfile):
            if "RandomID" in lvec:
                continue
            sample = lvec[6]
            patient = lvec[0]
            time = lvec[8]
            if time=="LongFollowUp":
                time = "T3"
            if time=="Index":
                time = "T2"
            if sample == "23521":
                time = "T2"
            if time in ["T1", "T2", "T3"]:
                if patient not in lists:
                    lists[patient] = []
                lists[patient].append(sample)
    for patient in lists:
        lists[patient].sort()
        codes[patient] = {}
        for n in range(len(lists[patient])):
            codes[patient][alpha[n]] = lists[patient][n]
    return codes

def readOptimInputs(codeToSampleMap):
    allsets = {}
    labelSamples = {}
    optimfiles = []
    for __, _, files in walk(optimdir):
        optimfiles += files
    for file in optimfiles:
        if "input" not in file:
            continue
        fvec = file.split("_")
        patient = fvec[0]
        if onlysomepatients and patient not in somepatients:
            continue
        tag = ""
        fvec = file.split(".")
        if len(fvec) > 1:
            tag = "_" + fvec[1]
        label = patient + tag
        labelSamples[label] = set()
        allsets[label] = {}
        #print(file, label)

        oneset = set()
        inputs = {}
        for line in open(optimdir + file, "r"):
            lvec = line.rstrip().split(" ")
            if len(lvec)==1:
                if len(oneset) > 0:
                    oneset = list(oneset)
                    oneset.sort()
                    oneset = tuple(oneset)
                    allsets[label][oneset] = inputs
                oneset = set()
                inputs = {}
                continue
            if "#" in line:
                continue
            if lvec[0] not in codeToSampleMap[patient]:
                assert(lvec[1] == "32" and lvec[2] == "799")
                lvec[0] = "c"
                lvec[3] = "x4"
            sample = codeToSampleMap[patient][lvec[0]]
            oneset.add(sample)
            labelSamples[label].add(sample)
            vaf = int(lvec[1])
            nmuts = int(lvec[2])
            tips = lvec[3]
            call = "1,1"
            if len(lvec)>4:
                call = lvec[4]
            tipvec = tips.split("_")
            tipset = set()
            for entry in tipvec:
                if "x" in entry:
                    tipset.add(entry)
            if sample not in inputs:
                inputs[sample] = {}
            tipset = list(tipset)
            tipset.sort()
            tipset = tuple(tipset)
            if tipset not in inputs[sample]:
                inputs[sample][tipset] = {}
            if call not in inputs[sample][tipset]:
                inputs[sample][tipset][call] = {}
                inputs[sample][tipset][call]["vaf"] = vaf
                inputs[sample][tipset][call]["nmuts"] = nmuts
            else:
                inputs[sample][tipset][call]["vaf"] = [inputs[sample][tipset][call]["vaf"], vaf]
                inputs[sample][tipset][call]["nmuts"] = [inputs[sample][tipset][call]["nmuts"], nmuts]
        if len(oneset) > 0:
            oneset = list(oneset)
            oneset.sort()
            oneset = tuple(oneset)
            allsets[label][oneset] = inputs
        if len(allsets[label]) == 1:
            print(label, "might not have spacing in it.")
    return allsets, labelSamples

def getOptimInputSplits(allsets):
    inputSplits = {}
    for label in allsets:
        for group in allsets[label]:
            for sample in allsets[label][group]:
                if len(allsets[label][group][sample]) > 1:
                    #assert len(allsets[label][group][sample]) == 2 ##Actually get 3 for 360!
                    if (label, group) not in inputSplits:
                        inputSplits[(label, group)] = {}
                    else:
                        print("Double splits for", label)
                    inputSplits[(label, group)][sample] = {}
                    for tips in allsets[label][group][sample]:
                        for call in allsets[label][group][sample][tips]:
                            vaf = allsets[label][group][sample][tips][call]["vaf"]
                            nmuts = allsets[label][group][sample][tips][call]["nmuts"]
                            if call not in inputSplits[(label, group)][sample]:
                                inputSplits[(label, group)][sample][call] = {}
                            inputSplits[(label, group)][sample][call][tips] = (vaf, nmuts)
        
    return inputSplits

def readSampleTrees(label):
    labelvec = label.split("_")
    treefile = treedir + labelvec[0] + ".newick"
    if len(labelvec)>1:
        treefile += "." + labelvec[1]
        assert(len(labelvec)==2)
    #print(label, "yields the newick file", treefile)
    newicks = []
    for line in open(treefile, "r"):
        if "#" in line:
            continue
        newicks.append(line.rstrip())
    halfindexes = range(round(len(newicks)/2))
    newicks = numpy.delete(newicks, halfindexes)
    trees = []
    for newick in newicks:
        tree = Tree(newick)
        for branch in tree.traverse():
            branch.dist = 0
        trees.append(tree)
    return trees

def isSplit(group, inputSplits, label):
    if (label, group) in inputSplits:
        return True
    return False

def getBranch(trees, tipnames):
    tipbranches = []
    for tree in trees:
        for tip in tree:
            for tipname in tipnames:
                if tipname in tip.name:
                    tipbranches.append(tip)
                elif "-1" in tipname:
                    if tip.name.split("_")[0] + "-1" == tipname:
                        tipbranches.append(tip)
        if len(tipbranches) == len(tipnames):
            return tree.get_common_ancestor(tipbranches)
        assert(len(tipbranches) == 0)
        tipbranches.clear()
    assert(False)
    return None

def getSigFile(branch, label, files):
    tipnames = []
    for tip in branch:
        tipnames.append(tip.name)
    tipnames.sort()
    sigfilename = label
    for tipname in tipnames:
        sigfilename += "_" + tipname
    sigfilename += ".mutations.tsv"
    if sigfilename not in files:
        sigfile = open(sigdir + sigfilename, "w")
        files[sigfilename] = sigfile
        sigfile.write("Chrom\tpos\talt\n")
    return files[sigfilename]

def closeFiles(files):
    for filename in files:
        files[filename].close()

def saveNonSplit(label, trees, group, allsets, files, chrom, pos, alt):
    if group not in allsets[label]:
        #At some point, we might save these somewhere, but for now they're just considered errors.
        return False
    tipnames = set()
    for sample in allsets[label][group]:
        for tips in allsets[label][group][sample]:
            for tip in tips:
                tipnames.add(combineSampleAndTipLabel(sample, tip))
    #print(label)
    assert(len(tipnames)>0)
    branch = getBranch(trees, tipnames)
    branch.dist += 1
    sigfile = getSigFile(branch, label, files)
    sigfile.write(chrom)
    sigfile.write("\t" + str(pos))
    sigfile.write("\t" + alt)
    sigfile.write("\n")
    return True

def combineSampleAndTipLabel(sample, tiplabel):
    nox = tiplabel.replace("x", "-")
    return sample + nox

def saveSingleSplit(allsets, label, group, keySample, keySampleTips, chrom, pos, alt, trees, files):
    tipnames = set()
    for sample in allsets[label][group]:
        if sample==keySample:
            for tip in keySampleTips:
                tipnames.add(combineSampleAndTipLabel(sample, tip))
        else:
            for tips in allsets[label][group][sample]:
                for tip in tips:
                    tipnames.add(combineSampleAndTipLabel(sample, tip))
    #print(label)
    assert(len(tipnames)>0)
    branch = getBranch(trees, tipnames)
    branch.dist += 1
    sigfile = getSigFile(branch, label, files)
    sigfile.write(chrom)
    sigfile.write("\t" + str(pos))
    sigfile.write("\t" + alt)
    sigfile.write("\n")
    return True
    

def saveForSplits(splits, doublesplits, group, sampleVAFs, chrom, pos, alt, CNVs, patient, splitPartitions):
    for keySample in splitPartitions:
        call = getCNVCall(patient, keySample, chrom, pos, CNVs)
        strcall = str(call[0]) + "," + str(call[1])
        if len(splitPartitions) > 1:
            #Save double splits separately.
            if group not in doublesplits:
                doublesplits[group] = {}
            if keySample not in doublesplits[group]:
                doublesplits[group][keySample] = {}
            if strcall not in doublesplits[group][keySample]:
                doublesplits[group][keySample][strcall] = []
            doublesplits[group][keySample][strcall].append((sampleVAFs[keySample], chrom, pos, alt))
        else:
            if group not in splits:
                splits[group] = {}
            if strcall not in splits[group]:
                splits[group][strcall] = []
            splits[group][strcall].append((sampleVAFs[keySample], chrom, pos, alt))

def divideLengthFromSplit(lengthOnly, tipsets, trees, keySample):
    totlen = 0
    branches = {}
    for tipset in tipsets:
        alltips = set()
        for tip in tipset:
            alltips.add(combineSampleAndTipLabel(keySample, tip))
        branch = getBranch(trees, alltips)
        totlen += branch.dist
        branches[tipset] = branch
    for tipset in branches:
        branch = branches[tipset]
        branch.dist += lengthOnly * (branch.dist/totlen)

def saveSplits(splits, trees, files, inputSplits, label):
    for group in splits:
        lengthOnly = 0
        tipsets = set()
        for call in splits[group]:
            vafvec = splits[group][call]
            splitinfo = inputSplits[(label, group)]
            keySample = list(splitinfo.keys())[0]
            if call not in splitinfo[keySample]:
                #We didn't say how to split this call, so the only thing we can do is 
                # add it to the length, but not the signature file
                lengthOnly += 1
                continue
            splitinfo = splitinfo[keySample][call]
#            if len(splitinfo)==1:
#                #Just save everything to the appropriate branch
#                for (vaf, chrom, pos, alt) in vafvec:
#                    saveSingleSplit(allsets, label, group, keySample, list(splitinfo.keys())[0], chrom, pos, alt, trees, files)
#                continue
            vnts = []
            muttot = 0
            for tips in splitinfo:
                tipsets.add(tips)
                (vaf, nmut) = splitinfo[tips]
                if isinstance(vaf, list):
                    #The call is unbalanced, leading to two different VAFs for the same set of tips.
                    for i in range(len(vaf)):
                        vnts.append((vaf[i], nmut[i], tips))
                        muttot += nmut[i]
                else:
                    vnts.append((vaf, nmut, tips))
                    muttot += nmut
            vafvec.sort()
            vnts.sort()
            currindex = 0
            currend = 0
            for (vaf, nmut, tips) in vnts:
                currend += round((nmut/muttot)*len(vafvec))
                while currindex < currend and currindex < len(vafvec):
                    (__, chrom, pos, alt) = vafvec[currindex]
                    saveSingleSplit(allsets, label, group, keySample, tips, chrom, pos, alt, trees, files)
                    currindex += 1
            if currend < len(vafvec):
                assert len(vafvec)-currend == 1
                #just in case we don't sum to 100%
                (__, chrom, pos, alt) = vafvec[-1]
                saveSingleSplit(allsets, label, group, keySample, tips, chrom, pos, alt, trees, files)
        divideLengthFromSplit(lengthOnly, tipsets, trees, keySample)

def divideLengthAmongBranches(divlen, branches, families):
    totlen = 0
    for tiplists in families:
        totlen += len(families[tiplists])
    for tiplists in families:
        branch = branches[tiplists]
        branch.dist += len(families[tiplists])
        branch.dist += round(divlen * len(families[tiplists])/totlen)
    

def saveDoubleSplitFamilies(families, trees, files, label, lengthOnlyPositions):
    alltips = set()
    tipsets = {}
    for tiplists in families:
        onetipset = set()
        for sample, tips in tiplists:
            for tip in tips:
                alltips.add(combineSampleAndTipLabel(sample, tip))
                onetipset.add(combineSampleAndTipLabel(sample, tip))
        tipsets[tiplists] = onetipset
    basebranch = getBranch(trees, alltips)
    badtiplists = []
    branches = {}
    for tiplists in tipsets:
        listbasebranch = getBranch(trees, tipsets[tiplists])
        if listbasebranch == basebranch and alltips != tipsets[tiplists]:
            lengthOnlyPositions += len(families[tiplists])
            badtiplists.append(tiplists)
        else:
            branches[tiplists] = listbasebranch
    for tiplists in badtiplists:
        del tipsets[tiplists]
        del families[tiplists]
    for tiplists in branches:
        branch = branches[tiplists]
        sigfile = getSigFile(branch, label, files)
        for (chrom, pos, alt) in families[tiplists]:
            sigfile.write(chrom)
            sigfile.write("\t" + str(pos))
            sigfile.write("\t" + alt)
            sigfile.write("\n")
    divideLengthAmongBranches(lengthOnlyPositions, branches, families)


def sortDoubleSplits(doublesplits, trees, files, inputSplits, label):
    #For the double splits, they're not stored as vectors, so we have to store them.
    for group in doublesplits:
        positions = {}
        lengthOnlyPositions = set()
        for keySample in doublesplits[group]:
            for call in doublesplits[group][keySample]:
                vafvec = doublesplits[group][keySample][call]
                splitinfo = inputSplits[(label, group)][keySample]
                if call not in splitinfo:
                    #We didn't say how to split this call, so the only thing we can do is 
                    # add it to the length, but not the signature file
                    for (vaf, chrom, pos, alt) in vafvec:
                        location = (chrom, pos, alt)
                        if location in positions:
                            del positions[location]
                        lengthOnlyPositions.add(location)
                    continue
                splitinfo = splitinfo[call]
                vnts = []
                muttot = 0
                for tips in splitinfo:
                    (vaf, nmut) = splitinfo[tips]
                    if isinstance(vaf, list):
                        #The call is unbalanced, leading to two different VAFs for the same set of tips.
                        for i in range(len(vaf)):
                            vnts.append((vaf[i], nmut[i], tips))
                            muttot += nmut[i]
                    else:
                        vnts.append((vaf, nmut, tips))
                        muttot += nmut
                vafvec.sort()
                vnts.sort()
                currindex = 0
                currend = 0
                for (vaf, nmut, tips) in vnts:
                    currend += round((nmut/muttot)*len(vafvec))
                    while currindex < currend and currindex < len(vafvec):
                        (__, chrom, pos, alt) = vafvec[currindex]
                        currindex += 1
                        location = (chrom, pos, alt)
                        if location in lengthOnlyPositions:
                            continue
                        if location not in positions:
                            positions[location] = []
                        positions[location].append((keySample, tips))
                if currend < len(vafvec):
                    assert len(vafvec)-currend == 1
                    #just in case we don't sum to 100%
                    (__, chrom, pos, alt) = vafvec[-1]
                    if location not in positions:
                        positions[location] = []
                    positions[location].append((keySample, tips))
        families = {}
        for location in positions:
            positions[location].sort()
            alltips = tuple(positions[location])
            if alltips not in families:
                families[alltips] = []
            families[alltips].append(location)
        saveDoubleSplitFamilies(families, trees, files, label, len(lengthOnlyPositions))
        if needDoubleSplitInfo:
            print("Double splits for", label, "group", str(group))
            tiplists = list(families.keys())
            tiplists.sort()
            for alltips in tiplists:
                print(str(alltips), "\t", str(len(families[alltips])))
            print("Unusable:\t", len(lengthOnlyPositions))

def writeTrees(trees, label):
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    print("Tree(s) for", label)
    outtrees = open(lengthdir + label + ".newick", "w")
    for tree in trees:
        line = tree.write(format=1)
        outtrees.write(line + "\n")
        if (showTrees):
            for branch in tree.traverse():
                if branch.dist != 0:
                    text_face = TextFace(str(int(branch.dist)))
                    branch.add_face(text_face, column=1, position="branch-top")
            tree.show(tree_style = ts)
    outtrees.close()
    

def sortMutations(mutations, allsets, inputSplits, deletions, CNVs, labelSamples):
    savesort = open(sigdir + "saved_or_skipped.tsv", "w")
    savesort.write("Tree\tsaved\tskipped\tsplit\n")
    for label in allsets:
        #print("Working on", label)
        trees = readSampleTrees(label)
        patient = label.split("_")[0]
        samples = labelSamples[label]
        muts = mutations[patient]
        splits = {}
        doublesplits = {}
        files = {}
        used = 0
        skipped = 0
        split = 0
        for chrom in muts:
            for pos in muts[chrom]:
                for alt in muts[chrom][pos]:
                    sampleVAFs = muts[chrom][pos][alt]
                    group = list(sampleVAFs.keys())
                    toRemove = []
                    for sample in group:
                        if sample not in samples:
                            toRemove.append(sample)
                    for sample in toRemove:
                        group.remove(sample)
                    if len(group) < 2:
                        continue
                    group.sort()
                    group = tuple(group)
                    if isSplit(group, inputSplits, label):
                        saveForSplits(splits, doublesplits, group, sampleVAFs, chrom, pos, alt, CNVs, patient, inputSplits[(label, group)])
                        split += 1
                    else:
                        if saveNonSplit(label, trees, group, allsets, files, chrom, pos, alt):
                            used += 1
                        else:
                            skipped += 1
        sortDoubleSplits(doublesplits, trees, files, inputSplits, label)
        saveSplits(splits, trees, files, inputSplits, label)
        closeFiles(files)
        #print("For", label, ",", str(used), "used mutations, and", str(skipped), "skipped mutations.")
        savesort.write(label)
        savesort.write("\t" + str(used))
        savesort.write("\t" + str(skipped))
        savesort.write("\t" + str(split))
        savesort.write("\n")
        writeTrees(trees, label)
    savesort.close()

codeToSampleMap = getCodeToSampleMap()

#Read in the 'optim' input files, and parse them to figure out where the splits are.
allsets, labelSamples = readOptimInputs(codeToSampleMap)
inputSplits = getOptimInputSplits(allsets)

#Read in the deletions/CNVs
(patientSampleMap, samplePatientMap) = lps.getPatientSampleMap()
deletions, CNVs = lps.loadDeletionsAndCNVs(samplePatientMap)

#Read in the mutations
mutations = readMutations()

sortMutations(mutations, allsets, inputSplits, deletions, CNVs, labelSamples)


#writeAllSampleVAFs(mutations, patientSampleMap, deletions)

