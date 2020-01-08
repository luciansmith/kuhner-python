#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:29:58 2019

@author: lpsmith
"""

import os
import ete3
import csv
from os import walk
from os import path
from os import mkdir

#outdir = "non17_phylogenies/"

onlysomepatients = False
#somepatients = [88, 130, 385, 541, 572, 729, 956, 42, 169, 279, 623, 728, 852, 951, 995]
somepatients = [422]

usepilot = True

#treedir = "phylip_input_noPilot/"
#outdir = "all_phylogenies_noPilot/"

treedir = "phylip_TS_analysis/"
outdir = "all_phylogenies/"

driverfile = "EA_drivers_by_sample.txt"
#driverfile = "EA_drivers_pilot.txt"

samplemapfile = "20181031_SampleCodeLucianTrees_agesDiffRounding_withPilot.txt"

checkDelete = False
delfilename = "check_for_deletions.tsv"

delresults = "found_deletions_or_not.txt"

displayDeletions = True
if displayDeletions:
    outdir = "all_phylogenies_withdels/"

svdir = "union.event.by.sample/"
displaySVs = True
if displaySVs:
    outdir = "all_phylogenies_SVs/"
    outfilename = "SV_matches.txt"


displayCombo = False
if displayCombo:
    outdir = "all_phylogenies_SVs_and_drivers/"

if displaySVs and not displayCombo:
    SV_out = open(outdir + outfilename, "w")
    SV_out.write("Patient\tSV_id\tMatch/Miss\n");


if not path.isdir(outdir):
    mkdir(outdir)

target_pids = [55, 59, 478, 635, 865, 126, 909, 381, 609, 184]



SVfiles = []
for __, _, files in walk(svdir):
    SVfiles += files



g_check_loss = set()

def readInEADrivers():
    drivers = {}
    gene2chrom = {}
    for line in open(driverfile, "r"):
        if "GeneName" in line:
            continue
        lvec = line.rstrip().split()
        (gene, chrom, patient, pos, alt, sample) = lvec[0:6]
        if gene not in gene2chrom:
            gene2chrom[gene] = chrom
        patient = int(patient)
        if patient not in drivers:
            drivers[patient] = {}
        if gene not in drivers[patient]:
            drivers[patient][gene] = {}
        if (pos, alt) not in drivers[patient][gene]:
            drivers[patient][gene][(pos, alt)] = []
        drivers[patient][gene][(pos, alt)].append(sample)
    return drivers, gene2chrom

def readInSVs():
    SVs = {}
    for file in SVfiles:
        if ".txt" not in file:
            continue
        patient = file.split(".")[0]
        patient = int(patient)
        labels = []
        with open(svdir + file, 'r') as csvfile:
            for lvec in csv.reader(csvfile):
                if int(lvec[0]) > 100:
                    try:
                        for n in range(len(lvec)):
                            x = int(lvec[n])
                            labels.append(lvec[n])
                    except:
                        continue
                svid = lvec[-1] + "_" + lvec[-3] + "_" + lvec[-2]
                if displayCombo:
                    if "dup" in svid:
                        continue
                    if "del" in svid:
                        continue
                samples = []
                for n in range(len(labels)):
                    if lvec[n] == "1":
                        samples.append(labels[n])
                if len(samples) == 0:
                    continue
                if patient not in SVs:
                    SVs[patient] = {}
                SVs[patient][svid] = samples
    return SVs

def readInDeletions():
    deletions = {}
    for line in open(delresults, "r"):
        if "Patient" in line:
            continue
        (patient, sample, chrom, pos, gene, incong, deleted) = line.rstrip().split()
        if deleted=="True":
            patient = int(patient)
            if patient not in deletions:
                deletions[patient] = {}
            if gene not in deletions[patient]:
                deletions[patient][gene] = {}
            if pos not in deletions[patient][gene]:
                deletions[patient][gene][pos] = set()
            deletions[patient][gene][pos].add(sample)
    return deletions

def getBranchOrBranchesFor(samples, tree, incongruous_samples, alt_samples, dels):
    #print(str(samples))
    tips = set()
    othertips = set()
    rettips = False
    for tip in tree:
        if tip.name in samples:
            tips.add(tip)
        elif tip.name in dels:
            continue
        else:
            othertips.add(tip)
            alt_samples.add(tip.name)
    common = tree.get_common_ancestor(tips)
    if len(samples)==1:
        return tips
    for alt in othertips:
        superset = tips.copy()
        superset.add(alt)
        if tree.get_common_ancestor(superset) == common:
            rettips = True
            incongruous_samples.add(alt.name)
#            return tips
    if rettips:
        return tips
    return [common]

def getSampleNameMap():
    sampleNameMap = {}
    for line in open(samplemapfile, "r"):
        if "Patient" in line:
            continue
        (__, sample, __, age, level, __, label, pnp, __, __, gej, os) = line.rstrip().split()
        level = int(level)
        gej = int(gej)
        level = gej-level
        level = str(level)
        if label=="Index":
            label = "T2"
        elif label=="LongFollowUp":
            label = "T3"
        if label=="GastricSample":
            label = "Gastric"
            sampleNameMap[sample] = label
        else:
            sampleNameMap[sample] = label + ", " + level + " cm"
    return sampleNameMap

def displayDriversOnPhylogenies(pid, t, drivers, deletions):
    if pid in drivers:
        for gene in drivers[pid]:
            append = ""
            if len(drivers[pid][gene]) > 1:
                append = "_a"
            for posalt in drivers[pid][gene]:
                text_face = ete3.TextFace(gene + append, fsize=vertical_margins[pid]-5)
                incongruous_samples = set()
                alt_samples = set()
                dels = set()
                if pid in deletions and gene in deletions[pid] and posalt[0] in deletions[pid][gene]:
                    dels = deletions[pid][gene][posalt[0]]
                branches = getBranchOrBranchesFor(drivers[pid][gene][posalt], t, incongruous_samples, alt_samples, dels)
                if len(branches) > 1:
                    text_face = ete3.TextFace(gene + append, fsize=vertical_margins[pid]-5, fgcolor="orangered")
                for branch in branches:
                    branch.add_face(text_face, column=1, position="branch-top")
                if checkDelete:
                    for sample in alt_samples:
                        del_look_file.write(str(pid))
                        del_look_file.write("\t" + sample)
                        del_look_file.write("\t" + gene2chrom[gene])
                        del_look_file.write("\t" + posalt[0])
                        del_look_file.write("\t" + gene)
                        del_look_file.write("\t" + str(sample in incongruous_samples))
                        del_look_file.write("\n")
                if displayDeletions and len(dels)>0 and len(drivers[pid][gene][posalt])>0:
                    text_face = ete3.TextFace(u'Δ ' + gene + append, fsize=vertical_margins[pid]-5)
                    #print(str(dels))
                    branches = getBranchOrBranchesFor(dels, t, set(), set(), set())
                    for branch in branches:
                        branch.add_face(text_face, column=1, position="branch-top")
                if append == "_a":
                    append = "_b"
                elif append == "_b":
                    append = "_c"
                elif append == "_c":
                    append = "_d"
                elif append == "_d":
                    append = "_e"
    else:
        print("No", pid, "in drivers")

def displaySVsOnPhylogenies(pid, t, SVs):
    if pid in SVs:
        for SV in SVs[pid]:
            SVname = SV
            text_face = ete3.TextFace(SVname, fsize=vertical_margins[pid]-5)
            branches = getBranchOrBranchesFor(SVs[pid][SV], t, set(), set(), set())
            if len(branches) > 1:
                text_face = ete3.TextFace(SVname, fsize=vertical_margins[pid]-5, fgcolor="orangered")
            for branch in branches:
                branch.add_face(text_face, column=1, position="branch-top")
            #Collect data for statistics
            samples = SVs[pid][SV]
            ntips = 0;
            for tip in t:
                ntips += 1
            if not displayCombo and len(samples) > 1 and len(samples) < ntips:
                SV_out.write(str(pid))
                SV_out.write("\t" + SV)
                if len(branches)>1:
                    SV_out.write("\tMiss")
                else:
                    SV_out.write("\tMatch")
                SV_out.write("\n")
    else:
        print("No", pid, "in SVs")

# Ingest the trees.  First, grab all the filenames/ paths of the tree files
tree_files = {'path': [], 'PatientID': []}
for path, _, files in os.walk(treedir):
    for file_ in files:
        if "b_" in file_:
            continue
        if file_[-8:] == 'tree.txt':
            tree_files['path'].append(os.path.join(path, file_))
            tree_files['PatientID'].append(file_.split('_')[0])
for x, pid in enumerate(tree_files['PatientID']):
    if pid == "891a":
        pid = 891
    tree_files['PatientID'][x] = int(pid)
    
vertical_margins = {
#        88:  25,
#        130: 20, 
#        385: 20, 
#        541: 20, 
#        572: 50, 
#        729: 55, 
#        956: 20, 
#        42:  25, 
#        169: 28, 
#        279: 20, 
#        623: 30, 
#        728: 20, 
#        852: 25, 
#        951: 25, 
#        995: 25,
#
#        55: 24, 
#        59: 18, 
#        126: 19, 
#        184: 12, 
#        381: 16, 
#        478: 12, 
#        609: 18, 
#        635: 20, 
#        865: 15, 
#        909: 15,   
#
#        286: 18, 
#        387: 26, 
#        521: 24, 
#        672: 33, 
}

drivers, gene2chrom = readInEADrivers()
deletions = readInDeletions()
sampleNameMap = getSampleNameMap()
SVs = readInSVs()

if checkDelete:
    del_look_file = open(delfilename, "w")
    del_look_file.write("Patient\tSample\tchrom\tpos\tgene\tincongruous?\n")

for x, path in enumerate(tree_files['path']):
    pid = tree_files['PatientID'][x]
    if pid not in vertical_margins:
        vertical_margins[pid] = 20
    if onlysomepatients and pid not in somepatients:
        continue
    print(pid)
    # And load the tree
    t = ete3.Tree(path)

    # Need to correlate leaf names to DNANums, so grab the leaf names
    #leaves = t.get_leaf_names()

    # Grab just the DNANums of the leaves
    #dnanums = [x.split('_')[0] for x in leaves if 'blood' != x[:5]]

    # Define the tree style
    tstyle = ete3.TreeStyle()
    # Make it a rectangular tree
    tstyle.mode = 'r'
    tstyle.branch_vertical_margin = vertical_margins[pid]
    tstyle.show_scale = False
    # Don't show the leaf names or print the branch lengths
    tstyle.show_leaf_name = False
    tstyle.show_branch_length = False
    for branch in t:
        if "N" in branch.name:
            branch.delete()
    for branch in t:
        if "blood" in branch.name:
            t.set_outgroup(branch)
            root = t.get_tree_root()
            root.dist = 2*branch.dist
            branch.delete()
            root.get_children()[0].delete()
    for branch in t.traverse():
        if "_" in branch.name:
            #Print the leaf names in our own font
            name_face = ete3.AttrFace("name", fsize=vertical_margins[pid]+10)
            branch.add_face(name_face, column=0, position="branch-right")
            branch.name = branch.name.split("_")[0]
#        else:
    if displaySVs:
        displaySVsOnPhylogenies(pid, t, SVs)
    else:
        displayDriversOnPhylogenies(pid, t, drivers, deletions)

    if displayCombo and displaySVs:
        displayDriversOnPhylogenies(pid, t, drivers, deletions)
        
    #        else:
    #            branch.add_feature(stuff="I'm internal!")
    #            upper = ete3.AttrFace("stuff", fsize=vertical_margins[pid])
    #            lower = ete3.TextFace("below", fsize=vertical_margins[pid])
    #            
    ##            # Set some attributes
    ##            hola.margin_top = 10
    ##            hola.margin_right = 10
    ##            hola.margin_left = 10
    ##            hola.margin_bottom = 10
    ##            hola.opacity = 0.5 # from 0 to 1
    ##            hola.inner_border.width = 1 # 1 pixel border
    ##            hola.inner_border.type = 1  # dashed line
    ##            hola.border.width = 1
    ##            hola.background.color = "LightGreen"
    #            
    #            branch.add_face(upper, column=0, position = "branch-top")
    #            branch.add_face(lower, column=1, position = "branch-bottom")
            
    # Stretch in x
    tstyle.scale = 1000

    t.render(outdir + str(pid) + "_sampleIDs_tree.png", tree_style=tstyle)
    t.write(outdir + str(pid) + "_tree.txt", format=1)
    for tip in t:
        tip.name = sampleNameMap[tip.name]

    t.render(outdir + str(pid) + "_timespaceIDs_tree.png", tree_style=tstyle)

if checkDelete:
    del_look_file.close()
print("done!")

if displaySVs and not displayCombo:
    SV_out.close()