#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 12:47:03 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
import re
from ete2 import Tree

BEAST_dir = "BEAST_output/"
ploidy_dir = "summary_stats/v6/ploidy/"
odds_file = "joint_summary_v6_odds.txt"

odds_cutoff = 0.25

def makeRatesNames(tree):
    rintro = re.compile(':\[\&rate')
    tree = rintro.sub("_rate", tree)
    equals = re.compile('=')
    tree = equals.sub("_", tree)
    point = re.compile('(_\d)\.')
    tree = point.sub(r'\1_', tree)
    point = re.compile('(_\d*)E-')
    tree = point.sub(r'\1En', tree)
    point = re.compile(']')
    tree = point.sub(r':', tree)
    return tree

def loadOdds():
    ofile = open(odds_file, "r")
    odds = {}
    for line in ofile:
        if line.find("patient") != -1:
            continue
        lvec = line.split()
        sample = lvec[1]
        oddsval = lvec[20]
        if (oddsval == "auto-diploid"):
            odds[sample] = 1
        elif (oddsval == "auto-tetraploid"):
            odds[sample] = 0
        else:
            odds[sample] = float(oddsval)
    return odds

def loadPloidies(patient):
    ploidies = {}
    pfile = open(ploidy_dir + "diploid/" + patient + "_fcn_ascat_ploidy.txt", "r")
    for line in pfile:
        if line.find("x") != -1:
            continue
        print line.split()
        (name, ploidy) = line.split()
        (patient, sample) = name.split("_")
        ploidies[sample, "diploid"] = float(ploidy)
    pfile = open(ploidy_dir + "tetraploid/" + patient + "_fcn_ascat_ploidy.txt", "r")
    for line in pfile:
        if line.find("x") != -1:
            continue
        (name, ploidy) = line.split()
        (patient, sample) = name.split("_")
        ploidies[sample, "tetraploid"] = float(ploidy)
    return ploidies

def getAncestorLength(node):
    l = 0
    for parent in node.get_ancestors():
        l += parent.dist
    return l
#test = "(1:[&rate=1.2340667318115055E-6]12.219250631991326,(3:[&rate=0.05056188882428441]5.8063728510850865,(4:[&rate=0.05056188882428441]3.0702802218454357,2:[&rate=0.05056188882428441]1.0602805103399353):[&rate=0.05056188882428441]2.736092629239651):[&rate=0.05056188882428441]8.42287749241174);"

#test2 = "(1:12.219250631991326,(3:5.8063728510850865,(4:3.0702802218454357,2:1.0602805103399353):2.736092629239651):8.42287749241174);"

#Tree(test2)

#tree = makeRatesNames(test)
#Tree(tree, format=1)

filenames = []
for (_, _, f) in walk(BEAST_dir):
    filenames += f
    break

odds = loadOdds()


for f in filenames:
    if f.find(".tree") == -1:
        continue
    patient = f.split("_")[0]
    ploidies = loadPloidies(patient)
    #print "Reading", f
    bfile = open(BEAST_dir + f, "r")
    best_tree = ""
    best_posterior = float("-inf")
    namesamplemap = {}
    for line in bfile:
        if line[0:4] != "tree":
            lvec = line.split()
            if len(lvec) == 2 and len(lvec[0]) == 1:
                namesamplemap[lvec[0]] = lvec[1].replace(",","")
            continue
        fullvec = line.rstrip().split(" ")
        posterior = re.split("[=,]+", fullvec[2])[1]
        posterior = float(posterior)
        if posterior > best_posterior:
            best_posterior = posterior
            best_tree = fullvec[5]
    best_tree = makeRatesNames(best_tree)
    tree = Tree(best_tree, format=1)
    
    tetraploids = []
    tip_tetraploids = []
    all_tetraploids = []
    all_diploids = []
    rates = {}
    for node in tree.traverse("postorder"):
        n = node.name
        nvec= n.split("_")
        rate = 0
        if len(nvec) == 4:
            name = nvec[0]
            postpoint = nvec[3].replace("n", "-")
            rate = nvec[2] + "." + postpoint
            rate = float(rate)
        rates[node] = rate;
        print nvec, rate
        if name != "":
            sample = namesamplemap[name]
            samp_odds = odds[sample]
            ploidy = -1
            if (samp_odds < odds_cutoff):
                ploidy = ploidies[sample, "tetraploid"]
            elif samp_odds > 1-odds_cutoff:
                ploidy = ploidies[sample, "diploid"]
            elif f.find("tetra") != -1:
                #'tetra' in the filename means this run had the tetraploid version
                ploidy = ploidies[sample, "tetraploid"]
            else:
                ploidy = ploidies[sample, "diploid"]
            if ploidy > 3.0:
                tetraploids.append(node)
                tip_tetraploids.append(node)
                all_tetraploids.append(node)
            else:
                all_diploids.append(node)
        else:
            #internal node
            all_tet = True
            all_dip = True
            for child in node.children:
                if child not in tetraploids:
                    all_tet = False
                if child not in all_diploids:
                    all_dip = False
            if all_tet:
                for child in node.children:
                    tetraploids.remove(child)
                tetraploids.append(child)
                all_tetraploids.append(node)
            if all_dip:
                all_diploids.append(node)
    print len(tip_tetraploids), "total tetraploids in tree, with", len(tetraploids), "parsimonious nodes."
    print len(all_diploids), "total diploids nodes in tree."
    
    lengths = {}
    for node in tree.traverse("preorder"):
        start = getAncestorLength(node)
        end = start + node.dist
        lengths[node] = (start, end)
             
    all_rates = [[], [], []]
    for node in tree.traverse("postorder"):
        if (rates[node] == 0):
            continue
        if node in all_diploids:
            all_rates[0].append((rates[node], lengths[node]))
        elif node in all_tetraploids:
            all_rates[1].append((rates[node], lengths[node]))
        else:
            all_rates[2].append((rates[node], lengths[node]))
