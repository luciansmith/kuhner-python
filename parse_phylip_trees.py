#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 12:47:03 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os.path import isfile
import numpy
from ete3 import Tree


def readProgressorOrNot():
    progressorMap = {}
    pfile = open("Patient_status.txt", "r")
    for line in pfile:
        if "Patient" in line:
            continue
        (patient, pstat, sex) = line.rstrip().split()
        progressorMap[patient] = pstat
    return progressorMap

progressorMap = readProgressorOrNot()

phylipdir  = "phylip_input_lowVAF/"
summaryfile = open(phylipdir + "tree_types.txt", "w")
summaryfile.write("Patient")
summaryfile.write("\tT1")
summaryfile.write("\tT2")
summaryfile.write("\tProgressionStatus")
summaryfile.write("\tAge Type")
summaryfile.write("\tLevel Type")
summaryfile.write("\n")

filenames = []
for (_, _, f) in walk(phylipdir):
    filenames += f
    break

for f in filenames:
    if "outtree" not in f:
        continue
    (patient, __) = f.split("_")
    if "a" in patient or "b" in patient or "c" in patient:
        patient = patient[0:3]
    print("Reading", f)
    tree = Tree(phylipdir + f)
    bloodstr = "blood";
    if progressorMap[patient] == "NP":
        bloodstr += "_non"
    else:
        bloodstr += "_prog"
    tree.set_outgroup(bloodstr)
    #print(tree)
    subpopfreqs = 0
    nodelist = {}
    ages = set()
    for node in tree.traverse("levelorder"):
        n = node.name
        if "blood" in n:
            continue
        if "_" not in n:
            continue
#            sample = n
#            code = "I01"
        (sample, code) = n.split("_")
        nodelist[n] = {}
        nodelist[n]["node"] = node
        nodelist[n]["sample"] = sample
        level = code[0]
        age = code[1:3]
        age = int(age)
        nodelist[n]["age"] = age
        nodelist[n]["level"] = level
        ages.add(age)
    agelist = list(ages)
    agelist.sort()
    agelist = agelist[0:2]
    if patient=="391":
        agelist = [60, 66]
    if patient=="611":
        agelist = [69, 73]
    print(agelist)
    keynodes = []
    oldtree = tree
    prunelist = [bloodstr]
    S_list = []
    T_list = []
    for n in nodelist:
        if nodelist[n]["age"] in agelist:
            if "N" in nodelist[n]["node"].name:
                continue
            keynodes.append(n)
            prunelist.append(nodelist[n]["node"].name)
    tree.prune(prunelist)
    for name in prunelist:
        if "blood" in name:
            continue
        if "S" in name:
            S_list.append(name)
        else:
            assert("T" in name)
            T_list.append(name)
    space = "equal"
    if len(S_list) > len(T_list):
        space = "three S"
    elif len(S_list) < len(T_list):
        space = "three T"
    if len(S_list)==0 or len(T_list)==0:
        space = "all same"
    print("Pruned version for patient", patient, tree)
    parents = {}
    for n1 in range(0,len(keynodes)-1):
        for n2 in range(n1+1, len(keynodes)):
            node1 = keynodes[n1]
            node2 = keynodes[n2]
            parents[(node1, node2)] = tree.get_common_ancestor(node1, node2)
    unique_parents = []
    all_parents = []
    for nodepair in parents:
        parent = parents[nodepair]
        if parent in unique_parents and parent in all_parents:
            unique_parents.remove(parent)
        elif parent not in all_parents:
            unique_parents.append(parent)
        all_parents.append(parent)
    age_treetype = "unknown"
    level_treetype = "unknown"
    if len(unique_parents) == 2:
        #Either paired coherent or paired incoherent
        for nodepair in parents:
            if parents[nodepair] in unique_parents:
                if nodelist[nodepair[0]]["age"] == nodelist[nodepair[1]]["age"]:
                    age_treetype = "age paired coherent"
                else:
                    age_treetype = "age paired incoherent"
                if space=="equal":
                    if nodelist[nodepair[0]]["level"] == nodelist[nodepair[1]]["level"]:
                        level_treetype = "level paired coherent"
                    else:
                        level_treetype = "level paired incoherent"
                else:
                    level_treetype = "two pairs, " + space
                break
    else:
        #Four options: paired t1, paired t2, incoherent t1 out, inocoherent t2 out
        for nodepair in parents:
            parent = parents[nodepair]
            if parent in unique_parents:
                if nodelist[nodepair[0]]["age"] == nodelist[nodepair[1]]["age"]:
                    if nodelist[nodepair[0]]["age"] == agelist[0]:
                        age_treetype = "age paired T1 only"
                    else:
                        age_treetype = "age paired T2 only"
                    if nodelist[nodepair[0]]["level"] == "S":
                        level_treetype = "age paired Stomach only"
                    else:
                        level_treetype = "age paired Teeth only"
                else:
                    for node in nodelist:
                        if node in nodepair:
                            continue
                        if parent.up == nodelist[node]["node"].up:
                            if nodelist[node]["age"] == agelist[0]:
                                age_treetype = "age incoherent, T2 out"
                            else:
                                age_treetype = "age incoherent, T1 out"
                if space == "equal":
                    if nodelist[nodepair[0]]["level"] == nodelist[nodepair[1]]["level"]:
                        if nodelist[nodepair[0]]["level"] == "S":
                            level_treetype = "level paired Stomach only"
                        else:
                            level_treetype = "level paired Teeth only"
                    else:
                        for node in nodelist:
                            if node in nodepair:
                                continue
                            if parent.up == nodelist[node]["node"].up:
                                if nodelist[node]["level"] == "S":
                                    level_treetype = "level incoherent, Teeth out"
                                else:
                                    level_treetype = "level incoherent, Stomach out"
                else:
                    if not(nodelist[nodepair[0]]["level"] == nodelist[nodepair[1]]["level"]):
                        level_treetype = "stairs tips incoherent, " + space
                    else:
                        for node in nodelist:
                            if node in nodepair:
                                continue
                            if parent.up == nodelist[node]["node"].up:
                                if nodelist[node]["level"] == nodelist[nodepair[0]]["level"]:
                                    level_treetype = "stairs coherent, " + space
                                else:
                                    level_treetype = "stairs middle incoherent, " + space
                    
                    
                break
    summaryfile.write(patient)
    summaryfile.write("\t" + str(agelist[0]))
    summaryfile.write("\t" + str(agelist[1]))
    summaryfile.write("\t" + progressorMap[patient])
    summaryfile.write("\t" + age_treetype)
    summaryfile.write("\t" + level_treetype)
    summaryfile.write("\n")
    
summaryfile.close();