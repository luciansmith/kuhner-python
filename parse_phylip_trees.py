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

phylipdir  = "phylip_input/"


filenames = []
for (_, _, f) in walk(phylipdir):
    filenames += f
    break

for f in filenames:
    if "outtree" not in f:
        continue
    (patient, __) = f.split("_")
    print("Reading", f)
    tree = Tree(phylipdir + f)
    tree.set_outgroup("blood")
    print(tree)
    subpopfreqs = 0
    nodelist = {}
    ages = set()
    for node in tree.traverse("levelorder"):
        n = node.name
#        if (n=="blood"):
#            sample = n
#            code = "B00"
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
    print(agelist)
    keynodes = []
    for n in nodelist:
        if nodelist[n]["age"] in agelist:
            keynodes.append(n)
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
    treetype = "unknown"
    if len(unique_parents) == 2:
        #Either paired coherent or paired incoherent
        for nodepair in parents:
            if parents[nodepair] in unique_parents:
                if nodelist[nodepair[0]]["age"] == nodelist[nodepair[1]]["age"]:
                    treetype = "paired coherent"
                else:
                    treetype = "paired incoherent"
                break
    else:
        #Four options: paired t1, paired t2, incoherent t1 out, inocoherent t2 out
        for nodepair in parents:
            parent = parents[nodepair]
            if parent in unique_parents:
                if nodelist[nodepair[0]]["age"] == nodelist[nodepair[1]]["age"]:
                    if nodelist[nodepair[0]]["age"] == agelist[0]:
                        treetype = "paired T1 only"
                    else:
                        treetype = "paired T2 only"
                else:
                    for node in nodelist:
                        if node in nodepair:
                            continue
                        if parent.up == nodelist[node]["node"].up:
                            if nodelist[node]["age"] == agelist[0]:
                                treetype = "incoherent, T2 out"
                            else:
                                treetype = "incoherent, T1 out"
                break
    print(treetype)