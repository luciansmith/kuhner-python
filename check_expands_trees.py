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
from ete2 import Tree

directories = ["expands_full_input_alldupes_noLOH/results/"]
for n in range(0,100):
    directories.append("shuffled_output/expands_full_input_shuffled_noLOH_" + str(n) + "/results/")

for directory in directories:
    filenames = []
    for (_, _, f) in walk(directory):
        filenames += f
        break

    for f in filenames:
        if f.find(".tree") == -1:
            continue
        (patient, sample, __) = f.split("_")
        #print "Reading", f
        tree = Tree(directory + f)
        tree.set_outgroup("Germline")
        balanced = "Balanced"
        testable = "Untestable"
        #print tree
        subpopfreqs = 0
        for node in tree.traverse("levelorder"):
            n = node.name
            if n.find("_") == -1:
                continue
            (id, freq) = n.split("_")
            try:
                freq = float(freq)
            except:
                continue
            subpopfreqs += freq
            childfreqs = 0
            for child in node.children:
                cn = child.name
                print cn
                if cn.find("_") == -1:
                    print "A child node with an odd name:", cn, "in", f
                    continue
                (__, cfreq) = cn.split("_")
                try:
                    cfreq = float(cfreq)
                except:
                    continue
                testable = "Testable"
                childfreqs += cfreq
            if (childfreqs > freq):
                balanced = "Unbalanced"
        print directory + "\t" + f + "\t" + balanced + "\t" + testable + "\t" + str(len(tree)) + "\t" + str(subpopfreqs)
