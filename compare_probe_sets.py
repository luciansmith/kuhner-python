
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:11:30 2017

@author: lpsmith
"""

from __future__ import division
import lucianSNPLibrary as lsl
from os import walk
#from os import mkdir
#from os import path
#import string

m1labels, m1rev_labels = lsl.getSNPLabels(True, True)
m25labels, m25rev_labels = lsl.getSNPLabels(False, True)

combined_labels = {}
renamed_labels = {}
repositioned_labels = []
missing_labels = []
nzerozero = 0
noneone = 0
nzero2one = 0
none2zero = 0

transferred = 0
new_snps = 0
missing_zero = 0
missing_gone = 0

for label in m1labels:
    m1l = m1labels[label]
    if label in m25labels:
        m25l = m25labels[label]
        if m1labels[label] == m25labels[label]:
            combined_labels[label] = m1labels[label]
            transferred += 1
        else:
#            print label, "is different in 1m vs. 2.5m"
#            print "\t1m:", m1labels[label]
#            print"\t2.5m:", m25labels[label]
            if m1l[0] == "0" or m1l[1]=="0":
                if m25l[0] == "0" or m25l[0] == "0":
                    nzerozero += 1
                else:
                    nzero2one += 1
            else:
                if m25l[0] == "0" or m25l[0] == "0":
                    none2zero += 1
                else:
                    noneone += 1
            combined_labels[label] = m25labels[label]
            repositioned_labels.append(label)
    else:
        if m1labels[label] in m25rev_labels:
            combined_labels[label] = m1labels[label]
            renamed_labels[label] = m25rev_labels[m1labels[label]]
            new_snps += 1
        else:
            #print label, "not found in 2.5m:", m1labels[label]
            missing_labels.append(label)
            if m1l[0] == "0" or m1l[1]=="0":
                missing_zero += 1
            else:
                missing_gone += 1
        
print "0 -> 0:", nzerozero
print "0 -> 1:", nzero2one
print "1 -> 0:", none2zero
print "1 -> 1:", noneone
print "Successful transfer", transferred
print "New SNP at same position", new_snps
print "Missing, but was zero", missing_zero
print "Missing; gone entirely", missing_gone