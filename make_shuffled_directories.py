#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 12:34:23 2017

@author: lpsmith
"""

#shuffleall = open("run_shuffle_all.bat", "w")
for r in range(100):
#    print "cp -r expands_full_input expands_full_input_shuffled_noLOH_" + str(r)
#    print "perl -p -i -e 's/full_input/full_input_shuffled_noLOH_" + str(r) + "/' expands_full_input_shuffled_noLOH_" + str(r) + "/*sge"

    print "cd expands_full_input_shuffled_noLOH_" + str(r)
    print "mkdir results"
    print "./runall.bat"
    print "cd .."

#    shuffle = open("run_validate_shuffled_noLOH_" + str(r) + ".sge", "w")
#    shuffle.write("cd R/\n")
#    shuffle.write("python validate_expands_results.py expands_full_input_shuffled_" + str(r) + "/results\n")
#    shuffleall.write("qsub run_validate_shuffled_" + str(r) + ".sge\n")