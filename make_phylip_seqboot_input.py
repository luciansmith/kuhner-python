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

import numpy
import math
import matplotlib.pyplot as plt
import csv

import lucianSNPLibrary as lsl


mutation_file = "snv_plus_indels.twoPlus.20181030.csv"
ploidy_file = "calling_evidence_odds.tsv"
seqboot_dir = "phylip_seqboot_lowVAF_analysis/"
phylip_dir = "phylip_input_lowVAF/"
samplefile = "20181031_SampleCodeLucianTrees_agesDiffRounding_withPilot.txt"

if not path.isdir(seqboot_dir):
    mkdir(seqboot_dir)


seqboot_batch = open(seqboot_dir + "run_seqboot.bat", "w")

def getOutgroupFrom(patient):
    orig_in = open(phylip_dir + patient + "_infile.txt", "r")
    orig_in.readline()
    orig_in.readline()
    orig_in.readline()
    line = orig_in.readline()
    line = line.rstrip()
    return line
    

flist = []
for (_, _, f) in walk(seqboot_dir):
    flist += f
    
for f in flist:
    if "phylip_in" not in f:
        continue
    patient = f.split('_')[0]
    seqboot_in = open(seqboot_dir + patient + "_seqboot_in.txt", "w")
    seqboot_batch.write("./seqboot < " + patient + "_seqboot_in.txt\n")
    seqboot_in.write(patient + "_phylip_in.txt\n")
    seqboot_in.write("I\n")
    seqboot_in.write("R\n")
    seqboot_in.write("1000\n")
    seqboot_in.write("Y\n")
    seqboot_in.write("1001\n")
    seqboot_batch.write("mv outfile " + patient + "_dnapars_data.txt\n")
    seqboot_batch.write("./dnapars < " + patient + "_dnapars_in.txt\n")
    dnapars_in = open(seqboot_dir + patient + "_dnapars_in.txt", "w")
    dnapars_in.write(patient + "_dnapars_data.txt\n")
    dnapars_in.write("I\n")
    dnapars_in.write("O\n")
    dnapars_in.write(getOutgroupFrom(patient) + "\n")
    dnapars_in.write("M\n")
    dnapars_in.write("D\n")
    dnapars_in.write("1000\n")
    dnapars_in.write("1001\n")
    dnapars_in.write("1\n")
    dnapars_in.write("Y\n")
    seqboot_batch.write("mv outfile " + patient + "_dnapars_out.txt\n")
    seqboot_batch.write("mv outtree " + patient + "_dnapars_trees.txt\n")
    seqboot_batch.write("./consense < " + patient + "_consense_in.txt\n")
    consense_in = open(seqboot_dir + patient + "_consense_in.txt", "w")
    consense_in.write(patient + "_dnapars_trees.txt\n")
    consense_in.write("Y\n")
    seqboot_batch.write("mv outfile " + patient + "_consense_out.txt\n")
    seqboot_batch.write("mv outtree " + patient + "_consense_tree.txt\n")
    seqboot_in.close()
    dnapars_in.close()
    consense_in.close()


seqboot_batch.close()