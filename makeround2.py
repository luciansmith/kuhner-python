#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 13:23:46 2017

@author: lpsmith
"""

halted = {}
noex = open("halted.txt", "r")
for l in noex:
    lvec = l.split("_")
    if len(lvec)>=3:
        patient = lvec[0]
        type = lvec[1]
        if not patient in halted:
            halted[patient] = []
        halted[patient].append(type)
        
runall = open("run_all_round3.bat", "w")
for patient in halted:
    fname = "run_" + patient + "_25M_round3.sge"
    runall.write("qsub -l mfree=4G run_" + patient + "_25M_round3.sge\n")
    runpatient = open(fname, "w")
    runpatient.write("module load modules modules-init modules-gs gmp/5.0.2 mpfr/latest mpc/0.8.2 gcc/latest R/latest java_jdk/latest\n")
    runpatient.write("cd R/pASCAT_output_25M_round3\n")
    for type in halted[patient]:
        if type=="uncon":
            dirname = "unconstrained"
            num = "0"
        elif type=="dip":
            dirname = "diploid"
            num="2"
        elif type=="tet":
            dirname = "tetraploid"
            num="4"
        else:
            "Unknown type", type, "for patient", patient
            continue
        runpatient.write("R CMD BATCH '--args " + patient + " " + num + "' ascat_on_segments.R Rout/" + patient + "_" + type + "_out_round2.txt\n")
        runpatient.write("mv " + patient + "*raw* " + patient + "_failed_arrays.txt " + patient + "_fcn_* " + dirname  + "\n")
        

