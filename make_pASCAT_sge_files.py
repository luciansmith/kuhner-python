d p#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:36:41 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir

for tag in ["combined_avSNPs", "combined_avSNPs_only1M", "combined_avSNPs_only25M"]:

    pASCAT_dir = "pASCAT_input_" + tag + "/"
    pASCAT_out = "pASCAT_input_" + tag + "/"

    if not(path.isdir(pASCAT_out)):
        mkdir(pASCAT_out)

    if not(path.isdir(pASCAT_out + "unconstrained/")):
        mkdir(pASCAT_out + "unconstrained/")

    if not(path.isdir(pASCAT_out + "diploid/")):
        mkdir(pASCAT_out + "diploid/")

    if not(path.isdir(pASCAT_out + "tetraploid/")):
        mkdir(pASCAT_out + "tetraploid/")

    if not(path.isdir(pASCAT_out + "Rout/")):
        mkdir(pASCAT_out + "Rout/")


    filenames = []
    for (_, _, f) in walk(pASCAT_dir):
        filenames += f
        break


    allfile = open(pASCAT_out + "run_all.bat", "w")
    for f in filenames:
        if f.find("logR") == -1:
            continue
        (patient, __) = f.split("_")
        outfname = "run_" + patient + "_" + tag + ".sge"
        outfile = open(pASCAT_out + outfname, "w")
        outfile.write("module load modules modules-init modules-gs gmp/5.0.2 mpfr/latest mpc/0.8.2 gcc/latest R/latest java_jdk/latest\n")
        outfile.write("cd R/pASCAT_output_" + tag + ";\n")
        outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + "' segment.R Rout/" + patient + "_segout.txt;\n")
        outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 0' ascat_on_segments.R Rout/" + patient + "_uncon_out.txt;\n")
        outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*" + " unconstrained;\n")
        outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 2' ascat_on_segments.R Rout/" + patient + "_dip_out.txt;\n")
        outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*" + " diploid;\n")
        outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 4' ascat_on_segments.R Rout/" + patient + "_tet_out.txt;\n")
        outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*" + " tetraploid;\n")
        outfile.close()
        allfile.write("qsub -l mfree=6G run_" + patient + "_" + tag + ".sge\n")

    allfile.close()

