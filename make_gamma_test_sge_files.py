#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:36:41 2017

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
from shutil import copytree

gamma_dir = "gamma_template/"
gamma_out = "gamma_test_nonans_pilot/"
if not(path.isdir(gamma_out)):
    mkdir(gamma_out)

allfile = open(gamma_out + "run_all.bat", "w")
makelinks = open(gamma_out + "make_links.bat", "w")

onlysomepatients = True
eightplus = False
#somepatients = ["17", "42", "55", "59", "74", "43", "163", "396", "1047"]
#somepatients = ["184"]
#somepatients = ["42", "43", "55", "59", "74"]
somepatients = ["1005", "222", "391", "422", "43", "551", "575", "59", "611", "619", "639", "672", "686", "728", "88", "915", "954"]

infiles = []
patients = []
for (_, _, f) in walk(gamma_dir):
    infiles += f
for file in infiles:
    if file.find("_logR.txt") != -1:
        patients.append(file.split("_")[0])


#for gamma in [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 2000, 2500, 3000]:
#for patient in firstpatients:
#for gamma in [3000, 1000, 400, 250, 100]:
for gamma in [1000, 400, 250, 100, 2000, 200, 600, 800, 150, 300, 350, 450, 500, 700, 900, 1200, 1400, 1600, 2500]:
#for gamma in [450]:
    pASCAT_dir = "pASCAT_input_g" + str(gamma) + "/"

    if not(path.isdir(gamma_out + pASCAT_dir)):
        mkdir(gamma_out + pASCAT_dir)

#    copytree(gamma_dir, gamma_out + pASCAT_dir)

    if not(path.isdir(gamma_out + pASCAT_dir + "diploid/")):
        mkdir(gamma_out + pASCAT_dir + "diploid/")

    if not(path.isdir(gamma_out + pASCAT_dir + "tetraploid/")):
        mkdir(gamma_out + pASCAT_dir + "tetraploid/")

    if eightplus:
        if not(path.isdir(gamma_out + pASCAT_dir + "eight/")):
            mkdir(gamma_out + pASCAT_dir + "eight/")

        if not(path.isdir(gamma_out + pASCAT_dir + "eight_high/")):
            mkdir(gamma_out + pASCAT_dir + "eight_high/")

    if not(path.isdir(gamma_out + pASCAT_dir + "Rout/")):
        mkdir(gamma_out + pASCAT_dir + "Rout/")


    makelinks.write("cd " + pASCAT_dir + "\n")
    makelinks.write("cp ../../" + gamma_dir + "*.R .\n")
    for patient in patients:
        if onlysomepatients and patient not in somepatients:
            continue
#    for patient in patients: #["266", "303", "360"]:
        outfname = "run_" + patient + "_g" + str(gamma) + ".sge"
        outfile = open(gamma_out + pASCAT_dir + outfname, "w")
        outfile.write("module load modules modules-init modules-gs gmp/latest mpfr/latest mpc/latest gcc/4.9.1 R/latest java_jdk/latest\n")
        outfile.write("cd R/" + gamma_out + pASCAT_dir + ";\n")
        outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " " + str(gamma) + "' segment.R Rout/" + patient + "_segout.txt;\n")
        #outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 0' ascat_on_segments.R Rout/" + patient + "_uncon_out.txt;\n")
        #outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*txt" + " unconstrained;\n")
        if not eightplus:
            outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 2' ascat_on_segments.R Rout/" + patient + "_dip_out.txt;\n")
            outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*" + " diploid;\n")
            outfile.write("rm diploid/*" + patient + "_*.RData\n")
            outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 4' ascat_on_segments.R Rout/" + patient + "_tet_out.txt;\n")
            outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*" + " tetraploid;\n")
            outfile.write("rm tetraploid/*" + patient + "_*.RData\n")
            outfile.write("rm *" + patient + "_*.RData\n")
        else:
            outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 8' ascat_on_segments.R Rout/" + patient + "_eight_out.txt;\n")
            outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*" + " eight;\n")
            outfile.write("rm eight/*" + patient + "_*.RData\n")
            outfile.write("R CMD BATCH '--no-save --no-restore --args " + patient + " 8plus' ascat_on_segments.R Rout/" + patient + "_eight_high_out.txt;\n")
            outfile.write("mv " + patient + "*raw* "+ patient + "_failed_arrays.txt " + patient + "_fcn_*" + " eight_high;\n")
            outfile.write("rm eight_high/*" + patient + "_*.RData\n")
        outfile.close()
        allfile.write("qsub -l mfree=6G -l h_rt=24:0:0 pASCAT_input_g" + str(gamma) + "/run_" + patient + "_g" + str(gamma) + ".sge\n")

        makelinks.write("ln -s ../../" + gamma_dir + patient + "_BAF.txt " + patient + "_BAF.txt\n")
        makelinks.write("ln -s ../../" + gamma_dir + patient + "_Normal_BAF.txt " + patient + "_Normal_BAF.txt\n")
        makelinks.write("ln -s ../../" + gamma_dir + patient + "_logR.txt " + patient + "_logR.txt\n")
    #    makelinks.write("mkdir Rout\n")
    #    makelinks.write("mkdir unconstrained\n")
    #    makelinks.write("mkdir diploid\n")
    #    makelinks.write("mkdir tetraploid\n")
    makelinks.write("cd ..\n")

allfile.close()
makelinks.close()

