#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 15:29:33 2019

@author: lpsmith
"""

ntips = 50
x2tips = ntips*2

runall = open("run_all.bat", "w")
for i in ["0"]:
    for j in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
        for k in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
            ijk = i+j+k
            simout = open("sim_thinned_" + ijk + ".sge", "w")
            runall.write("qsub -l mfree=4G sim_thinned_" + ijk + ".sge\n")
            simout.write("module load modules modules-init modules-gs gcc/latest python/latest numpy/1.16.0 scipy/1.2.0 glibc/2.14\n")
            simout.write("cd tsinfer_2xr/" + str(ntips) + "tip_thinned/\n")
            simout.write("rm reccount\n")
            simout.write("rm migcount\n")
            simout.write("rm totreccount\n")
            simout.write("rm RFresults\n")
            simout.write("# ./ms " + str(x2tips) + " 1 -T -t 550.000000 -r 19.800000 100000 | tail -n +4 | grep -v // > truetree" + ijk + "\n")
            simout.write("python3 truncate_ms_tree.py truetree" + ijk + " treefile" + ijk + "\n")
            simout.write("python thinned.py treefile" + ijk + " thinnedtree" + ijk + " " + str(ntips) + "\n")
            simout.write("./seq-gen -mF84 -l 100000 -s 0.005500 -p 1000 thinnedtree" + ijk + " > dnafile" + ijk + "\n")
            simout.write("python3 dna_to_bool.py dnafile" + ijk + " translated_seqgen_dna" + ijk + "\n")
            simout.write("python3 parse_seqgen.py translated_seqgen_dna" + ijk + " treedump" + ijk + " inferredtree" + ijk + "\n")
            #simout.write("mv timefile timefile" + ijk + "\n")
            simout.write("# compare converted original tree to inferred tree\n")
            simout.write("python ts_argcompare.py treefile" + ijk + " inferredtree" + ijk + "\n")
            simout.write("rm seqgen_dna translated_seqgen_dna treedump" + ijk + "\n")
            simout.write("rm reccount migcount totreccount\n")
            simout.close()
runall.close()

initout = open("sim_init_trees.sge", "w")
initout.write("module load modules modules-init modules-gs gcc/latest python/latest numpy/1.16.0 scipy/1.2.0 glibc/2.14\n")
initout.write("cd tsinfer_2xr/" + str(ntips) + "tip_thinned/\n")
initout.write("for i in 0\n")
initout.write("do\n")
initout.write("for j in 0 1 2 3 4 5 6 7 8 9\n")
initout.write("do\n")
initout.write("for k in 0 1 2 3 4 5 6 7 8 9\n")
initout.write("do\n")
initout.write("  ./ms " + str(x2tips) + " 1 -T -t 550.000000 -r 19.800000 100000 | tail -n +4 | grep -v // > truetree$i$j$k\n")
initout.write("done\n")
initout.write("done\n")
initout.write("done\n")
initout.close()