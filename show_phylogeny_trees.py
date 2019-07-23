#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:29:58 2019

@author: lpsmith
"""

import os
import ete3

outdir = "patty_phylogenies/"

# Ingest the trees.  First, grab all the filenames/ paths of the tree files
tree_files = {'path': [], 'PatientID': []}
for path, _, files in os.walk('phylip_TS_analysis'):
    for file_ in files:
        if file_[-8:] == 'tree.txt':
            tree_files['path'].append(os.path.join(path, file_))
            tree_files['PatientID'].append(file_.split('_')[0])
for x, pid in enumerate(tree_files['PatientID']):
    if "891" in pid:
        pid = 891
    tree_files['PatientID'][x] = int(pid)
    
target_pids = [88, 130, 385, 541, 572, 729, 956, 42, 169, 279, 623, 728, 852, 951, 995]


for x, path in enumerate(tree_files['path']):
    pid = tree_files['PatientID'][x]
    if pid not in target_pids:
        continue

    # And load the tree
    t = ete3.Tree(path)

    # Need to correlate leaf names to DNANums, so grab the leaf names
    #leaves = t.get_leaf_names()

    # Grab just the DNANums of the leaves
    #dnanums = [x.split('_')[0] for x in leaves if 'blood' != x[:5]]

    # Define the tree style
    tstyle = ete3.TreeStyle()
    # Make it a rectangular tree
    tstyle.mode = 'r'
    tstyle.branch_vertical_margin = 20
    tstyle.show_scale = False
    # Show the leaf names, and the branch lenghts
    tstyle.show_leaf_name = False
    tstyle.show_branch_length = False
    for branch in t:
        if "_" in branch.name:
            name_face = ete3.AttrFace("name", fsize=30)
            branch.add_face(name_face, column=0, position="branch-right")
#        if "blood" in branch.name:
#            branch.name = "blood " + str(pid)
#        elif "_" in branch.name:
            branch.name = branch.name.split("_")[0]
    # Stretch in x
    tstyle.scale = 1000

    #print("Before displaying anything")
    # Before we display this tree, let's output the DNANum
    #display(Markdown('## Patient ' + str(pid)))

    #Somehow, '%%inline' isn't working any more; revert to just printing the basic tree with zero information (sad)
    print(t)
    #t.show()
    #print("['%%inline' not working for some reason; try later to get prettier trees]")
    t.render(outdir + str(pid) + "_tree.png", tree_style=tstyle)

print("done!")