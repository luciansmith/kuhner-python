#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:29:58 2019

@author: lpsmith
"""

import os
import ete3

#outdir = "non17_phylogenies/"
#target_pids = [88, 130, 385, 541, 572, 729, 956, 42, 169, 279, 623, 728, 852, 951, 995]

outdir = "all_phylogenies/"
target_pids = [55, 59, 478, 635, 865, 126, 909, 381, 609, 184]

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
    
vertical_margins = {
        88:  25,
        130: 20, 
        385: 20, 
        541: 20, 
        572: 50, 
        729: 55, 
        956: 20, 
        42:  25, 
        169: 28, 
        279: 20, 
        623: 30, 
        728: 20, 
        852: 25, 
        951: 25, 
        995: 25,

        55: 24, 
        59: 18, 
        126: 19, 
        184: 12, 
        381: 16, 
        478: 12, 
        609: 18, 
        635: 20, 
        865: 15, 
        909: 15, 

        286: 18, 
        387: 26, 
        521: 24, 
        672: 33, 
}

for x, path in enumerate(tree_files['path']):
    pid = tree_files['PatientID'][x]
    if pid not in vertical_margins:
        vertical_margins[pid] = 20
#    if pid not in target_pids:
#        continue

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
    tstyle.branch_vertical_margin = vertical_margins[pid]
    tstyle.show_scale = False
    # Don't show the leaf names or print the branch lengths
    tstyle.show_leaf_name = False
    tstyle.show_branch_length = False
    for branch in t:
        if "N" in branch.name:
            branch.delete()
    for branch in t:
        if "_" in branch.name:
            #Print the leaf names in our own font
            name_face = ete3.AttrFace("name", fsize=vertical_margins[pid]+10)
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
    print("\nTree for ", pid, ":")
    print(t)
    #t.show()
    #print("['%%inline' not working for some reason; try later to get prettier trees]")
    t.render(outdir + str(pid) + "_tree.png", tree_style=tstyle)

print("done!")