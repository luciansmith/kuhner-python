#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 11:38:28 2018

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy

import lucianSNPLibrary as lsl


mutcount_file = "tri_out/mutcount_all.tsv"

for line in open(mutcount_file, "r"):
    