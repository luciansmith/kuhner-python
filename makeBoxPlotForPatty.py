#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:32:32 2019

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import readlink
from os import mkdir
from os.path import isfile
from copy import deepcopy
from ete3 import Tree
import seaborn as sns
import pandas as pd

import numpy
import math
import matplotlib.pyplot as plt
import csv


drivers_file = "20191003_CountMutatedPutativeEADriversByPatient.txt"
data = pd.read_csv(drivers_file, sep="\t")

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
PS_colors = {
        'NP': 'blue',
        'Prog': 'orangered',
    }

sns.boxplot(x=data.Prog, y=data.Count,
              ax=ax, color="white");
sns.swarmplot(x=data.Prog, y=data.Count, ax=ax,
             palette=PS_colors)#, edgecolor='grey', alpha=0.99, linewidth=.8);
plt.savefig("FIGX_EAdrivers_PNP_count.png")
plt.show()
