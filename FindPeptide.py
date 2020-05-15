#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: FindPeptide.py
Author: webdev@tencent.com
Date: 2019/05/27 16:05:38
Brief: 
"""

import sys
import cPickle as pickle
import os
from GetPeptideMass import get_mass

query_sq = sys.argv[1].strip()
mass = get_mass(query_sq)
print("Mass", mass)

start_mass = int(mass / 100) * 100
end_mass = start_mass + 100
print("Mass range", start_mass, end_mass)

index_folder = "./index/test_decoy_dataset2/"
f = open(os.path.join(index_folder, '{}_{}.pkl'.format(start_mass, end_mass)), 'rb')
peps = pickle.load(f)

f.close()
f = open('./index/test_decoy_dataset2/protein.pkl', 'rb')
pros = pickle.load(f)
f.close()
f = open('./index/test_decoy_dataset2/modification.pkl', 'rb')
mods = pickle.load(f)
f.close()

for p in peps:
    pro = pros[p.pro_index]
    sq = pro.sq[p.start_pos:p.end_pos]
    if sq == query_sq:
        mod_str = ""
        mmm = p.mods
        for m in mmm:
            one_mod = mods[m.mod_index]
            mod_str += str(m.site) + ":" + one_mod.name + "\t"
        print("Find", sq, p.mass, pro.ac, mod_str)
