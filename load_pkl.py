#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: load_pkl.py
Author: webdev@tencent.com
Date: 2019/05/23 21:05:21
Brief: 
"""

import cPickle as pickle
import sys
import os

folder = "./index/test"
pkl_list = []
data = []
for p in os.listdir(folder):
    if p == "protein.pkl" and p == "modification.pkl":
        continue
    if p.endswith("_ind.pkl"): continue

    path = os.path.join(folder, p)
    print("P", path)
    f = open(path, 'rb')
    data = pickle.load(f)
    print("Len", len(data))
    data = []
    del data
    f.close()
