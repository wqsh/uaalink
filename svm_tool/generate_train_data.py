#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: generate_train_data.py
Author: webdev@tencent.com
Date: 2019/11/29 16:11:23
Brief: 
"""

import sys

all_data_path = "./input_data/all_data.txt"
fw = open("./input_data/train_data.txt", 'w')

topk, last_topk = 300, 100

all_lines = []
with open(all_data_path) as f:
    for line in f:
        all_lines.append(line.strip())

for i in range(topk):
    if i >= len(all_lines): break
    if len(all_lines[i]) == 0: continue
    #if all_lines[i][0] == '1':
    fw.write(all_lines[i] + '\n')

for i in range(last_topk):
    i2 = len(all_lines) - 1 - i
    if i2 < 0: break
    if len(all_lines[i2]) == 0: continue
    if all_lines[i2][0] == '-':
        fw.write(all_lines[i2] + '\n')

