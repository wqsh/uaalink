#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: get_one_mgf.py
Author: webdev@tencent.com
Date: 2019/11/10 16:11:29
Brief: 
"""

import sys
reload(sys)
sys.setdefaultencoding('utf8')

mgf_path = './dataset4/20190418_BING_K_IN_SOLUTION_20190418175056_HCDFT.mgf'
title = sys.argv[1]

fw = open('one.mgf', 'w')
line_list = []
flag = False
with open(mgf_path) as f:
    for line in f:
        line = line.strip()
        if line.startswith('TITLE='):
            t = line.split('=')[1]
            if t == title:
                flag = True

        if flag:
            line_list.append(line)
        if line.startswith("END IONS"):
            flag = False

fw.write("BEGIN IONS" + '\n')
for line in line_list:
    fw.write(line + '\n')
fw.close()
