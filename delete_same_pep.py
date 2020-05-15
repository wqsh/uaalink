#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: delete_same_pep.py
Author: webdev@tencent.com
Date: 2019/12/07 19:12:33
Brief: 
"""

import sys, os


score = 0.0
if len(sys.argv) > 1: score = float(sys.argv[1])

config_path = "./configure.txt"
if len(sys.argv) > 2: config_path = sys.argv[2]
mgf_file_list = []
folder = "./result"
with open(config_path) as f:
    for line in f:
        line = line.strip()
        index = line.find('#')
        if index >= 0: line = line[:index].strip()
        if line.startswith('PATH_MS2'):
            value = line.split('=')[1]
            if os.path.isdir(value):
                for p in os.listdir(value):
                    if not p.endswith(".mgf"): continue
                    mgf_file_list.append(os.path.join(value, p))
            else:
                mgf_file_list = value.split(';')
        elif line.startswith('PATH_RESULT_EXPORT'):
            folder = line.split('=')[1]

print("Mgf file list", mgf_file_list)
print("Result folder", folder)

file_list = []
plabel_file_list = []
mgf_set = set()
for file in mgf_file_list:
    if file.find('/') >= 0:
        file = file.split('/')[-1][:-4]
    else:
        file = file.split('\\')[-1][:-4]
    file1 = file + ".txt"
    if file.endswith("_HCDFT"):
        file2 = file[:-6] + ".plabel"
    else:
        file2 = file + ".plabel"
    mgf_set.add(file1)
    mgf_set.add(file2)

for file in os.listdir(folder):
    if not file.endswith(".txt"): continue
    if file in mgf_set: file_list.append(os.path.join(folder, file))

for file in os.listdir(folder):
    if not file.endswith(".plabel"): continue
    if file in mgf_set: plabel_file_list.append(os.path.join(folder, file))

print("plabel file:", plabel_file_list)

pep_score = {}
pep_line = {}
t_p_line = {}
head_line = ""
print("File", file_list)
title_score = {}
title_line = {}
title_plabel_line = {}
pep_plabel_line = {}

for file in file_list:
    with open(file) as f:
        line_num = 0
        for line in f:
            line = line.strip()
            line_num += 1
            if line_num <= 1: 
                head_line = line
                continue
            tmp = line.strip().split('\t')
            title = tmp[0]
            sq2 = tmp[2] + "@" + tmp[3] + "@" + str(tmp[4])
            ac_list = tmp[11].split(';')
            is_target = False
            for ac in ac_list:
                if not ac.startswith("REV_"):
                    is_target = True
                    break
            if not is_target: continue
            s = float(tmp[12])
            if s <= score: break
            title_score[title] = s
            title_line[title] = line
            if sq2 not in pep_score or pep_score[sq2] < s:
                pep_score[sq2] = s
                pep_line[sq2] = line
        print("Line", file, line_num)

for sq in pep_line:
    line = pep_line[sq]
    tmp = line.strip().split('\t')
    title = tmp[0]
    t_p_line[title] = 1

modificaton_lines = []
for file in plabel_file_list:
    with open(file) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if lines[i].startswith("name="):
            title = lines[i].strip().split('=')[1].strip()
            if title in title_line:
                line = lines[i+1].strip()
                title_plabel_line[title] = line
            if title in t_p_line:
                line = lines[i+1].strip()
                pep_plabel_line[title] = line
        elif lines[i].startswith("[Modification]"):
            modificaton_lines = []
            for j in range(i+1, len(lines)):
                ind = lines[j].find('=')
                if ind < 0: break
                modificaton_lines.append(lines[j].strip())

sorted_res = sorted(title_score.items(), key=lambda d: d[1], reverse=True)
with open('./res_psm.txt', 'w') as f:
    f.write(head_line + '\n')
    for (r, s) in sorted_res:
        one_line = title_line[r]
        f.write(one_line + '\n')

with open('./res_psm.plabel', 'w') as f:
    f.write("[FilePath]\n")
    f.write("File_Path=" + ";".join(mgf_file_list) + "\n")
    f.write("[Modification]\n")
    for mod_line in modificaton_lines:
        f.write(mod_line + '\n')
    f.write("[xlink]\n")
    f.write("xlink=UUA\n")
    f.write("[Total]\n")
    f.write("total=" + str(len(sorted_res)) + '\n')
    spec_id = 0
    for (r, s) in sorted_res:
        spec_id += 1
        one_line = title_plabel_line[r]
        f.write("[Spectrum" + str(spec_id) + "]\n")
        f.write("name=" + r + "\n")
        f.write(one_line + "\n")



sorted_res = sorted(pep_score.items(), key=lambda d: d[1], reverse=True)
with open('./res.txt', 'w') as f:
    f.write(head_line + '\n')
    for (r, s) in sorted_res:
        one_line = pep_line[r]
        f.write(one_line + '\n')

with open('./res.plabel', 'w') as f:
    f.write("[FilePath]\n")
    f.write("File_Path=" + ";".join(mgf_file_list) + "\n")
    f.write("[Modification]\n")
    for mod_line in modificaton_lines:
        f.write(mod_line + '\n')
    f.write("[xlink]\n")
    f.write("xlink=UUA\n")
    f.write("[Total]\n")
    f.write("total=" + str(len(sorted_res)) + '\n')
    spec_id = 0
    for (r, s) in sorted_res:
        spec_id += 1
        one_line_tmp = pep_line[r]
        tmp = one_line_tmp.strip().split('\t')
        title = tmp[0]
        f.write("[Spectrum" + str(spec_id) + "]\n")
        f.write("name=" + title + "\n")
        one_line = pep_plabel_line[title]
        f.write(one_line + "\n")

