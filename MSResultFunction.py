#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: MSResultFunction.py
Author: webdev@tencent.com
Date: 2019/11/28 16:11:56
Brief: 
"""

import os
from MSData import CModSite, CPeptide, CSinglePeptide, CLinkPeptide, CSingle_Score_Feature
import time

class CResultFunction:
    def __init__(self, conf):
        self.conf = conf

    def get_taskid(self, conf, need_time=True):
        if len(conf.mgf_file_list) > 0:
            mgf_name = conf.mgf_file_list[0]
        else:
            mgf_name = conf.mgf_folder
        if mgf_name.find('/') >= 0:
            mgf_name = mgf_name.split('/')[-1]
        else:
            mgf_name = mgf_name.split('\\')[-1]
        mgf_name = mgf_name[:mgf_name.rfind('.')]
        if need_time:
            taskid = mgf_name + "-" + str(int(time.time()))
        else:
            taskid = mgf_name
        return taskid

    def find_single_res_file(self, conf):
        mgf_file_list = conf.mgf_file_list
        mgf_set = set()
        for file in mgf_file_list:
            if file.find('/') >= 0:
                file = file.split('/')[-1][:-4]
            else:
                file = file.split('\\')[-1][:-4]
            mgf_set.add(file + "_singlePeptide.txt")
        folder = conf.result_folder
        file_list = []
        for file in os.listdir(folder):
            if file in mgf_set: file_list.append(os.path.join(folder, file))
        return file_list

    def find_res_file(self, conf, cand_flag = False, label_flag = False):
        mgf_file_list = conf.mgf_file_list
        mgf_set = set()
        for file in mgf_file_list:
            if file.find('/') >= 0:
                file = file.split('/')[-1][:-4]
            else:
                file = file.split('\\')[-1][:-4]
            if cand_flag: file += "_cand"
            if not cand_flag and label_flag:
                index = file.rfind('_')
                if index >= 0: file = file[:index]
            if not label_flag: file += ".txt"
            else: file += ".plabel"
            mgf_set.add(file)
        folder = conf.result_folder
        file_list = []
        for file in os.listdir(folder):
            if not label_flag and not file.endswith(".txt"): continue
            if label_flag and not file.endswith(".plabel"): continue
            if file in mgf_set: file_list.append(os.path.join(folder, file))
        return file_list

    def _get_mod_str(self, mods):
        mod_str = ""
        for m in mods:
            mod_str += str(m.site) + "@" + str(m.mod_name) + ";"
        return mod_str

    def _get_str_mod(self, mod_str, mod2mass):
        mod_list = []
        for m in mod_str.split(';'):
            if m.find('@') < 0: continue
            site, name = m.split('@')
            site = int(site)
            one_mod = CModSite(name, site, mod2mass[name])
            mod_list.append(one_mod)
        return mod_list

    def _get_mod2mass(self, mod2mod):
        mod2mass = {}
        for name in mod2mod:
            mod2mass[name] = mod2mod[name].mass
        return mod2mass

    def load_single_psm(self, path):
        psm_list = []
        with open(path, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num <= 1: continue
                title, mass, sq, mod_mass, site, score, delta_score, rank, is_target = line.strip().split('\t')
                mass = float(mass)
                mod_mass = float(mass)
                site = int(site)
                score = float(score)
                delta_score = float(delta_score)
                rank = int(rank)
                decoy_flag = (is_target != "True")
                pro_index = 0
                if decoy_flag: pro_index = -1
                peptide = CSinglePeptide(sq, site, mod_mass, pro_index, score)
                psm_list.append((title, peptide, decoy_flag, rank, delta_score))
        return psm_list


    def load_psm(self, path):
        psm_list = []
        mod2mass = self._get_mod2mass(self.conf.modification.modname2mod)
        with open(path, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num <= 1: continue
                tmp = line.strip().split('\t')
                score = float(tmp[12])
                delta_score = float(tmp[13])
                title = tmp[0]
                sq1 = tmp[2]
                mods1 = tmp[3]
                site1 = int(tmp[4]) - 1
                mass1 = float(tmp[5])
                sq2 = tmp[6]
                mods2 = tmp[7]
                site2 = int(tmp[8]) - 1
                mass2 = float(tmp[9])
                ac_list = tmp[11]
                decoy_flag = True
                for ac in ac_list.split(';'):
                    if not ac.startswith("REV_"): 
                        decoy_flag = False
                        break
                ranker = int(tmp[14])
                score1 = float(tmp[15])
                mean_fra_me1 = float(tmp[16])
                std_fra_me1 = float(tmp[17])
                ion_ratio1 = float(tmp[18])
                inten_ratio1 = float(tmp[19])
                continue_b_score1 = float(tmp[20])
                continue_y_score1 = float(tmp[21])
                continue_score1 = float(tmp[22])
                len1 = int(tmp[23])
                score2 = float(tmp[24])
                mean_fra_me2 = float(tmp[25])
                std_fra_me2 = float(tmp[26])
                ion_ratio2 = float(tmp[27])
                inten_ratio2 = float(tmp[28])
                continue_b_score2 = float(tmp[29])
                continue_y_score2 = float(tmp[30])
                continue_score2 = float(tmp[31])
                len2 = int(tmp[32])
                alpha_score = CSingle_Score_Feature(score1, mean_fra_me1, std_fra_me1, ion_ratio1, inten_ratio1, continue_b_score1, continue_y_score1, continue_score1)
                beta_score = CSingle_Score_Feature(score2, mean_fra_me2, std_fra_me2, ion_ratio2, inten_ratio2, continue_b_score2, continue_y_score2, continue_score2)
                pep1 = CPeptide(sq1, -1, mass1, self._get_str_mod(mods1, mod2mass), [])
                pep2 = CPeptide(sq2, -1, mass2, self._get_str_mod(mods2, mod2mass), [])
                link_peptide = CLinkPeptide(pep1, pep2, site1, site2, score, alpha_score, beta_score)
                psm_list.append((title, link_peptide, decoy_flag, ranker, delta_score))
        return psm_list
