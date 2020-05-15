#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: MSData.py
Author: webdev@tencent.com
Date: 2019/05/21 20:05:52
Brief: 
"""

import sys
#reload(sys)
import re

class CElement:
    def __init__(self, element_file):
        self.element2mass = {}
        with open(element_file, 'r') as f:
            for line in f:
                if line.find('=') < 0: continue
                if line.strip() == "" or line[0] == '@': continue
                line = line.strip()
                line = line[line.find('=')+1:]
                tmp = line.split('|')
                element_name = tmp[0]
                element_mass = [float(v) for v in tmp[1][:-1].split(',')]
                element_percentage = [float(v) for v in tmp[2][:-1].split(',')]
                assert(len(element_mass) == len(element_percentage))
                index = element_percentage.index(max(element_percentage))
                self.element2mass[element_name] = element_mass[index]

class CAA:
    def __init__(self, aa_file, elem2mass, aa_num):
        self.aa2mass = [0 for i in range(aa_num)]
        with open(aa_file, 'r') as f:
            for line in f:
                if line.find('=') < 0: continue
                if line.strip() == "" or line[0] == '@': continue
                name, value = line.strip().split('=')
                aa, elem_str = value[:-1].split('|')
                self.aa2mass[int(ord(aa[0]) - ord('A'))] = self.compute_mass_from_element(elem_str, elem2mass)

    def compute_mass_from_element(self, elem_str, elem2mass):
        tmp = re.split('\(|\)', elem_str)
        mass = 0.0
        for i in range(len(tmp)):
            if i % 2 != 0: continue
            if tmp[i].strip() == "": continue
            if tmp[i] not in elem2mass: continue
            mass += (float(elem2mass[tmp[i]]) * int(tmp[i+1]))
        return mass


class CMod:
    def __init__(self, name, mass, type, aa):
        self.name = name # modification name
        self.mass = mass # modification mass
        self.type = type # modification type (NORMAL, PEP_N, PEP_C, PRO_N, PRO_C)
        self.aa = aa # modification can be occured in which amino acids

class CModification:
    def __init__(self, mod_file):
        self.modname2mod = {}
        with open(mod_file, 'r') as f:
            for line in f:
                if line.find('=') < 0: continue
                if line.strip() == "" or line[0] == '@': continue
                if line.startswith("name"): continue
                name, value = line.strip().split('=')
                tmp = value.split(' ')
                mass = float(tmp[2])
                aa = tmp[0]
                type = tmp[1]
                self.modname2mod[name] = CMod(name, mass, type, aa)

class CModSite:
    def __init__(self, mod_name, site, mass):
        self.mod_name = mod_name # the modification name
        self.site = site # the modification occured in the index of peptide sequence
        self.mass = mass # the mass of modification

class CPeptide:
    def __init__(self, sq, pro_index, mass, mods=[], pro_index_list=[]):
        self.sq = sq # peptide sequence
        self.pro_index = pro_index # protein index
        self.mass = mass # mass of peptide sequence
        self.mods = mods # modifications in peptide sequence
        self.pro_index_list = pro_index_list

    def print_peptide(self):
        mod_str = ""
        for m in self.mods:
            mod_str += str(m.site) + ":" + str(m.mod_name) + "\t"
        print("Sq:{0}, Protein_index:{1}, Mass:{2}, Mods:{3}".format(self.sq, self.pro_index, self.mass, mod_str))

class CModSitePKL:
    def __init__(self, mod_index, site):
        self.mod_index = mod_index
        self.site = site

class CPeptidePKL: # In order to save storage space, use start and end position (int) to replace the sequence.
    def __init__(self, pro_index, start_pos, end_pos, mass, mods=[], gdm=0.0):
        self.pro_index = pro_index
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.mass = mass
        self.mods = mods
        self.gdm = gdm
        self.pro_index_list = [] # the peptides with same sequence and modifications (sites) should be combined in multiply proteins


class CProtein:
    def __init__(self, ac, de, sq):
        self.ac = ac
        self.de = de
        self.sq = sq


class CPeak:
    def __init__(self, mz, inten):
        self.mz = mz
        self.inten = inten

class CSpectrum:
    def __init__(self, title, mass, charge, peaks=[]):
        # the mass of precursor ion is added by (H2O+P)
        self.title = title
        self.mass = mass
        self.charge = charge
        self.peaks = peaks

class CSinglePeptide:
    def __init__(self, beta_peptide_sq, beta_peptide_site, beta_peptide_modmass, pro_index, score=0.0):
        self.peptide_sq = beta_peptide_sq
        self.peptide_site = beta_peptide_site
        self.peptide_modmass = beta_peptide_modmass
        self.pro_index = pro_index
        self.score = score
        

class CLinkPeptide:
    def __init__(self, alpha_peptide, beta_peptide, alpha_site, beta_site, score=0.0, alpha_score=None, beta_score=None):
        #alpha_score and beta_score are class objects (see "CSingle_Score_Feature")
        self.alpha_peptide = alpha_peptide
        self.beta_peptide = beta_peptide
        self.alpha_site = alpha_site
        self.beta_site = beta_site
        self.score = score
        self.alpha_score = alpha_score
        self.beta_score = beta_score

class CLink_Score_Feature:
    def __init__(self, pre_me, alpha_score, beta_score):
        #pre_me: precursor mass error
        #alpha_score: score feature of alpha peptide
        #beta_score: score feature of beta peptide
        self.pre_me = pre_me
        self.alpha_score = alpha_score
        self.beta_score = beta_score

class CSingle_Score_Feature:
    def __init__(self, score, mean_fra_me, std_fra_me, ion_ratio, inten_ratio, continue_b_score, continue_y_score, continue_score):
        #score: use pNovo score method
        #mean_fra_me: mean of fragment mass error
        #std_fra_me: std of fragment mass error
        #ion_ratio: percentage of matched ions
        #inten_ratio: percentage of matched peak intensities
        #continue_score: score of continue matched peaks
        self.score = score
        self.mean_fra_me = mean_fra_me
        self.std_fra_me = std_fra_me
        self.ion_ratio = ion_ratio
        self.inten_ratio = inten_ratio
        self.continue_b_score = continue_b_score
        self.continue_y_score = continue_y_score
        self.continue_score = continue_score

    def get_string(self, split_flag = '\t', fid_flag = False, fid_start = 0):
        if fid_flag:
            return str(fid_start) + ":" + str(self.score) + split_flag + str(fid_start+1) + ":" + str(self.mean_fra_me) + split_flag \
            + str(fid_start+2) + ":" + str(self.std_fra_me) + split_flag + str(fid_start+3) + ":" + str(self.ion_ratio) + split_flag \
            + str(fid_start+4) + ":" + str(self.inten_ratio) + split_flag + str(fid_start+5) + ":" + str(self.continue_b_score) + split_flag \
            + str(fid_start+6) + ":" + str(self.continue_y_score) + split_flag + str(fid_start+7) + ":" + str(self.continue_score)
        return str(self.score) + split_flag + str(self.mean_fra_me) + split_flag + str(self.std_fra_me) + split_flag \
        + str(self.ion_ratio) + split_flag + str(self.inten_ratio) + split_flag \
        + str(self.continue_b_score) + split_flag + str(self.continue_y_score) + split_flag \
        + str(self.continue_score)

    def _print(self):
        print("+++++++++print single score+++++++++++++")
        print("Score:", self.score)
        print("Mean fragment mass error:", self.mean_fra_me)
        print("Std fragment mass error:", self.std_fra_me)
        print("Ion ratio", self.ion_ratio)
        print("Inten ratio", self.inten_ratio)
        print("B continue score", self.continue_b_score)
        print("Y continue score", self.continue_y_score)
        print("B\Y continue score", self.continue_score)

class CPSM:
    def __init__(self, spec, peptide, peptide_list = [], decoy_flag=False):
        self.spec = spec
        self.peptide = peptide
        self.decoy_flag = decoy_flag
        if len(peptide_list) == 0:
            peptide_list.append(peptide)
        self.peptide_list = peptide_list
