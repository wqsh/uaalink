#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: MSTool.py
Author: webdev@tencent.com
Date: 2019/05/27 18:05:15
Brief: 
"""

import sys
#reload(sys)
try:
    import cPickle as pickle
except:
    import pickle
import os
import operator
from MSFunction import CPickleFunction
from MSData import CPeptide, CModSite, CSinglePeptide, CLinkPeptide, CSingle_Score_Feature
import time
import numpy as np
import math

def update_peptide_list(pep_list, pro_list):
    #update peptide sequence when the first sequence is decoy and delta_score is 0.0
    if len(pep_list) <= 1: return pep_list
    ac = pro_list[pep_list[0].alpha_peptide.pro_index].ac
    if not ac.startswith("REV_"): return pep_list
    if pep_list[0].score > pep_list[1].score: return pep_list
    for i in range(1, len(pep_list)):
        if pep_list[i].score < pep_list[0].score: break
        new_ac = pro_list[pep_list[i].alpha_peptide.pro_index].ac
        if not new_ac.startswith("REV_"): break # target protein
    if i < len(pep_list) and pep_list[i].score >= pep_list[0].score-0.0001:
        tmp = pep_list[0]
        pep_list[0] = pep_list[i]
        pep_list[i] = tmp
    return pep_list

def link_site_valid(C_str, N_str, sq, site):
    if site == len(sq)-1 and sq[site] in C_str:
        return False
    if site == 0 and sq[site] in N_str:
        return False
    return True

def global_match_for_single_peptide(input_data):
    spec, beta_pep, conf = input_data
    beta_site = beta_pep.sq.find(conf.beta_uua)
    peak_index, all_peak_inten = create_peak_index(spec.peaks, conf.multi)  # create the peak index
    beta_score, beta_match_info, beta_score_f = match_beta_PSM(spec, peak_index, all_peak_inten, beta_pep, beta_site,
                                                               conf)
    return CSinglePeptide(beta_pep.sq, beta_site, (spec.mass - conf.ATOM_MASS_P - beta_pep.mass), beta_pep.pro_index, beta_score)


def global_match_for_one_link_peptide(input_data):
    spectrum_list, pkl_file, index_pkl_file, beta_pep, conf = input_data
    multi = conf.multi
    cand_num = conf.cand_num # number of candidate peptides for each spectrum, (e.g., 10000)
    print("[Info]#Spectra for {0} is {1}".format(pkl_file, len(spectrum_list)))
    #if not pkl_file.endswith('1500_1600.pkl'): return []
    f = open(pkl_file, 'rb')
    peptide_list = pickle.load(f)
    f.close()
    if len(peptide_list) == 0: return []
    f = open(index_pkl_file, 'rb')
    peptide_index = pickle.load(f)
    f.close()
    f = open(os.path.join(conf.index_out_folder, conf.PRO_PKL_FILE), 'rb')
    protein_list = pickle.load(f)
    f.close()
    f = open(os.path.join(conf.index_out_folder, conf.MOD_PKL_FILE), 'rb')
    modification_list = pickle.load(f)
    f.close()
    beta_mass = beta_pep.mass
    beta_site = beta_pep.sq.find(conf.beta_uua)
    pickleFunc = CPickleFunction()
    pep_start_mass, pep_end_mass = pickleFunc.get_start_end_mass(pkl_file)
    link_aa_dict = conf.link_aa_dict
    res = []
    C_str, N_str = "", ""
    for i in range(len(conf.enzyme_aa)):
        if conf.enzyme_flag[i] == "N":
            N_str += conf.enzyme_aa[i]
        elif conf.enzyme_flag[i] == "C":
            C_str += conf.enzyme_aa[i]
    for spec in spectrum_list:
        start_time = time.time()
        peak_index, all_peak_inten = create_peak_index(spec.peaks, multi) # create the peak index
        end_time = time.time()
        print("[Debug]Time create peak index for spectrum {0} is {1}".format(spec.title, (end_time - start_time)))
        start_time = time.time()
        alpha_mass = spec.mass - beta_mass - conf.ATOM_MASS_P
        if conf.is_ppm_pre: tol_mass = spec.mass * conf.tol_pre
        else: tol_mass = conf.tol_pre
        start_mass =  alpha_mass - tol_mass
        end_mass = alpha_mass + tol_mass
        if start_mass > peptide_list[-1].mass: 
            res.append([])
            continue
        if start_mass < 0: start_mass = 0.0
        if end_mass > peptide_list[-1].mass: end_mass = peptide_list[-1].mass
        s_start_ind = peptide_index[int((start_mass - pep_start_mass) * multi)]
        while s_start_ind < len(peptide_list) and peptide_list[s_start_ind].mass < start_mass:
            s_start_ind += 1
        s_end_ind = peptide_index[int((end_mass - pep_start_mass) * multi)]
        while s_end_ind < len(peptide_list) and peptide_list[s_end_ind].mass <= end_mass:
            s_end_ind += 1
        candidate_peptide_list = peptide_list[s_start_ind: s_end_ind]
        if len(candidate_peptide_list) == 0: 
            res.append([])
            continue
        ori_size = len(candidate_peptide_list)
        candidate_peptide_list = candidate_peptide_list[:cand_num]
        new_size = len(candidate_peptide_list)
        end_time = time.time()
        print("[Debug]Time getting candidate peps for spectrum {0} is {1}".format(spec.title, (end_time - start_time)))
        print("[Debug]#Candidate peptides for spectrum {0} is from {1} to {2}".format(spec.title, ori_size, new_size))
        start_time = time.time()
        link_peptide_list = []
        beta_score, beta_match_info, beta_score_f = match_beta_PSM(spec, peak_index, all_peak_inten, beta_pep, beta_site, conf)
        if conf.is_filter_beta and len(beta_match_info) < conf.beta_filter_v: 
            res.append([])
            continue
        for p in candidate_peptide_list:
            pro = protein_list[p.pro_index]
            sq = pro.sq[p.start_pos:p.end_pos]
            mods = p.mods
            alpha_mods = []
            mod_site_dict = set()
            for m in mods:
                one_mod = modification_list[m.mod_index]
                mod_site_dict.add(m.site-1)
                alpha_mods.append(CModSite(one_mod.name, m.site, one_mod.mass))
            alpha_pep = CPeptide(sq, p.pro_index, p.mass, alpha_mods, p.pro_index_list)
            alpha_sites = []
            for s in range(len(sq)):
                if sq[s] not in link_aa_dict: continue
                if s in mod_site_dict: continue
                if not link_site_valid(C_str, N_str, sq, s): continue
                alpha_sites.append(s)
            if len(alpha_sites) == 0: continue
            for s in range(len(alpha_sites)):
                link_peptide = CLinkPeptide(alpha_pep, beta_pep, alpha_sites[s], beta_site)
                #link_peptide.score, _ = match_PSM(spec, peak_index, all_peak_inten, link_peptide, conf)
                alpha_score, _, alpha_score_f = match_alpha_PSM(spec, peak_index, all_peak_inten, link_peptide, conf)
                #link_peptide.score = beta_score + alpha_score
                link_peptide.score = alpha_score
                link_peptide.alpha_score = alpha_score_f
                link_peptide.beta_score = beta_score_f
                link_peptide_list.append(link_peptide)
        end_time = time.time()
        print("[Debug]Time match for spectrum {0} is {1}".format(spec.title, (end_time-start_time)))
        cmpfun = operator.attrgetter('score')
        link_peptide_list.sort(key=cmpfun, reverse=True)
        link_peptide_list = update_peptide_list(link_peptide_list, protein_list)
        link_peptide_list = link_peptide_list[:conf.top_k_result]
        res.append(link_peptide_list)
    if len(res) != len(spectrum_list): print("[Error]#Results and #Spectra in {0} are not same, are {1} and {2}, respectively".format(pkl_file, len(res), len(spectrum_list)))
    assert(len(res) == len(spectrum_list))
    return res

def create_peak_index(peak_list, multi):
    end_mass = peak_list[-1].mz
    num = int(end_mass * multi) + 1
    index_list = [-1 for i in range(num)]
    all_peak_inten = 0.0
    for i in range(len(peak_list)):
        all_peak_inten += peak_list[i].inten
        p = peak_list[i]
        index = int(p.mz * multi)
        if index >= num: index = num - 1
        if index_list[index] == -1:
            index_list[index] = i

    end_val = len(peak_list)
    for i in range(num)[::-1]:
        if index_list[i] == -1: index_list[i] = end_val
        else: end_val = index_list[i]
    return index_list, all_peak_inten

def create_aa_list(sq, mods, aa2mass, beta_uua = 'a', beta_uua_mass = 0.0):
    # create the mass list of peptide sequence
    # e.g., peptide is "ACEDFK with C+57 modification"
    # return [M(A), M(C)+57, M(E), ..., M(K)] in which M(X) is the mass of X
    aa = [0.0 for i in range(len(sq))]
    for m in mods:
        aa[m.site-1] += m.mass
    for i, p in enumerate(sq):
        if p <'A' or p > 'Z': continue
        if p == beta_uua: 
            aa[i] += beta_uua_mass
        else:
            aa[i] += aa2mass[int(ord(p)-ord('A'))]
    return aa

def get_mod_str(mods):
    mod_str = ""
    for m in mods:
        mod_str += str(m.site) + "@" + str(m.mod_name) + ";"
    return mod_str

def get_str_mod(mod_str, mod2mass):
    mod_list = []
    for m in mod_str.split(';'):
        if m.find('@') < 0: continue
        site, name = m.split('@')
        site = int(site)
        one_mod = CModSite(name, site, mod2mass[name])
        mod_list.append(one_mod)
    return mod_list

def create_beta_peptide_theory_ions(pep2, site2, aa2mass, mass_P, mass_H2O, mass1):
    #mass1 is the mass of alpha peptide
    #pep2 is beta peptide
    sq2, mass2, mods2 = pep2.sq, pep2.mass, pep2.mods
    aa2 = create_aa_list(sq2, mods2, aa2mass)

    b2_mz, y2_mz = [], []
    b2_flag, y2_flag = [], []
    tmp_mass = 0.0
    for i in range(len(sq2) - 1):
        tmp_mass += aa2[i]
        if i < site2:
            b_charge1 = tmp_mass + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        else:
            b_charge1 = tmp_mass + mass1 + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        b2_mz.append(b_charge1)
        b2_mz.append(b_charge2)
        b2_flag.append("βb" + str(i+1) + "+")
        b2_flag.append("βb" + str(i+1) + "++")
    tmp_mass = 0.0
    for i in range(1, len(sq2))[::-1]:
        tmp_mass += aa2[i]
        if i > site2:
            y_charge1 = tmp_mass + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        else:
            y_charge1 = tmp_mass + mass1 + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        y2_mz.append(y_charge1)
        y2_mz.append(y_charge2)
        y_pos = len(sq2) - i
        y2_flag.append("βy" + str(y_pos) + "+")
        y2_flag.append("βy" + str(y_pos) + "++")
    return b2_mz + y2_mz, b2_flag + y2_flag

def create_alpha_peptide_theory_ions(link_peptide, aa2mass, mass_P, mass_H2O):
    pep1 = link_peptide.alpha_peptide
    pep2 = link_peptide.beta_peptide
    site1 = link_peptide.alpha_site
    sq1, mass1, mods1 = pep1.sq, pep1.mass, pep1.mods
    mass2 = pep2.mass
    aa1 = create_aa_list(sq1, mods1, aa2mass)
    
    b1_mz, y1_mz = [], []
    b1_flag, y1_flag = [], []
    tmp_mass = 0.0
    for i in range(len(sq1) - 1):
        tmp_mass += aa1[i]
        if i < site1:
            b_charge1 = tmp_mass + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        else:
            b_charge1 = tmp_mass + mass2 + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        b1_mz.append(b_charge1)
        b1_mz.append(b_charge2)
        b1_flag.append("αb" + str(i+1) + "+")
        b1_flag.append("αb" + str(i+1) + "++")
    tmp_mass = 0.0
    for i in range(1, len(sq1))[::-1]:
        tmp_mass += aa1[i]
        if i > site1:
            y_charge1 = tmp_mass + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        else:
            y_charge1 = tmp_mass + mass2 + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        y1_mz.append(y_charge1)
        y1_mz.append(y_charge2)
        y_pos = len(sq1) - i
        y1_flag.append("αy" + str(y_pos) + "+")
        y1_flag.append("αy" + str(y_pos) + "++")
    return b1_mz + y1_mz, b1_flag + y1_flag

def create_link_peptide_theory_ions(link_peptide, aa2mass, beta_uua, beta_uua_mass, mass_P, mass_H2O):
    # create the b and y ions (+1 and +2)
    pep1 = link_peptide.alpha_peptide
    pep2 = link_peptide.beta_peptide
    site1 = link_peptide.alpha_site
    site2 = link_peptide.beta_site
    sq1, mass1, mods1 = pep1.sq, pep1.mass, pep1.mods
    sq2, mass2, mods2 = pep2.sq, pep2.mass, pep2.mods
    aa1 = create_aa_list(sq1, mods1, aa2mass)
    aa2 = create_aa_list(sq2, mods2, aa2mass, beta_uua, beta_uua_mass)
    #print("=======START==============")
    #print("SQ1", sq1, get_mod_str(mods1))
    #print("SQ2", sq2, get_mod_str(mods2))
    #print("AA1", aa1)
    #print("AA2", aa2)
    #print("Sites", site1, site2)
    
    b1_mz, y1_mz, b2_mz, y2_mz = [], [], [], []
    b1_flag, y1_flag, b2_flag, y2_flag = [], [], [], []
    tmp_mass = 0.0
    for i in range(len(sq1) - 1):
        tmp_mass += aa1[i]
        if i < site1:
            b_charge1 = tmp_mass + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        else:
            b_charge1 = tmp_mass + mass2 + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        b1_mz.append(b_charge1)
        b1_mz.append(b_charge2)
        b1_flag.append("αb" + str(i+1) + "+")
        b1_flag.append("αb" + str(i+1) + "++")
    tmp_mass = 0.0
    for i in range(1, len(sq1))[::-1]:
        tmp_mass += aa1[i]
        if i > site1:
            y_charge1 = tmp_mass + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        else:
            y_charge1 = tmp_mass + mass2 + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        y1_mz.append(y_charge1)
        y1_mz.append(y_charge2)
        y_pos = len(sq1) - i
        y1_flag.append("αy" + str(y_pos) + "+")
        y1_flag.append("αy" + str(y_pos) + "++")

    tmp_mass = 0.0
    for i in range(len(sq2) - 1):
        tmp_mass += aa2[i]
        if i < site2:
            b_charge1 = tmp_mass + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        else:
            b_charge1 = tmp_mass + mass1 + mass_P
            b_charge2 = (b_charge1 + mass_P) / 2
        b2_mz.append(b_charge1)
        b2_mz.append(b_charge2)
        b2_flag.append("βb" + str(i+1) + "+")
        b2_flag.append("βb" + str(i+1) + "++")
    tmp_mass = 0.0
    for i in range(1, len(sq2))[::-1]:
        tmp_mass += aa2[i]
        if i > site2:
            y_charge1 = tmp_mass + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        else:
            y_charge1 = tmp_mass + mass1 + mass_P + mass_H2O
            y_charge2 = (y_charge1 + mass_P) / 2
        y2_mz.append(y_charge1)
        y2_mz.append(y_charge2)
        y_pos = len(sq2) - i
        y2_flag.append("βy" + str(y_pos) + "+")
        y2_flag.append("βy" + str(y_pos) + "++")
    
    #print("Ions b1", b1_mz)
    #print("Ions y1", y1_mz)
    #print("Ions b2", b2_mz)
    #print("Ions y2", y2_mz)
    return b1_mz + y1_mz + b2_mz + y2_mz, b1_flag + y1_flag + b2_flag + y2_flag

def _weight(len_norm, f, K1=0.08, constant_b=0.95):
    return (K1 + 1) * f / (K1 * (1 - constant_b + constant_b * len_norm) + f)

def match_between_peak_and_ions(peak_index, peaks, all_peak_inten, ions, ion_flags, is_ppm_fra, tol_fra, multi, norm=1.0, need_compute_conscore=True):
    # ions: b1+\b1++\b2+\b2++...y1+\y1++\y2+\y2++
    # need_compute_conscore: whether to compute the continue matched score
    # use score in pNovo 3
    match_num = 0
    match_inten = 0.0
    match_peak_index, match_peak_me = [], []
    match_N_flag = [0 for _ in range(int(len(ions)/4))]
    match_C_flag = [0 for _ in range(int(len(ions)/4))]
    match_score = 0.0
    max_me = 20.0
    if is_ppm_fra: max_me = tol_fra * 1e6
    else: max_me = tol_fra
    for mi,mz in enumerate(ions):
        if is_ppm_fra: tol_mass = mz * tol_fra
        else: tol_mass = tol_fra
        start_mass = mz - tol_mass
        end_mass = mz + tol_mass
        if end_mass <= 0: continue
        if start_mass > peaks[-1].mz: continue
        if start_mass < 0: start_mass = 0.0
        if end_mass > peaks[-1].mz: end_mass = peaks[-1].mz
        s_start_ind = peak_index[int(start_mass * multi)]
        while s_start_ind < len(peaks) and peaks[s_start_ind].mz < start_mass:
            s_start_ind += 1
        if int(end_mass * multi) >= len(peak_index): s_end_ind = peak_index[-1]
        else: s_end_ind = peak_index[int(end_mass * multi)]
        while s_end_ind < len(peaks) and peaks[s_end_ind].mz <= end_mass:
            s_end_ind += 1
        if s_start_ind >= s_end_ind: continue
        match_num += 1
        max_inten = 0.0
        max_inten_index = -1
        for i in range(s_start_ind, s_end_ind):
            if peaks[i].inten > max_inten: 
                max_inten = peaks[i].inten
                max_inten_index = i
        #compute Score
        
        if is_ppm_fra: 
            one_me = (peaks[max_inten_index].mz - mz) * 1e6 / mz
        else: 
            one_me = peaks[max_inten_index].mz - mz
            
        abs_one_me = abs(one_me)
        one_score = _weight(norm, math.log(max_inten)) * math.cos(abs_one_me / max_me * 3.1415926 / 2)
        match_score += one_score
        match_inten += max_inten
        match_peak_index.append((max_inten_index, ion_flags[mi]))
        match_peak_me.append(one_me)
        if mi < int(len(ions)/2):
            if mi % 2 == 0: match_N_flag[int(mi/2)] = one_score
        else:
            if mi % 2 == 0: match_C_flag[int((len(ions)-2-mi)/2)] = one_score
    match_ion_num_ratio = match_num * 1.0 / len(ions)
    match_inten_ratio = match_inten / all_peak_inten
    mean_frag_me, std_frag_me = max_me, 100.0
    if len(match_peak_me) > 0: mean_frag_me = abs(np.mean(match_peak_me))
    if len(match_peak_me) > 1: std_frag_me = np.std(match_peak_me)
    b_s, y_s, by_s = 0.0, 0.0, 0.0
    if compute_continue_score:
        b_s, y_s, by_s = compute_continue_score(match_N_flag, match_C_flag)
    single_score = CSingle_Score_Feature(match_score, mean_frag_me, std_frag_me, match_ion_num_ratio, match_inten_ratio, b_s, y_s, by_s)
    return single_score, match_peak_index

def compute_continue_score(match_N_flag, match_C_flag):
    b_s, y_s, by_s = 0, 0, 0
    TAG_LEN = 4
    for i in range(len(match_N_flag)-TAG_LEN+1):
        bCont = True
        for j in range(TAG_LEN):
            if match_N_flag[i + j] == 0.0: 
                bCont = False
                break
        if bCont:
            lfSum = 0.0
            for j in range(TAG_LEN): lfSum += match_N_flag[i + j]
            b_s += (lfSum / TAG_LEN)
        
        bCont = True
        for j in range(TAG_LEN):
            if match_C_flag[i + j] == 0.0:
                bCont = False
                break
        if bCont:
            lfSum = 0.0
            for j in range(TAG_LEN): lfSum += match_C_flag[i + j]
            y_s += (lfSum / TAG_LEN)

        bCont = True
        for j in range(TAG_LEN):
            if match_N_flag[i + j] == 0.0 and match_C_flag[i + j] == 0.0:
                bCont = False
                break
        if bCont:
            lfSum = 0.0
            for j in range(TAG_LEN): lfSum += (match_N_flag[i + j] + match_C_flag[i + j])
            by_s += (lfSum / TAG_LEN)
    return b_s, y_s, by_s


def match_PSM(spec, peak_index, all_peak_inten, link_peptide, conf):
    theory_ion_list, theory_ion_flag_list = create_link_peptide_theory_ions(link_peptide, conf.aa.aa2mass, conf.beta_uua, conf.beta_uua_mass, conf.ATOM_MASS_P, conf.MOLECULE_MASS_H2O)
    single_score, match_info = match_between_peak_and_ions(peak_index, spec.peaks,
                all_peak_inten, theory_ion_list, theory_ion_flag_list, conf.is_ppm_fra, conf.tol_fra, conf.multi, 1.0, False)
    return single_score.ion_ratio * single_score.inten_ratio * 100.0, match_info

def match_alpha_PSM(spec, peak_index, all_peak_inten, link_peptide, conf):
    alpha_ion_list, alpha_ion_flag_list = create_alpha_peptide_theory_ions(link_peptide, conf.aa.aa2mass, conf.ATOM_MASS_P, conf.MOLECULE_MASS_H2O)
    sq_len = float(len(link_peptide.alpha_peptide.sq))
    avrgd = link_peptide.alpha_peptide.mass / conf.AVRG_AA_MASS 
    norm = (sq_len * sq_len) / (avrgd * avrgd)
    single_score, match_info = match_between_peak_and_ions(peak_index, spec.peaks, all_peak_inten, alpha_ion_list, alpha_ion_flag_list, conf.is_ppm_fra, conf.tol_fra, conf.multi, norm)
    return single_score.score + single_score.continue_b_score + single_score.continue_y_score + single_score.continue_score, match_info, single_score
    #return single_score.ion_ratio, single_score.inten_ratio, match_info

def match_beta_PSM(spec, peak_index, all_peak_inten, beta_peptide, site2, conf, mass1 = None):
    if mass1 is None: #mass1 is the mass of alpha peptide
        mass1 = spec.mass - conf.ATOM_MASS_P - beta_peptide.mass
    beta_ion_list, beta_ion_flag_list = create_beta_peptide_theory_ions(beta_peptide, site2, conf.aa.aa2mass, conf.ATOM_MASS_P, conf.MOLECULE_MASS_H2O, mass1)
    sq_len = float(len(beta_peptide.sq))
    avrgd = beta_peptide.mass / conf.AVRG_AA_MASS
    norm = (sq_len * sq_len) / (avrgd * avrgd)
    #print("Beta_sq", beta_peptide.sq)
    single_score, match_info = match_between_peak_and_ions(peak_index, spec.peaks, all_peak_inten, beta_ion_list, beta_ion_flag_list, conf.is_ppm_fra, conf.tol_fra, conf.multi, norm)
    return single_score.score + single_score.continue_b_score + single_score.continue_y_score + single_score.continue_score, match_info, single_score
    #return single_score.ion_ratio, single_score.inten_ratio, match_info

def draw_one_psm(p):
    spec, title, link_peptide, conf = p
    peak_index, all_peak_inten = create_peak_index(spec.peaks, 1000)
    match_score, match_peak_info = match_PSM(spec, peak_index, all_peak_inten, link_peptide, conf)  
    #print("Score", match_score)
