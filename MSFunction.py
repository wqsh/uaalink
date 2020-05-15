#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: MSFunction.py
Author: webdev@tencent.com
Date: 2019/05/23 11:05:47
Brief: 
"""

import platform
import sys
#reload(sys)
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
try:
    import cPickle as pickle
except:
    import pickle
import operator
import os
import math
import copy

from MSData import CElement, CAA, CModification, CProtein, CPeptide, CModSite, CPeptidePKL, CModSitePKL, CPeak, CSpectrum, CLinkPeptide

class CConfigFunction:
    def config2file(self, path, config):
        with open(path, 'w') as f:
            f.write("#[data]\n")
            f.write("TYPE_MS2=" + config.ms2_type + "\n")
            f.write("PATH_MS2=" + config.mgf_file_list_str + "\n")
            f.write("PATH_FASTA=" + config.fasta_file + "\n")
            f.write("PATH_FASTA_EXPORT=" + config.index_out_folder + "\n")
            f.write("PATH_RESULT_EXPORT=" + config.result_folder + "\n")
            f.write("\n")
            f.write("#[biology]\n")
            f.write("NAME_ENZYME=" + config.enzyme_str + "\n")
            f.write("TYPE_DIGEST=" + config.enzyme_type + "\n")
            f.write("NUMBER_MAX_MISS_CLV=" + config.enzyme_miss_cle_num + "\n")
            f.write("\n")
            f.write("NAME_MOD_FIX=" + config.fix_mod_str + "\n")
            f.write("NAME_MOD_VAR=" + config.var_mod_str + "\n")
            f.write("NUMBER_MAX_MOD=" + config.var_mod_site_num + "\n")
            f.write("\n")
            f.write("UAA_SEQ=" + config.protein_beta + "\n")
            f.write("UAA_AA=" + config.beta_uua + "\n")
            f.write("UAA_LEN_LOW=" + config.beta_pep_len_min + "\n")
            f.write("UAA_LEN_UP=" + config.beta_pep_len_max + "\n")
            f.write("UAA_NAME_MOD_FIX=" + config.beta_fix_mod_str + "\n")
            f.write("UAA_NAME_MOD_VAR=" + config.beta_var_mod_str + "\n")
            f.write("UAA_COM=" + config.beta_uua_str + "\n")
            f.write("UAA_NAME_ENZYME=" + config.enzyme_str_beta + "\n")
            f.write("UAA_TYPE_DIGEST=" + config.enzyme_type_beta + "\n")
            f.write("UAA_NUMBER_MAX_MISS_CLV=" + config.enzyme_miss_cle_num_beta + "\n")
            f.write("UAA_LINKED_AA=" + config.link_aa + "\n")
            f.write("\n")
            f.write("#[mass spectrometry]\n")
            f.write("PPM_TOL_PRECURSOR=" + config.tol_pre_str + "\n")
            f.write("PPM_TOL_FRAGMENT=" + config.tol_fra_str + "\n")
            f.write("TYPE_ACTIVATION=" + config.activation_type + "\n")
            f.write("\n")
            f.write("#[performance]\n")
            f.write("NUMBER_THREAD=" + config.thread_num + "\n")
            f.write("TYPE_THREAD=" + config.thread_type + "\n")
            f.write("NUMBER_SELECT_PEAK=" + config.select_peak_num + "\n")
            f.write("NUMBER_SPECTRUM=" + config.spectrum_num + "\n")
            f.write("LEN_MAX_PROTEIN=" + config.protein_len_max + "\n")
            f.write("MASS_PEP_LOW=" + config.pep_mass_min + "\n")
            f.write("MASS_PEP_UP=" + config.pep_mass_max + "\n")
            f.write("LEN_PEP_LOW=" + config.pep_len_min + "\n")
            f.write("LEN_PEP_UP=" + config.pep_len_max + "\n")
            f.write("INDEX_SPLIT_MASS=" + config.index_split_mass + "\n")
            f.write("NUMBER_TOP_RESULT=" + config.top_k_result + "\n")
            f.write("\n")
            f.write("MULTI_MASS=" + config.multi + "\n")
            f.write("TYPE_TASK=" + config.task_type + "\n")
            f.write("TYPE_FILTER_BETA=" + config.is_filter_beta + "\n")
            f.write("NUMBER_PEAK_BETA=" + config.beta_filter_v + "\n")
            f.write("\n")
            f.write("OPEN_SEARCH_SINGLE=" + config.need_open_search + "\n")
            f.write("MASS_WINDOW_BETA=" + config.window_size + "\n")
            f.write("PATH_PFIND_RESULT=" + config.pfind_result_path + "\n")
            f.write("\n")
            f.write("#[filter]\n")
            f.write("FDR_PSM=" + config.fdr_value + "\n")
            f.write("\n")
            f.write("#[ini]\n")
            f.write("PATH_INI_ELEMENT=" + config.element_file + "\n")
            f.write("PATH_INI_AA=" + config.aa_file + "\n")
            f.write("PATH_INI_MOD=" + config.modification_file + "\n")

    def file2config(self, file_path, config):
        config.using_decoy = False
        config.cand_num = 10000
        config.enzyme_aa_beta = []
        sys_flag = platform.system()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.find('=') < 0: continue
                if line.find('#') >= 0:
                    line = line[:line.find('#')].strip()
                tmp = line.split('=')
                if len(tmp) < 2: continue
                name = tmp[0]
                value = tmp[1]
                print(name, ":", value)
                try:
                    if name == "PATH_INI_ELEMENT":
                        config.element_file = value
                        config.element = CElement(config.element_file)
                    elif name == "PATH_INI_AA":
                        try:
                            config.aa_file = value
                            config.aa = CAA(config.aa_file, config.element.element2mass, config.AA_NUM)
                        except Exception as err:
                            print("[Error] The element_file must be in front of AA_file.", err)
                            sys.exit(-1)
                    elif name == "PATH_INI_MOD":
                        config.modification_file = value
                        config.modification = CModification(config.modification_file)
                    elif name == "PATH_FASTA":
                        if sys_flag == "Linux": value = value.replace('\\', '/')
                        else: value = value.replace('/', '\\')
                        config.fasta_file = value
                    elif name == "PATH_FASTA_EXPORT":
                        if sys_flag == "Linux": value = value.replace('\\', '/')
                        else: value = value.replace('/', '\\')
                        config.index_out_folder = value
                        if not os.path.exists(config.index_out_folder):
                            os.makedirs(config.index_out_folder)
                    elif name == "PATH_RESULT_EXPORT":
                        if sys_flag == "Linux": value = value.replace('\\', '/')
                        else: value = value.replace('/', '\\')
                        config.result_folder = value
                        if not os.path.exists(config.result_folder):
                            os.makedirs(config.result_folder)
                    elif name == "INDEX_SPLIT_MASS":
                        config.index_split_mass = float(value)
                    elif name == "NUMBER_TOP_RESULT":
                        config.top_k_result = int(value)
                    elif name == "NUMBER_CAND":
                        config.cand_num = int(value)
                    elif name == "NAME_ENZYME":
                        value_list = value.split(';')
                        config.enzyme_name = []
                        config.enzyme_aa = []
                        config.enzyme_flag = []
                        for v in value_list:
                            if len(v.split(' ')) < 3: continue
                            tn, ta, tf = v.split(' ')
                            config.enzyme_name.append(tn)
                            config.enzyme_aa.append(ta)
                            config.enzyme_flag.append(tf)
                    elif name == "TYPE_DIGEST":
                        config.enzyme_type = int(value)
                    elif name == "NUMBER_MAX_MISS_CLV":
                        config.enzyme_miss_cle_num = int(value)
                    elif name == "MASS_PEP_LOW":
                        config.pep_mass_min = float(value)
                    elif name == "MASS_PEP_UP":
                        config.pep_mass_max = float(value)
                    elif name == "LEN_PEP_LOW":
                        config.pep_len_min = int(value)
                    elif name == "LEN_PEP_UP":
                        config.pep_len_max = int(value)
                    elif name == "NAME_MOD_FIX":
                        config.fix_mod_str = value
                    elif name == "NAME_MOD_VAR":
                        config.var_mod_str = value
                    elif name == "NUMBER_MAX_MOD":
                        config.var_mod_site_num = int(value)
                    elif name == "FDR_PSM":
                        config.fdr_value = float(value)
                        config.using_decoy = True
                    elif name == "PPM_TOL_PRECURSOR":
                        config.is_ppm_pre = True
                        if value.endswith("ppm") or value.endswith("PPM"):
                            config.tol_pre = float(value[:-3]) * 1e-6
                        elif value.endswith("da") or value.endswith("DA") or value.endswith("Da"):
                            config.is_ppm_pre = False
                            config.tol_pre = float(value[:-2])
                        else:
                            config.tol_pre = 20e-6 # default is 20ppm
                    elif name == "PPM_TOL_FRAGMENT":
                        config.is_ppm_fra = True
                        if value.endswith("ppm") or value.endswith("PPM"):
                            config.tol_fra = float(value[:-3]) * 1e-6
                        elif value.endswith("da") or value.endswith("DA") or value.endswith("Da"):
                            config.is_ppm_fra = False
                            config.tol_fra = float(value[:-2])
                        else:
                            config.tol_fra = 20e-6 # default is 20ppm
                    elif name == "PATH_PFIND_RESULT":
                        config.pfind_result_path = value.strip()
                    elif name == "NUMBER_THREAD":
                        config.thread_num = int(value)
                    elif name == "TYPE_THREAD":
                        config.thread_type = int(value)
                    elif name == "PATH_MS2":
                        if sys_flag == "Linux": value = value.replace('\\', '/')
                        else: value = value.replace('/', '\\')
                        if os.path.isdir(value):
                            config.mgf_file_list = []
                            config.mgf_folder = value
                        else:
                            config.mgf_file_list = value.split(';')
                    elif name == "TYPE_MS2":
                        config.ms2_type = value
                    elif name == "NUMBER_SELECT_PEAK":
                        config.select_peak_num = int(value)
                    elif name == "NUMBER_SPECTRUM":
                        config.spectrum_num = int(value)
                    elif name == "TYPE_TASK":
                        config.task_type = int(value)
                    elif name == "MULTI_MASS":
                        config.multi = int(value)
                    elif name == "TYPE_FILTER_BETA":
                        config.is_filter_beta = (int(value) == 1)
                    elif name == "NUMBER_PEAK_BETA":
                        config.beta_filter_v = int(value)
                    elif name == "OPEN_SEARCH_SINGLE":
                        config.need_open_search = int(value)
                    elif name == "MASS_WINDOW_BETA":
                        config.window_size = float(value)
                    elif name == "THRESHOLD_SCORE":
                        config.score_t = float(value)
                    elif name == "UAA_SEQ":
                        config.protein_beta = value
                    elif name == "UAA_AA":
                        config.beta_uua = value[0]
                    elif name == "UAA_COM":
                        config.beta_uua_str = value
                    elif name == "UAA_LEN_LOW":
                        config.beta_pep_len_min = int(value)
                    elif name == "UAA_LEN_UP":
                        config.beta_pep_len_max = int(value)
                    elif name == "UAA_NAME_MOD_FIX":
                        config.beta_fix_mod_str = value
                    elif name == "UAA_NAME_MOD_VAR":
                        config.beta_var_mod_str = value
                    elif name == "UAA_NAME_ENZYME":
                        value_list = value.split(';')
                        config.enzyme_name_beta = []
                        config.enzyme_aa_beta = []
                        config.enzyme_flag_beta = []
                        for v in value_list:
                            if len(v.split(' ')) < 3: continue
                            tn, ta, tf = v.split(' ')
                            config.enzyme_name_beta.append(tn)
                            config.enzyme_aa_beta.append(ta)
                            config.enzyme_flag_beta.append(tf)
                    elif name == "UAA_TYPE_DIGEST":
                        config.enzyme_type_beta = int(value)
                    elif name == "UAA_NUMBER_MAX_MISS_CLV":
                        config.enzyme_miss_cle_num_beta = int(value)
                    elif name == "UAA_LINKED_AA":
                        config.link_aa = value
                        config.link_aa_dict = {}
                        for one_link_aa in value:
                            config.link_aa_dict[one_link_aa] = 1
                except Exception as err:
                    print("[Error] occured in configure.txt", err)
                    sys.exit(-1)
        ret = self._init_mod(config)
        self._init_prime(config)
        self._init_beta(config)
        self._init_ms2(config)
        if ret != 0:
            print("[Error] occured in _init_mod")
            sys.exit(-1)

    def _init_ms2(self, config):
        if len(config.mgf_file_list) > 0 or config.mgf_folder is None:
            return
        for p in os.listdir(config.mgf_folder):
            if not p.endswith("." + config.ms2_type): continue
            config.mgf_file_list.append(os.path.join(config.mgf_folder, p))

    def _init_beta(self, config):
        aa = config.aa
        elem2mass = config.element.element2mass
        try:
            config.beta_uua_mass = float(config.beta_uua_str)
        except:
            config.beta_uua_mass = aa.compute_mass_from_element(config.beta_uua_str, elem2mass)

    def _init_fix_mod(self, config, fix_mod):
        fix_N_aa, fix_C_aa, fix_aa = {}, {}, {} #PEP_N, PEP_C, NORMAL
        fix_PN_aa, fix_PC_aa = {}, {} #PRO_N, PRO_C
        m2mod = config.modification.modname2mod
        for mod_name in fix_mod:
            if mod_name not in m2mod: 
                print("[Error] {} not in modification.ini".format(mod_name))
                return -1
            mod = m2mod[mod_name]
            if mod.type == "NORMAL":
                for a in mod.aa:
                    if a in fix_aa: return -1
                    fix_aa[a] = mod
            elif mod.type == "PEP_N":
                for a in mod.aa:
                    if a in fix_N_aa: return -1
                    fix_N_aa[a] = mod
            elif mod.type == "PEP_C":
                for a in mod.aa:
                    if a in fix_C_aa: return -1
                    fix_C_aa[a] = mod
            elif mod.type == "PRO_N":
                for a in mod.aa:
                    if a in fix_PN_aa: return -1
                    fix_PN_aa[a] = mod
            elif mod.type == "PRO_C":
                for a in mod.aa:
                    if a in fix_PC_aa: return -1
                    fix_PC_aa[a] = mod
        return fix_N_aa, fix_C_aa, fix_aa, fix_PN_aa, fix_PC_aa

    def _init_var_mod(self, config, var_mod):
        var_N_aa, var_C_aa, var_aa = {}, {}, {} #PEP_N, PEP_C, NORMAL
        var_PN_aa, var_PC_aa = {}, {} #PRO_N, PRO_C
        m2mod = config.modification.modname2mod
        for mod_name in var_mod:
            if mod_name not in m2mod:
                print("[Error] {} not in modification.ini".format(mod_name))
                return -1
            mod = m2mod[mod_name]
            if mod.type == "NORMAL":
                for a in mod.aa:
                    if a not in var_aa: var_aa[a] = []
                    var_aa[a].append(mod)
            elif mod.type == "PEP_N":
                for a in mod.aa:
                    if a not in var_N_aa: var_N_aa[a] = []
                    var_N_aa[a].append(mod)
            elif mod.type == "PEP_C":
                for a in mod.aa:
                    if a not in var_C_aa: var_C_aa[a] = []
                    var_C_aa[a].append(mod)
            elif mod.type == "PRO_N":
                for a in mod.aa:
                    if a not in var_PN_aa: var_PN_aa[a] = []
                    var_PN_aa[a].append(mod)
            elif mod.type == "PRO_C":
                for a in mod.aa:
                    if a not in var_PC_aa: var_PC_aa[a] = []
                    var_PC_aa[a].append(mod)
        return var_N_aa, var_C_aa, var_aa, var_PN_aa, var_PC_aa

    def _init_mod(self, config):
        # set fixed modifications and variable modifications
        if config.fix_mod_str.strip() == "":
            fix_mod = []
        else:
            fix_mod = config.fix_mod_str.split(';')
        if config.var_mod_str.strip() == "":
            var_mod = []
        else:
            var_mod = config.var_mod_str.split(';')

        if config.beta_fix_mod_str.strip() == "":
            beta_fix_mod = []
        else:
            beta_fix_mod = config.beta_fix_mod_str.split(';')
        if config.beta_var_mod_str.strip() == "":
            beta_var_mod = []
        else:
            beta_var_mod = config.beta_var_mod_str.split(';')

        # For self.fix_**: the key is the amino acid, value is the Mod object
        # For self.var_**: the key is the amino acid, value is the list of Mod objects
        config.fix_N_aa, config.fix_C_aa, config.fix_aa, config.fix_PN_aa, config.fix_PC_aa = self._init_fix_mod(config, fix_mod)
        config.var_N_aa, config.var_C_aa, config.var_aa, config.var_PN_aa, config.var_PC_aa = self._init_var_mod(config, var_mod)
        
        config.beta_fix_N_aa, config.beta_fix_C_aa, config.beta_fix_aa, config.beta_fix_PN_aa, config.beta_fix_PC_aa = self._init_fix_mod(config, beta_fix_mod)
        config.beta_var_N_aa, config.beta_var_C_aa, config.beta_var_aa, config.beta_var_PN_aa, config.beta_var_PC_aa = self._init_var_mod(config, beta_var_mod)
        
        #print("Fix", config.fix_N_aa, config.fix_C_aa, config.fix_aa, config.fix_PN_aa, config.fix_PC_aa)
        #print("Var", config.var_N_aa, config.var_C_aa, config.var_aa, config.var_PN_aa, config.var_PC_aa)
        #print("Fix-2", config.beta_fix_N_aa, config.beta_fix_C_aa, config.beta_fix_aa, config.beta_fix_PN_aa, config.beta_fix_PC_aa)
        #print("Var-2", config.beta_var_N_aa, config.beta_var_C_aa, config.beta_var_aa, config.beta_var_PN_aa, config.beta_var_PC_aa)
        
        # initial the modification index
        m2mod = config.modification.modname2mod
        config.modname2index = {}
        config.mod_list = []
        config.modname2type = {} # 0 is for fixed modification and 1 is for variable modification

        for mod_name in fix_mod:
            if mod_name in config.modname2index: return -1
            config.modname2index[mod_name] = len(config.mod_list)
            config.mod_list.append(m2mod[mod_name])
            config.modname2type[mod_name] = 0
        for mod_name in var_mod:
            if mod_name in config.modname2index: return -1
            config.modname2index[mod_name] = len(config.mod_list)
            config.mod_list.append(m2mod[mod_name])
            config.modname2type[mod_name] = 1
        return 0

    def _init_prime(self, config):
        var_mod = config.var_mod_str.split(';')
        prime_num = config.pep_len_max + config.var_mod_site_num + 1
        aa_num = config.AA_NUM + len(var_mod) + 1
        config.G_matrix = [] # gdm matrix (aa_num * prime_num)
        for i in range(aa_num):
            tmp = [0.0 for j in range(prime_num)]
            config.G_matrix.append(tmp)
        config.prime_list = [0 for j in range(prime_num)] 
        config.prime_list[0] = 2
        for i in range(1, prime_num):
            start_v = config.prime_list[i-1] + 1
            while True:
                is_prime = True
                v = int(math.sqrt(start_v))
                if v == 1: v = 2
                for k in range(v, start_v):
                    if start_v % k == 0:
                        is_prime = False
                        break
                if is_prime: break
                start_v += 1
            config.prime_list[i] = start_v

        for i in range(aa_num):
            for j in range(prime_num):
                config.G_matrix[i][j] = (i+1)*1.0*math.log(config.prime_list[j])
                #print(i,j,config.G_matrix[i][j])

def process_one_pkl(input_one_pkl):
    pkl_file, conf = input_one_pkl
    start_mass, end_mass = 0.0, 0.0
    if pkl_file.find('/') >= 0:
        start_end_str = pkl_file[pkl_file.rfind('/')+1:]
    else:
        start_end_str = pkl_file[pkl_file.rfind('\\')+1:]
    start_end_str = start_end_str[:start_end_str.find('.')]
    index_str = start_end_str.find('_')
    start_mass = int(start_end_str[:index_str]) * 1.0
    end_mass = int(start_end_str[index_str+1:]) * 1.0

    all_peptides = []
    f = open(pkl_file, 'rb')
    while True:
        try:
            one_peptide = pickle.load(f)
            all_peptides.append(one_peptide)
        except Exception: # the end of the file
            break
    f.close()

    cmpfun = operator.attrgetter('mass','gdm')
    all_peptides.sort(key=cmpfun)
    old_num = len(all_peptides)
    all_peptides = update_peptide(all_peptides)
    print("[Info]#Peptides in {0} file: {1} to {2}".format(pkl_file, old_num, len(all_peptides)))
    index_list = create_peptide_or_spectrum_index(all_peptides, start_mass, end_mass, conf.multi)
    
    fw = open(pkl_file, 'wb')
    pickle.dump(all_peptides, fw)
    fw.close()
    index_pkl_file = pkl_file[:pkl_file.rfind('.')] + "_ind.pkl"
    fw = open(index_pkl_file, 'wb')
    pickle.dump(index_list, fw)
    fw.close()

def update_peptide(all_peptides):
    # all_peptides is the list of peptides sorted by masses and gdm values
    # return the list of peptides after removing the peptides with the same gdm values and update the protein indexes
    start_i = 0
    new_all_peptides = []
    while start_i < len(all_peptides):
        new_all_peptides.append(all_peptides[start_i])
        cur_p = all_peptides[start_i]
        cur_mass = cur_p.mass
        cur_gdm = cur_p.gdm
        start_j = start_i + 1
        cur_pro_index_list = []
        cur_pro_index_list.append(cur_p.pro_index)
        while start_j < len(all_peptides) and all_peptides[start_j].mass == cur_mass and all_peptides[start_j].gdm == cur_gdm:
            cur_pro_index_list.append(all_peptides[start_j].pro_index)
            start_j += 1
        new_all_peptides[-1].pro_index_list = cur_pro_index_list
        start_i = start_j
    return new_all_peptides
    
def create_peptide_or_spectrum_index(all_peptides, start_mass, end_mass, mul=1): #mul=1000
    # all_peptides is the list of peptides sorted by masses and gdm values
    # return the list of mass index to the all_peptides
    num = int((end_mass - start_mass) * mul) + 1 # the number of list
    index_list = [-1 for i in range(num)] # index_list[X] is the min index in all_peptides whose mass is >= X
    for i in range(len(all_peptides)):
        p = all_peptides[i]
        index = int((p.mass - start_mass) * mul)
        if index >= num: index = num - 1
        if index_list[index] == -1:
            index_list[index] = i

    end_val = len(all_peptides)
    for i in range(num)[::-1]:
        if index_list[i] == -1: index_list[i] = end_val
        else: end_val = index_list[i]
    return index_list

def multi_process_pkl(pkl_file_list, thread_num, thread_type=0):
    if thread_type == 0:
        pool = Pool(processes=thread_num) # multiply processes
    else:
        pool = ThreadPool(thread_num) # multiply threads
    pool.map(process_one_pkl, pkl_file_list)
    pool.close()
    pool.join()

class CEnzymeFunction:
    def __init__(self, conf):
        self.conf = conf
        self.enzyme_name=conf.enzyme_name
        self.aa = conf.enzyme_aa
        self.split_flag = conf.enzyme_flag
        self.split_type = conf.enzyme_type
        self.miss_cleav_num = conf.enzyme_miss_cle_num
        self.pep_mass_min = conf.pep_mass_min
        self.pep_mass_max = conf.pep_mass_max
        self.pep_len_min = conf.pep_len_min
        self.pep_len_max = conf.pep_len_max
        for flag in self.split_flag:
            if flag != "N" and flag != "C":
                print("[Error] Enzyme split flag can only be N or C.")
                sys.exit(-1)
        if self.split_type < 0 or self.split_type > 2:
            print("[Error] Enzyme type can only be 0, 1, 2.")
            sys.exit(-1)
        if self.miss_cleav_num < 0:
            print("[Error] Enzyme miss cleavage number must be >= 0.")
            sys.exit(-1)
        if self.pep_mass_min < 0:
            print("[Error] Peptide mass minimum must be >= 0.")
            sys.exit(-1)
        if self.pep_len_min < 0:
            print("[Error] Peptide length minimum must be >= 0.")
            sys.exit(-1)

    def _find_split_site_ori(self, protein_seq):
        # find all split sites in protein sequence
        sites = []
        sites.append(-1)
        for i in range(len(protein_seq)):
            p = protein_seq[i]
            if self.aa.find(p) < 0: continue
            if self.split_flag == 'C': 
                if i != len(protein_seq)-1:
                    sites.append(i)
            elif self.split_flag == 'N':
                if i != 0:
                    sites.append(i-1)
        sites.append(len(protein_seq)-1)
        return sites

    def _find_split_site(self, protein_seq):
        sites = []
        site_dict = {}
        aa_split_type = {}
        for i, aa in enumerate(self.aa):
            for one_aa in aa:
                flag = self.split_flag[i]
                if flag == "C":
                    if one_aa not in aa_split_type: aa_split_type[one_aa] = 0
                    aa_split_type[one_aa] += 1
                elif flag == "N":
                    if one_aa not in aa_split_type: aa_split_type[one_aa] = 0
                    aa_split_type[one_aa] -= 1
        sites.append(-1)
        site_dict[-1] = 1
        for i in range(len(protein_seq)):
            p = protein_seq[i]
            if p not in aa_split_type: continue
            if aa_split_type[p] >= 0:
                if i != len(protein_seq) - 1 and i not in site_dict:
                    sites.append(i)
                    site_dict[i] = 1
            if aa_split_type[p] <= 0:
                if i != 0 and i-1 not in site_dict:
                    sites.append(i-1)
                    site_dict[i-1] = 1
        if len(protein_seq)-1 not in site_dict:
            sites.append(len(protein_seq)-1)
            site_dict[len(protein_seq)-1] = 1
        return sites

    def _specific_split(self, split_sites, protein_seq):
        # specific cut 
        # For example: protein sequence is "ACDEFKGGGRPSQED"
        # return: "ACDEFK", "GGGR", "PSQED"; [0,6],[6,10],[10,15]
        peptide_seq_list = []
        peptide_pos_list = [] #
        for i in range(len(split_sites)-1):
            start_site = split_sites[i] + 1
            for j in range(1, self.miss_cleav_num+2):
                if i+j >= len(split_sites): break
                end_site = split_sites[i+j]
                if not self._check_length(start_site, end_site+1): continue
                peptide_seq = protein_seq[start_site:end_site+1]
                peptide_seq_list.append(peptide_seq)
                peptide_pos_list.append([start_site, end_site+1])
        # split the protein when the first amino acid is "M"
        cur_len = len(peptide_pos_list)
        for i in range(cur_len):
            if peptide_pos_list[i][0] == 0 and peptide_seq_list[i].startswith("M"):
                end_site = peptide_pos_list[i][1]
                if not self._check_length(1, end_site): continue
                peptide_seq_list.append(peptide_seq_list[i][1:])
                peptide_pos_list.append([1, end_site])
        return peptide_seq_list, peptide_pos_list

    def _semi_specific_split(self, split_sites, protein_seq):
        peptide_seq_list = []
        peptide_pos_list = []
        return peptide_seq_list, peptide_pos_list

    def _non_specific_split(self, split_sites, protein_seq):
        peptide_seq_list = []
        peptide_pos_list = []
        return peptide_seq_list, peptide_pos_list

    def _check_mass(self, mass):
        # check the mass of peptide sequence is in [min_mass, max_mass]
        conf = self.conf
        return mass >= conf.pep_mass_min and mass <= conf.pep_mass_max

    def _check_length(self, start_site, end_site):
        # check the length of peptide sequence is in [min_len, max_len]
        len_peptide_seq = end_site - start_site
        return len_peptide_seq >= self.pep_len_min and len_peptide_seq <= self.pep_len_max

    def split_protein(self, protein_seq, print_flag=False):
        # split the protein sequence to list of CPeptide
        split_sites = self._find_split_site(protein_seq)
        if self.split_type == 0:
            peptide_seq_list, peptide_pos_list = self._specific_split(split_sites, protein_seq)
        elif self.split_type == 1:
            peptide_seq_list, peptide_pos_list = self._semi_specific_split(split_sites, protein_seq)
        else:
            peptide_seq_list, peptide_pos_list = self._non_specific_split(split_sites, protein_seq)

        if print_flag:
            print("[Info]Protein sq", protein_seq)
            for p in peptide_seq_list:
                print("[Info]Peptide sq", p)
        return peptide_seq_list, peptide_pos_list

class CProteinFunction:
    def __init__(self, enzymeFunc):
        self.enzymeFunc = enzymeFunc

    def load_protein(self, path):
        conf = self.enzymeFunc.conf
        using_decoy = conf.using_decoy
        protein_list = []
        ac, de, sq = "", "", ""
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "": continue
                if line[0] == '>':
                    if ac != "":
                        pro = CProtein(ac, de, sq)
                        protein_list.append(pro)
                        if using_decoy:
                            protein_list.append(self._generate_decoy_protein(ac, sq))
                    line = line[1:]
                    ind = line.find(' ')
                    if ind >= 0:
                        ac = line[:ind]
                    else:
                        ac = line
                    de = line[ind+1:]
                    sq = ""
                else:
                    sq += line
        if ac != "":
            pro = CProtein(ac, de, sq)
            protein_list.append(pro)
            if using_decoy:
                protein_list.append(self._generate_decoy_protein(ac, sq))
        return protein_list

    def _generate_decoy_protein(self, ac, sq):
        new_ac = "REV_" + ac
        new_sq = sq[::-1]
        decoy_protein = CProtein(new_ac, "", new_sq)
        return decoy_protein

    def _generate_modified_pep(self, protein_len, peptide_seq_list, peptide_pos_list, protein_index, fw_list, pickleFunc):
        conf = self.enzymeFunc.conf
        for i,pep_seq in enumerate(peptide_seq_list):
            pep_pos = peptide_pos_list[i]
            pep = self._generate_fix_one_pep(protein_len, pep_seq, pep_pos, protein_index)
            var_pep_list = self._generate_var_one_pep(protein_len, pep, pep_pos, protein_index)
            #print("======PEP=========", pep_seq)
            for p in var_pep_list:
                #p.print_peptide() # print the peptide information
                mass = p.mass
                fw_index = int((mass - conf.pep_mass_min) / conf.index_split_mass)
                new_mods = []
                for m in p.mods:
                    mod_name = m.mod_name
                    mod_ind = conf.modname2index[mod_name]
                    new_mods.append(CModSitePKL(mod_ind, m.site))
                gdm_value = self._get_peptide_gdm(p.sq, p.mods)
                new_p = CPeptidePKL(p.pro_index, pep_pos[0], pep_pos[1], p.mass, new_mods, gdm_value)
                pickleFunc.save_pep_to_pkl(new_p, fw_list[fw_index])

    def _dfs_var_site(self, var_sites, max_var_num, v_index, cur_sites, all_candidate):
        # var_sites is the list of variable sites
        # max_var_num is the maximum number of variable modifications in one peptide (=3)
        # all_candidate is the result of all peptides with and without variable modifications
        if v_index >= len(var_sites):
            all_candidate.append(cur_sites)
            return
        self._dfs_var_site(var_sites, max_var_num, v_index+1, cur_sites, all_candidate)
        if len(cur_sites) < max_var_num:
            new_sites = copy.deepcopy(cur_sites)
            new_sites.append(var_sites[v_index])
            self._dfs_var_site(var_sites, max_var_num, v_index+1, new_sites, all_candidate)

    def _dfs_var_pep(self, cand, var_mod, i, cur_mod, var_mod_list):
        if i == len(cand):
            var_mod_list.append(cur_mod)
            return
        mod_list = var_mod[cand[i]]
        for m in mod_list:
            mod_site = CModSite(m.name, cand[i], m.mass)
            new_mod = copy.deepcopy(cur_mod)
            new_mod.append(mod_site)
            self._dfs_var_pep(cand, var_mod, i+1, new_mod, var_mod_list)

    def _generate_var_peps(self, all_candidate, var_mod, fixed_peptide):
        var_peptide_list = []
        for cand_i, cand in enumerate(all_candidate):
            # all candidate of variable sites
            # For example: sq = "MVLTIYPDELVQIVSDK"
            # all_candidate is [[], [1]]
            if len(cand) == 0:
                non_var_pep = copy.deepcopy(fixed_peptide)
                if self.enzymeFunc._check_mass(non_var_pep.mass):
                    var_peptide_list.append(non_var_pep)
                continue
            var_mod_list = []
            
            self._dfs_var_pep(cand, var_mod, 0, [], var_mod_list)
            # consider Acetyl[AnyN-term], Oxidation[M] and Dioxidation[M]
            # we can get three peptides with modification:
            # MVLTIYPDELVQIVSDK with Ace, Oxi and Dioxi
            for one_mod_list in var_mod_list:
                new_var_pep = copy.deepcopy(fixed_peptide)
                new_var_pep.mods += one_mod_list
                self._update_peptide_mass(new_var_pep, one_mod_list)
                if self.enzymeFunc._check_mass(new_var_pep.mass):
                    var_peptide_list.append(new_var_pep)
        return var_peptide_list

    def _generate_var_one_pep(self, protein_len, fixed_peptide, peptide_pos, protein_index):
        # According to the original peptide sequence to generate multiply peptides with variable modifications
        conf = self.enzymeFunc.conf
        var_aa = conf.var_aa
        var_N_aa, var_C_aa = conf.var_N_aa, conf.var_C_aa
        var_PN_aa, var_PC_aa = conf.var_PN_aa, conf.var_PC_aa
        pep_seq = fixed_peptide.sq
        var_site = [0 for i in range(len(pep_seq)+2)]
        var_mod = [[] for i in range(len(var_site))]
        
        if len(var_N_aa) > 0:
            if pep_seq[0] in var_N_aa:
                mod_list = var_N_aa[pep_seq[0]]
                var_site[0] = 1
                var_mod[0] = copy.deepcopy(mod_list)
        if len(var_PN_aa) > 0:
            if peptide_pos[0] == 0 and pep_seq[0] in var_PN_aa:
                mod_list = var_PN_aa[pep_seq[0]]
                var_site[0] = 1
                var_mod[0] = copy.deepcopy(mod_list)
        if len(var_C_aa) > 0:
            if pep_seq[-1] in var_C_aa:
                mod_list = var_C_aa[pep_seq[-1]]
                var_site[-1] = 1
                var_mod[-1] = copy.deepcopy(mod_list)
        if len(var_PC_aa) > 0:
            if peptide_pos[1] == protein_len and pep_seq[-1] in var_PC_aa:
                mod_list = var_PC_aa[pep_seq[-1]]
                var_site[-1] = 1
                var_mod[-1] = copy.deepcopy(mod_list)
        if len(var_aa) > 0:
            for i in range(len(pep_seq)):
                a = pep_seq[i]
                if a not in var_aa: continue
                mod_list = var_aa[a]
                var_site[i + 1] = 1
                var_mod[i + 1] = copy.deepcopy(mod_list)
        # If the site is used for fixed modification, then this site cannot be used for variable modification.
        fixed_mods = fixed_peptide.mods
        for m in fixed_mods:
            var_site[m.site] = 0
            var_mod[m.site] = []
        # N terminal modification and the variable modification occured on the 1st amino acid are the same.
        var_site[1] = var_site[0]
        var_site[0] = 0
        var_mod[1] += var_mod[0]
        var_mod[0] = []
        # C terminal modification is similar.
        var_site[-2] = var_site[-1]
        var_site[-1] = 0
        var_mod[-2] += var_mod[-1]
        var_mod[-1] = []

        all_var_site_list = []
        for i in range(len(var_site)):
            if var_site[i] == 0: continue
            all_var_site_list.append(i)
        all_candidate = []
        self._dfs_var_site(all_var_site_list, conf.var_mod_site_num, 0, [], all_candidate)

        return self._generate_var_peps(all_candidate, var_mod, fixed_peptide)

    def _generate_fix_one_pep(self, protein_len, peptide_seq, peptide_pos, protein_index):
        # According to the original peptide sequence to generate a new peptide with fixed modifications
        conf = self.enzymeFunc.conf
        fix_aa = conf.fix_aa
        fix_N_aa, fix_C_aa = conf.fix_N_aa, conf.fix_C_aa
        fix_PN_aa, fix_PC_aa = conf.fix_PN_aa, conf.fix_PC_aa
        # generate N-terminal, C-terminal and normal fixed modification
        fix_mods = []
        if len(fix_N_aa) > 0:
            if peptide_seq[0] in fix_N_aa:
                mod = fix_N_aa[peptide_seq[0]]
                modSite = CModSite(mod.name, 0, mod.mass)
                fix_mods.append(modSite)
        if len(fix_PN_aa) > 0:
            if peptide_pos[0] == 0 and peptide_seq[0] in fix_PN_aa:
                mod = fix_PN_aa[peptide_seq[0]]
                modSite = CModSite(mod.name, 0, mod.mass)
                fix_mods.append(modSite)
        if len(fix_C_aa) > 0:
            if peptide_seq[-1] in fix_C_aa:
                mod = fix_C_aa[peptide_seq[-1]]
                modSite = CModSite(mod.name, len(peptide_seq) + 1, mod.mass)
                fix_mods.append(modSite)
        if len(fix_PC_aa) > 0:
            if peptide_pos[1] == protein_len and peptide_seq[-1] in fix_PC_aa:
                mod = fix_PC_aa[peptide_seq[-1]]
                modSite = CModSite(mod.name, len(peptide_seq) + 1, mod.mass)
                fix_mods.append(modSite)
        if len(fix_aa) > 0:
            for i in range(len(peptide_seq)):
                a = peptide_seq[i]
                if a not in fix_aa: continue
                mod = fix_aa[a]
                modSite = CModSite(mod.name, i + 1, mod.mass)
                fix_mods.append(modSite)
        mass = self._get_peptide_mass(peptide_seq, fix_mods)
        peptide_res = CPeptide(peptide_seq, protein_index, mass, fix_mods)
        return peptide_res

    def _get_peptide_mass(self, peptide_seq, mods):
        conf = self.enzymeFunc.conf
        aa = conf.aa.aa2mass # mass of amino acid (A-Z)
        mass = conf.MOLECULE_MASS_H2O
        for p in peptide_seq:
            if p < 'A' or p > 'Z': continue
            mass += aa[int(ord(p)-ord('A'))]
        for m in mods:
            mass += m.mass
        return mass

    def _update_peptide_mass(self, peptide, var_mods):
        for m in var_mods:
            peptide.mass += m.mass

    def _get_peptide_gdm(self, peptide_seq, mods):
        # get the gdm value of peptide sequence with variable modifications
        # The same peptide sequence with same modifications and sites have the same gdm values
        conf = self.enzymeFunc.conf
        m2type = conf.modname2type
        m2index = conf.modname2index
        aa_num = conf.AA_NUM
        var_mods_dict = {}
        for m in mods:
            if m2type[m.mod_name] == 0: continue #fix modification
            var_mods_dict[m.site] = m2index[m.mod_name]
        new_sq_list = []
        for i in range(len(peptide_seq)):
            if i in var_mods_dict:
                new_sq_list.append(aa_num + var_mods_dict[i])
            new_sq_list.append(int(ord(peptide_seq[i]) - ord('A')))

        gdm_value = 0.0
        for i in range(len(new_sq_list)):
            ind = i
            if ind >= len(conf.G_matrix[0]): ind = len(conf.G_matrix[0]) - 1
            ind2 = new_sq_list[i]
            if ind2 >= len(conf.G_matrix): ind2 = len(conf.G_matrix) - 1
            gdm_value += conf.G_matrix[ind2][ind]
        return gdm_value

    def split_protein(self, protein_list, print_flag=False):
        conf = self.enzymeFunc.conf
        index_folder = conf.index_out_folder
        for p in os.listdir(index_folder):
            if p.endswith(".pkl"): os.remove(os.path.join(index_folder, p))
        pkl_num = int((conf.pep_mass_max - conf.pep_mass_min) / conf.index_split_mass) + 1

        fw_list, pkl_file_list = [], []
        for i in range(pkl_num):
            start_mass = conf.pep_mass_min + i * conf.index_split_mass
            end_mass = conf.pep_mass_min + (i + 1) * conf.index_split_mass
            if end_mass > conf.pep_mass_max: end_mass = conf.pep_mass_max
            if start_mass == end_mass: continue
            pkl_file = os.path.join(index_folder, "{0}_{1}.pkl".format(int(start_mass), int(end_mass)))
            fw = open(pkl_file, "wb")
            fw_list.append(fw)
            pkl_file_list.append((pkl_file, conf))
        pickleFunc = CPickleFunction()

        for pro_index, protein in enumerate(protein_list):
            protein_sq = protein.sq
            peptide_seq_list, peptide_pos_list = self.enzymeFunc.split_protein(protein_sq)
            if print_flag:
                print("[Info]#Peptides in {}".format(protein.ac), len(peptide_seq_list))
                for i,p in enumerate(peptide_seq_list):
                    print("[Info]Peptide sq", p, peptide_pos_list[i])

            self._generate_modified_pep(len(protein_sq), peptide_seq_list, peptide_pos_list, pro_index, fw_list, pickleFunc)
            if pro_index % 1000 == 0:
                print("[Info]Process: {}%".format(pro_index*100/len(protein_list)))
        for p in fw_list:
            p.close()
       
        # Sort the peptides by masses and gdm values
        # Merge the same peptide (same modifications) and update the protein indexes
        # Create the peptide index file (*_*_ind.pkl file)
        multi_process_pkl(pkl_file_list, conf.thread_num, conf.thread_type)

class CPickleFunction:
    def __init__(self):
        pass

    def save_pep_to_pkl(self, peptide, f):
        pickle.dump(peptide, f)

    def save_pro_to_pkl(self, protein_list, path):
        f = open(path, 'wb')
        pickle.dump(protein_list, f)
        f.close()

    def save_mod_to_pkl(self, mod_list, path):
        f = open(path, 'wb')
        pickle.dump(mod_list, f)
        f.close()

    def save_psm_to_pkl(self, psm_list, path):
        f = open(path, 'wb')
        pickle.dump(psm_list, f)
        f.close()

    def load_pro_from_pkl(self, path):
        f = open(path, 'rb')
        protein_list = pickle.load(f)
        f.close()
        return protein_list

    def get_start_end_mass(self, path):
        # From pkl file path to get the start and end mass
        if path.find('/') >= 0:
            filename = path[path.rfind('/')+1:]
        else:
            filename = path[path.rfind('\\')+1:]
        filename = filename[:filename.find('.')]
        index = filename.find('_')
        start_mass = float(filename[:index])
        end_mass = float(filename[index+1:])
        return start_mass, end_mass

    def get_all_pkl_files(self, folder, pro_pkl_file, mod_pkl_file):
        pkl_file_list, ind_file_list = [], []
        for p in os.listdir(folder):
            if not p.endswith('.pkl'): continue
            if p == pro_pkl_file: continue
            if p == mod_pkl_file: continue
            if p.endswith("_ind.pkl"):
                #ind_file_list.append(os.path.join(folder, p))
                pass
            else:
                pkl_file_list.append(os.path.join(folder, p))
                ind_file_list.append(os.path.join(folder, p[:-4] + "_ind.pkl"))
        assert(len(pkl_file_list) == len(ind_file_list))
        return pkl_file_list, ind_file_list

class COutputFileFunction:
    def __init__(self):
        pass

    def write_label(self, psm_list, mgf_file_list, folder, conf):
        mgf_name_to_psm_list = {}
        for psm in psm_list:
            title = psm.spec.title
            title = title[:title.find('.')]
            if title not in mgf_name_to_psm_list:
                mgf_name_to_psm_list[title] = []
            mgf_name_to_psm_list[title].append(psm)

        for mgf_path in mgf_file_list:
            if mgf_path.find('/') >= 0:
                mgf_name = mgf_path[mgf_path.rfind('/')+1:]
            else:
                mgf_name = mgf_path[mgf_path.rfind('\\')+1:]
            mgf_name = mgf_name[:mgf_name.rfind('.')]
            if mgf_name.find("_") >= 0:
                mgf_name = mgf_name[:mgf_name.rfind('_')]
            self.write_one_label(mgf_name_to_psm_list[mgf_name], mgf_path, folder, conf)

    def write_one_label(self, psm_list, mgf_path, folder, conf, xlink_flag="UAA"):
        if mgf_path.find('/') >= 0:
            mgf_name = mgf_path[mgf_path.rfind('/')+1:]
        else:
            mgf_name = mgf_path[mgf_path.rfind('\\')+1:]
        mgf_name = mgf_name[:mgf_name.rfind('.')]
        if mgf_name.find("_") >= 0:
            mgf_name = mgf_name[:mgf_name.rfind('_')]
        path = os.path.join(folder, mgf_name+".plabel")
        with open(path, 'w') as fw:
            fw.write("[FilePath]\n")
            fw.write("File_Path="+mgf_path+'\n')
            fw.write("[Modification]\n")
            mod_list = conf.mod_list
            mod_to_index = {}
            for i, mod in enumerate(mod_list):
                modname = mod.name
                mod_to_index[modname] = i+1
                fw.write(str(i+1) + '=' + modname + '\n')
            fw.write("[xlink]\n")
            fw.write("xlink=" + xlink_flag + '\n')
            fw.write("[Total]\n")
            fw.write("total=" + str(len(psm_list)) + '\n')
            for i, psm in enumerate(psm_list):
                pep = psm.peptide
                sq1, mods1, site1 = pep.alpha_peptide.sq, pep.alpha_peptide.mods, pep.alpha_site
                sq2, mods2, site2 = pep.beta_peptide.sq, pep.beta_peptide.mods, pep.beta_site
                site1 += 1
                site2 += 1
                fw.write("[Spectrum" + str(i+1) + "]\n")
                fw.write("name=" + psm.spec.title + '\n')
                fw.write("pep1=3 " + str(site1) + ' ' + str(site2) + ' ' + sq1 + ' ' + str(pep.score) + ' ' + sq2 + ' 1 ')
                mod_str = self._get_mod_str2(mods1, 0, mod_to_index)
                mod_str += self._get_mod_str2(mods2, len(sq1) + 3, mod_to_index)
                fw.write(mod_str + '\n')

    def write_single_psms(self, psm_list, path, conf, is_write_candidate=False):
        fw = open(path, 'w')
        fw.write("Title\tMass\tSingle peptide sequence\tSingle peptide modification mass\tSingle peptide site\t" +
                 "Score\tDelta score\tRank\tTarget-Decoy\n")
        for p in psm_list:
            title = p.spec.title
            mass = p.spec.mass
            delta_score = p.peptide.score
            if len(p.peptide_list) > 1:
                delta_score -= p.peptide_list[1].score
            is_target = (p.peptide.pro_index == 0)
            fw.write(title + '\t' + str(mass) + '\t' + p.peptide.peptide_sq + '\t' + str(p.peptide.peptide_modmass) +
                     '\t' + str(p.peptide.peptide_site + 1) + '\t' + str(p.peptide.score) + '\t' + str(delta_score) + '\t1\t' +
                     str(is_target) + '\n')
        fw.close()

    def write_psms(self, psm_list, path, conf, is_write_candidate=False):
        fw = open(path, 'w')
        fw.write("Title\tMass\tAlpha peptide sequence\tAlpha peptide modification\tAlpha site\tAlpha mass\t" + 
                "Beta peptide sequence\tBeta peptide modification\tBeta site\tBeta mass\t" + 
                "Match tolerance (ppm/Da)\tProteins\tScore\tDelta score\tRank\t" +
                 "Alpha score\tAlpha mean fragment mass tol\t" +
                 "Alpha std fragment mass tol\tAlpha ion ratio\tAlpha inten ratio\t" +
                 "Alpha b-score\tAlpha y-score\tAlpha b/y-score\tAlpha sq len\t" +
                 "Beta score\tBeta mean fragment mass tol\t" +
                 "Beta std fragment mass tol\tBeta ion ratio\tBeta inten ratio\t" +
                 "Beta b-score\tBeta y-score\tBeta b/y-score\tBeta sq len\n")
        self.mass_P = conf.ATOM_MASS_P
        pickleFunc = CPickleFunction()
        self.protein_list = pickleFunc.load_pro_from_pkl(os.path.join(conf.index_out_folder, conf.PRO_PKL_FILE))
        for p in psm_list:
            title = p.spec.title
            mass = p.spec.mass
            delta_score = p.peptide.score
            if len(p.peptide_list) > 1:
                delta_score -= p.peptide_list[1].score
            self._write_one_psm(title, mass, p.peptide, conf.is_ppm_pre, fw, delta_score)
            if is_write_candidate:
                for i in range(1,len(p.peptide_list)):
                    self._write_one_psm(title, mass, p.peptide_list[i], conf.is_ppm_pre, fw, 0.0, i+1)
        fw.close()

    def _get_mod_str2(self, mods, start_site, mod_2_index):
        mod_str = ""
        for m in mods:
            site = m.site + start_site
            name = m.mod_name
            index = mod_2_index[name]
            if name.find('N-term') >= 0: site -= 1
            elif name.find('C-term') >= 0: site += 1
            mod_str += str(site) + ',' + str(index) + ' '
        return mod_str

    def _get_mod_str(self, mods):
        mod_str = ""
        for m in mods:
            mod_str += str(m.site) + "@" + str(m.mod_name) + ";"
        if len(mod_str) > 0 and mod_str[-1] == '\t':
            mod_str = mod_str[:-1]
        return mod_str

    def _write_one_psm(self, title, mass, pep, is_ppm_pre, fw, delta_score = 0.0, rank=1):
        sq1, mods1, site1, mass1, score_f1 = pep.alpha_peptide.sq, pep.alpha_peptide.mods, pep.alpha_site, pep.alpha_peptide.mass, pep.alpha_score
        sq2, mods2, site2, mass2, score_f2 = pep.beta_peptide.sq, pep.beta_peptide.mods, pep.beta_site, pep.beta_peptide.mass, pep.beta_score
        mods1 = self._get_mod_str(mods1)
        mods2 = self._get_mod_str(mods2)
        site1 += 1 # start from 1 not 0
        site2 += 1
        ac_list = [self.protein_list[i].ac for i in pep.alpha_peptide.pro_index_list]
        ac_str = ";".join(ac_list)
        if is_ppm_pre: match_tolerance = (mass - mass1 - mass2 - self.mass_P) * 1e6 / mass
        else: match_tolerance = (mass - mass1 - mass2 - self.mass_P)
        fw.write(title + '\t' + str(mass) + '\t' + sq1 + '\t' + mods1 + '\t' + str(site1) + '\t' + str(mass1) + '\t' +
                 sq2 + '\t' + mods2 + '\t' + str(site2) + '\t' + str(mass2) + '\t' +
                 str(match_tolerance) + '\t' + ac_str + '\t' + str(pep.score) + '\t' + str(delta_score) + '\t' + str(rank) + '\t' + score_f1.get_string() + '\t' + str(len(sq1)) + '\t' + score_f2.get_string() + '\t' + str(len(sq2)) + '\n')

class CSpectrumFunction:
    def __init__(self, conf):
        self.conf = conf

    def _load_pfind_spectrum(self, pfind_path, q_value_t=0.01):
        Q_INDEX = 4
        AC_INDEX = 12
        FLAG_INDEX = 15
        T_INDEX = 0
        title_set = set()
        if pfind_path == "" or not os.path.exists(pfind_path):
            return title_set
        with open(pfind_path) as f:
            lines = f.readlines()
        for i in range(1, len(lines)):
            tmp = lines[i].strip().split('\t')
            q_value = float(tmp[Q_INDEX])
            if q_value > q_value_t: break
            is_target = (tmp[FLAG_INDEX] == "target")
            if is_target: title_set.add(tmp[T_INDEX])
        return title_set

    def _load_single_spectrum(self, single_res_path):
        title_set = set()
        if single_res_path == "" or not os.path.exists(single_res_path):
            return title_set
        with open(single_res_path) as f:
            lines = f.readlines()
        for i in range(1, len(lines)):
            tmp = lines[i].strip().split('\t')
            title = tmp[0]
            title_set.add(title)
        return title_set

    def load_spectrum(self, path):
        conf = self.conf
        spectrum_list = []
        title_set = self._load_pfind_spectrum(conf.pfind_result_path)
        print("[Info]#PSMs in pFind is {}".format(len(title_set)))
        #title_set2 = self._load_single_spectrum(conf.single_res_path)
        #print("[Info]#PSMs in this workflow is {}".format(len(title_set2)))
        #title_set = title_set | title_set2
        title = ""
        mass = 0.0
        charge = 0
        peak_list = []
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "": continue
                if line[0] == 'T':
                    title = line[6:]
                elif line[0] == 'P':
                    mass = float(line[8:])
                elif line[0] == 'C':
                    if line[-1] == '+':
                        line = line[:-1]
                    charge = int(line[7])
                elif line[0] >= '1' and line[0] <= '9':
                    ind = line.find(' ')
                    p_mz = float(line[:ind])
                    p_inten = float(line[ind+1:])
                    peak_list.append(CPeak(p_mz, p_inten))
                elif line[0] == 'E':
                    mass = mass * charge - (charge - 1) * conf.ATOM_MASS_P
                    if len(peak_list) > conf.select_peak_num:
                        cmpfun = operator.attrgetter('inten')
                        peak_list.sort(key=cmpfun, reverse=True)
                        peak_list = peak_list[:conf.select_peak_num]
                    cmpfun = operator.attrgetter("mz")
                    peak_list.sort(key=cmpfun)
                    spec = CSpectrum(title, mass, charge, peak_list)
                    if len(spec.peaks) > 0 and title not in title_set:
                        spectrum_list.append(spec)
                    title = ""
                    mass = 0.0
                    charge = 0
                    peak_list = []
        return spectrum_list

class CBetaProteinFunction:
    def __init__(self, conf):
        self.protein_seq = conf.protein_beta
        self.enzyme_aa_list = conf.enzyme_aa_beta
        self.enzyme_flag_list = conf.enzyme_flag_beta
        self.enzyme_type = conf.enzyme_type_beta
        self.enzyme_miss_cle_num = conf.enzyme_miss_cle_num_beta
        self.conf = conf

    def _find_split_site(self):
        sq = self.protein_seq
        sites = []
        site_dict = {}
        aa_split_type = {}
        for i, aa in enumerate(self.enzyme_aa_list):
            for one_aa in aa:
                flag = self.enzyme_flag_list[i]
                if flag == "C":
                    if one_aa not in aa_split_type: aa_split_type[one_aa] = 0
                    aa_split_type[one_aa] += 1
                elif flag == "N":
                    if one_aa not in aa_split_type: aa_split_type[one_aa] = 0
                    aa_split_type[one_aa] -= 1
        print("AA", aa_split_type)
        sites.append(-1)
        site_dict[-1] = 1
        for i in range(len(sq)):
            p = sq[i]
            if p not in aa_split_type: continue
            if aa_split_type[p] >= 0:
                if i != len(sq) - 1 and i not in site_dict:
                    sites.append(i)
                    site_dict[i] = 1
            if aa_split_type[p] <= 0:
                if i != 0 and i-1 not in site_dict:
                    sites.append(i-1)
                    site_dict[i-1] = 1
        if len(sq)-1 not in site_dict:
            sites.append(len(sq)-1)
            site_dict[len(sq)-1] = 1
        return sites

    def _contain_link_site(self, pep_sq):
        return pep_sq.find(self.conf.beta_uua) >= 0

    def _specific_split(self, split_sites):
        print("[Info]Split sites", split_sites)
        sq = self.protein_seq
        peptide_seq_list, peptide_pos_list = [], []
        conf = self.conf
        for i in range(len(split_sites)-1):
            start_site = split_sites[i] + 1
            for j in range(1, self.enzyme_miss_cle_num+2):
                if i+j >= len(split_sites): break
                end_site = split_sites[i+j]
                peptide_seq = sq[start_site:end_site+1]
                if self._contain_link_site(peptide_seq) and len(peptide_seq) >= conf.beta_pep_len_min and len(peptide_seq) <= conf.beta_pep_len_max:
                    peptide_seq_list.append(peptide_seq)
                    peptide_pos_list.append([start_site, end_site+1])
        return peptide_seq_list, peptide_pos_list

    def _semi_specific_split(self, split_sites):
        sq = self.protein_seq
        peptide_seq_list = []
        peptide_pos_list = []
        return peptide_seq_list, peptide_pos_list

    def _non_specific_split(self, split_sites=None):
        sq = self.protein_seq
        target_seq_list, _ = self._non_specific_split_by_sq(sq)
        decoy_seq_list, _ = self._non_specific_split_by_sq(sq[::-1])
        return target_seq_list, decoy_seq_list 

    def _check_site(self, protein_sq, start_pos, end_pos, N_site_aa, C_site_aa):
        site_num = 0 # the number of sites which is consistent with the enzymes
        if protein_sq[start_pos] in N_site_aa: site_num += 1
        if start_pos <= 0 or protein_sq[start_pos - 1] in C_site_aa: site_num += 1
        if protein_sq[end_pos] in C_site_aa: site_num += 1
        if end_pos + 1 >= len(protein_sq) or protein_sq[end_pos + 1] in N_site_aa: site_num += 1
        return (site_num >= 1)

    def _non_specific_split_by_sq(self, sq, need_one_site=True):
        #need_one_site is True: one site is enzyme split site
        if need_one_site:
            N_site_aa, C_site_aa = set(), set()
            for i in range(len(self.conf.enzyme_aa_beta)):
                flag = self.conf.enzyme_flag_beta[i]
                aa = self.conf.enzyme_aa_beta[i]
                if flag == "N":
                    for c in aa:
                        N_site_aa.add(c)
                elif flag == "C":
                    for c in aa:
                        C_site_aa.add(c)

        peptide_seq_list = []
        peptide_pos_list = []
        index = sq.find(self.conf.beta_uua)
        if index < 0: return peptide_seq_list, peptide_pos_list
        pep_min = self.conf.beta_pep_len_min
        pep_max = self.conf.beta_pep_len_max
        start_index = max(0, index - pep_max)
        end_index = min(len(sq), index + pep_max)
        for i in range(start_index, index + 1):
            for j in range(index, end_index):
                if j - i + 1 > pep_max or j - i + 1 < pep_min: continue
                pep_sq = sq[i: j + 1]
                if need_one_site:
                    if self._check_site(sq, i, j, N_site_aa, C_site_aa):
                        peptide_seq_list.append(pep_sq)
                        peptide_pos_list.append([i, j+1])
                else:
                    peptide_seq_list.append(pep_sq)
                    peptide_pos_list.append([i, j+1])
        return peptide_seq_list, peptide_pos_list

    def _get_peptide_mass(self, peptide_seq, mods=[]):
        aa = self.conf.aa.aa2mass
        mass = self.conf.MOLECULE_MASS_H2O
        for p in peptide_seq:
            if p < 'A' or p > 'Z': continue
            if p == self.conf.beta_uua: # the mass of U and linker
                mass += self.conf.beta_uua_mass
            else:
                mass += aa[int(ord(p)-ord('A'))]
        for m in mods:
            mass += m.mass
        return mass

    def split_protein_non(self, delete_u_mass=False):
        peptide_seq_list_target, peptide_seq_list_decoy = self._non_specific_split()
        print("[Info]#Peptides in beta protein (non-spec) is {}".format(len(peptide_seq_list_target)))
        #for p in peptide_seq_list_target:
        #    print("[Info]Beta peptide (non-spec):", p)

        peptide_list = []
        for p in peptide_seq_list_target:
            mods = []  # don't support the modifications in beta peptide
            mass = self._get_peptide_mass(p, mods)
            if delete_u_mass:
                mass -= self.conf.beta_uua_mass # delete
            pep = CPeptide(p, 0, mass, mods)
            peptide_list.append(pep)
        for p in peptide_seq_list_decoy:
            mods = []  # don't support the modifications in beta peptide
            mass = self._get_peptide_mass(p, mods)
            if delete_u_mass:
                mass -= self.conf.beta_uua_mass # delete
            pep = CPeptide(p, -1, mass, mods)
            peptide_list.append(pep)
        return peptide_list

    def _generate_fix_one_pep(self, protein_len, peptide_seq, peptide_pos=None):
        conf = self.conf
        fix_aa = conf.beta_fix_aa
        fix_N_aa, fix_C_aa = conf.beta_fix_N_aa, conf.beta_fix_C_aa
        fix_PN_aa, fix_PC_aa = conf.beta_fix_PN_aa, conf.beta_fix_PC_aa
        fix_mods = []
        if peptide_pos != None:
            if len(fix_N_aa) > 0:
                if peptide_seq[0] in fix_N_aa:
                    mod = fix_N_aa[peptide_seq[0]]
                    modSite = CModSite(mod.name, 0, mod.mass)
                    fix_mods.append(modSite)
            if len(fix_PN_aa) > 0:
                if peptide_pos[0] == 0 and peptide_seq[0] in fix_PN_aa:
                    mod = fix_PN_aa[peptide_seq[0]]
                    modSite = CModSite(mod.name, 0, mod.mass)
                    fix_mods.append(modSite)
            if len(fix_C_aa) > 0:
                if peptide_seq[-1] in fix_C_aa:
                    mod = fix_C_aa[peptide_seq[-1]]
                    modSite = CModSite(mod.name, len(peptide_seq) + 1, mod.mass)
                    fix_mods.append(modSite)
            if len(fix_PC_aa) > 0:
                if peptide_pos[1] == protein_len and peptide_seq[-1] in fix_PC_aa:
                    mod = fix_PC_aa[peptide_seq[-1]]
                    modSite = CModSite(mod.name, len(peptide_seq) + 1, mod.mass)
                    fix_mods.append(modSite)
        if len(fix_aa) > 0:
            for i in range(len(peptide_seq)):
                a = peptide_seq[i]
                if a not in fix_aa: continue
                mod = fix_aa[a]
                modSite = CModSite(mod.name, i + 1, mod.mass)
                fix_mods.append(modSite)
        mass = self._get_peptide_mass(peptide_seq, fix_mods)
        peptide_res = CPeptide(peptide_seq, -1, mass, fix_mods)
        return peptide_res

    def split_protein(self):
        if self.protein_seq.find(';') >= 0:
            peptide_list = []
            for p in self.protein_seq.split(';'):
                #mods = []
                #mass = self._get_peptide_mass(p, mods)
                #pep = CPeptide(p, 0, mass, mods)
                pep = self._generate_fix_one_pep(len(self.protein_seq), p)
                peptide_list.append(pep)
            print("[Info]#Peptides in beta protein is {}".format(len(peptide_list)))
            for p in peptide_list:
                print("[Info]Beta peptide:", p.sq)
            return peptide_list
        split_sites = self._find_split_site()
        if self.enzyme_type == 0:
            peptide_seq_list, peptide_pos_list = self._specific_split(split_sites)
        elif self.enzyme_type == 1:
            peptide_seq_list, peptide_pos_list = self._semi_specific_split(split_sites)
        else: #need modify
            peptide_seq_list, peptide_pos_list = [], []
        #    peptide_seq_list = self._non_specific_split(split_sites)

        print("[Info]#Peptides in beta protein is {}".format(len(peptide_seq_list)))
        for p in peptide_seq_list:
            print("[Info]Beta peptide:", p)
        peptide_list = []
        for pi,p in enumerate(peptide_seq_list):
            #mods = [] # don't support the modifications in beta peptide
            #mass = self._get_peptide_mass(p, mods)
            #pep = CPeptide(p, 0, mass, mods)
            pep = self._generate_fix_one_pep(len(self.protein_seq), p, peptide_pos_list[pi])
            peptide_list.append(pep)
        return peptide_list

