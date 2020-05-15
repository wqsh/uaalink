#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: Config.py
Author: webdev@tencent.com
Date: 2019/05/21 17:05:51
Brief: 
"""

import sys
#reload(sys)
import os
from MSData import CElement, CAA, CModification
import math

class CConfig:
    def __init__(self):
        self.ATOM_MASS_C = 12.0
        self.ATOM_MASS_H = 1.007825035
        self.ATOM_MASS_P = 1.00727647012
        self.ATOM_MASS_O = 15.99491463
        self.ATOM_MASS_N = 14.003074
        self.ATOM_MASS_C13 = 13.00335
        self.ATOM_MASS_S = 31.972072
        self.MOLECULE_MASS_H2O = 18.0105647
        self.MOLECULE_MASS_NH3 = 17.026549105
        self.AA_NUM = 26 # 26 amino acids
        
        self.AVRG_AA_MASS = 110.0

        self.PRO_PKL_FILE = "protein.pkl"
        self.MOD_PKL_FILE = "modification.pkl"

        self.INPUT_SVM_FOLDER = "./svm_tool/input_data/"
        self.MODEL_PATH = "./svm_tool/models/v1"
        
        self.SINGLE_PSM = "res_single_psm"
        self.SINGLE_PEP = "res_single"
        self.XLINK_PSM = "res_psm"
        self.XLINK_PEP = "res"
        self.single_res_path = ""
        
        #initial configure default value in configure.txt
        self.ms2_type = "mgf"
        self.mgf_file_list_str = " #VIP"
        self.fasta_file = " #VIP"
        self.index_out_folder = "./protein_index/ #VIP"
        self.result_folder = "./result/ #VIP"
        self.enzyme_str = "trypsin KR C #use ';' to set multiply enzymes"
        self.enzyme_type = "0 #0 for specific; 1 for semi-specific; 2 for non-specific" 
        self.enzyme_miss_cle_num = "3"
        self.fix_mod_str = " #VIP, use ';' to set multiply fixed modifications"
        self.var_mod_str = "Oxidation[M] #VIP, use ';' to set multiply variable modifications"
        self.var_mod_site_num = "3 #Maximum of variable modification in one peptide sequence (not consider the fixed modifications)"
        self.protein_beta = " #VIP"
        self.beta_uua = "U"
        self.beta_pep_len_min = "4"
        self.beta_pep_len_max = "20"
        self.beta_fix_mod_str = ""
        self.beta_var_mod_str = ""
        self.beta_uua_str = " #VIP"
        self.enzyme_str_beta = "trypsin KR C # use ';' to set multiply enzymes"
        self.enzyme_type_beta = "0 #0 for specific; 1 for semi-specific; 2 for non-specific"
        self.enzyme_miss_cle_num_beta = "0"
        self.link_aa = " #VIP, beta peptide is linked with which amino acids in alpha peptide, it can be multiply amino acids (e.g., ACDEF)"
        self.tol_pre_str = "20ppm"
        self.tol_fra_str = "20ppm"
        self.activation_type = "HCD"
        self.thread_num = "8"
        self.thread_type = "0 #0 is for multi-process (high speed but use more memory); 1 is for multi-thread (low speed but use less memory)"
        self.select_peak_num = "200"
        self.spectrum_num = "10000"
        self.protein_len_max = "100000"
        self.pep_mass_min = "400"
        self.pep_mass_max = "10000"
        self.pep_len_min = "6"
        self.pep_len_max = "100"
        self.index_split_mass = "100 #create one pkl file for each 100Da ([0, 100], [100, 200], ..., [9900, 10000])"
        self.top_k_result = "10 #output top-10 peptides for each spectrum"
        self.multi = "1 #use mass hash (mass*MUTLI_MASS) to retrive peptide, spectrum or peak, this value of creating peptide index and searching mgf must be same."
        self.task_type = "0 #0 is create peptide index for (fasta), search mgf files, rerank results and computing fdr; 1 is only create peptide index; \
 2 is only search mgf files, rerank results and computing fdr; 3 is rerank results and computing fdr; 4 is only computing fdr"
        self.is_filter_beta = "1 #whether to filter spectrum which has not matched ion when matching with beta peptide (default is 1)"
        self.beta_filter_v = "1 #when 'TYPE_FILTER_BETA' is set as 1, then this value is valid, \
 it will filter spectrum where the number of matched beta ions is less than 'NUMBER_PEAK_BETA' (default is 1)"
        self.need_open_search = "0 # open search for single peptide (beta), 0 means it don't support open search while 1 means it supports"
        self.window_size = "300 #when open search only for single beta peptide, the mass window size of open search is 300 Da"
        self.pfind_result_path = " #path of pfind result file; if not exists, it can be empty"
        self.fdr_value = "0.05"
        self.element_file = "./ini/element.ini"
        self.aa_file = "./ini/aa.ini"
        self.modification_file = "./ini/modification.ini"