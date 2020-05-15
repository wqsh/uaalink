#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: MSSingleSearchTask.py
Author: webdev@tencent.com
Date: 2019/05/24 10:05:37
Brief:
"""


import sys
#reload(sys)
import os
import operator
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from MSFunction import create_peptide_or_spectrum_index, CBetaProteinFunction, CSpectrumFunction, CPickleFunction, COutputFileFunction
from MSResultFunction import CResultFunction
from MSTool import global_match_for_single_peptide
from MSData import CPSM

class CSingleSearchTask:
    def __init__(self):
        pass

    def search(self, conf):
        self._search_mgf(conf)

    def _search_mgf(self, conf):
        beta_peptide_list = self._process_beta_protein(conf)
        mgf_file_list = conf.mgf_file_list
        spectrumFunc = CSpectrumFunction(conf)
        pickleFunc = CPickleFunction()
        outFunc = COutputFileFunction()
        for mgf_file in mgf_file_list:
            spectrum_list = spectrumFunc.load_spectrum(mgf_file)
            cmpfun = operator.attrgetter('mass')
            spectrum_list.sort(key=cmpfun)
            print("[Info]#Spectra in {0} is {1}".format(mgf_file, len(spectrum_list)))
            spectrum_result_list = [[] for i in range(len(spectrum_list))]

            spec_index_whole = create_peptide_or_spectrum_index(spectrum_list, 0, spectrum_list[-1].mass, conf.multi)
            for beta_index, beta_peptide in enumerate(beta_peptide_list):
                s_list, spec_start_ind, spec_end_ind = self._preprocess_data(spectrum_list, spec_index_whole, beta_peptide, conf)
                for i, r in enumerate(s_list):
                    spec_index = spec_start_ind + i
                    spectrum_result_list[spec_index].append(r)
                    if len(spectrum_result_list[spec_index]) > 1:
                        cmpfun = operator.attrgetter('score')
                        spectrum_result_list[spec_index].sort(key=cmpfun, reverse=True)
                        spectrum_result_list[spec_index] = spectrum_result_list[spec_index][:conf.top_k_result]

            psm_list = []
            for i in range(len(spectrum_result_list)):
                p_list = spectrum_result_list[i]
                if len(p_list) == 0: continue
                best_p = p_list[0]
                cur_spec = spectrum_list[i]
                cur_spec.peaks = []
                psm_list.append(CPSM(cur_spec, best_p, p_list))
            cmpfun = operator.attrgetter('peptide.score')
            psm_list.sort(key=cmpfun, reverse=True)
            if mgf_file.find('/') >= 0:
                result_file_name = mgf_file[mgf_file.rfind('/') + 1:]
            else:
                result_file_name = mgf_file[mgf_file.rfind('\\') + 1:]
            result_file_name = result_file_name[:result_file_name.find('.')]
            txt_file_name = result_file_name + "_singlePeptide.txt"
            print("[Info]Results for {0} are saved to {1}.".format(mgf_file, txt_file_name))
            outFunc.write_single_psms(psm_list, os.path.join(conf.result_folder, txt_file_name), conf)
            #outFunc.write_one_label(psm_list, mgf_file, conf.result_folder, conf)

    def _preprocess_data(self, spectrum_list, spec_index, beta_peptide, conf):
        #According to the mass of beta peptide, find the spectra whose mass is within [-200 Da, 200 Da]
        beta_mass = beta_peptide.mass
        spec_start_mass, spec_end_mass = beta_mass + conf.ATOM_MASS_P - conf.window_size, beta_mass + conf.ATOM_MASS_P + conf.window_size
        if spec_start_mass > spectrum_list[-1].mass: return [], 0, -1
        if spec_start_mass < 0: spec_start_mass = 0.0
        if spec_end_mass > spectrum_list[-1].mass: spec_end_mass = spectrum_list[-1].mass

        multi = conf.multi
        s_start_ind = spec_index[int(spec_start_mass * multi)]
        while s_start_ind < len(spectrum_list) and spectrum_list[s_start_ind].mass < spec_start_mass:
            s_start_ind += 1
        s_end_ind = spec_index[int(spec_end_mass * multi)]
        while s_end_ind < len(spectrum_list) and spectrum_list[s_end_ind].mass <= spec_end_mass:
            s_end_ind += 1

        input_data_list = []
        for spec_index in range(s_start_ind, s_end_ind):
            input_data_list.append((spectrum_list[spec_index], beta_peptide, conf))

        #for one_data in input_data_list:
        #    global_match_for_single_peptide(one_data)


        if conf.thread_type == 0:
            pool = Pool(processes=conf.thread_num)
        else:
            pool = ThreadPool(conf.thread_num)
        res = pool.map(global_match_for_single_peptide, input_data_list)
        pool.close()
        pool.join()

        return res, s_start_ind, s_end_ind

    def _process_beta_protein(self, conf):
        betaProteinFunc = CBetaProteinFunction(conf)
        return betaProteinFunc.split_protein_non(delete_u_mass=True)

class CSingleComputeFDRTask:
    def __init__(self):
        pass

    def _load_result(self, path):
        resFunc = CResultFunction(self.conf)
        psm_list = resFunc.load_single_psm(path)
        return psm_list

    def _compute_fdr_score(self):
        conf = self.conf
        res_fun = CResultFunction(conf)
        cand_file_list = res_fun.find_single_res_file(conf)
        print("[Info] Candidate file list", cand_file_list)
        psm_list = []
        for file_path in cand_file_list:
            psm_list += self._load_result(file_path)
        psm_list.sort(key=lambda k: k[1].score, reverse=True) #按照score排序
        target_num, decoy_num, all_num = 0, 0, len(psm_list)
        for (title, pep, decoy, ranker, delta_score) in psm_list:
            if decoy:
                decoy_num += 1
            else:
                target_num += 1
        print("[Info] #PSMs is {0}, #Target is {1} and #Decoy is {2}".format(all_num, target_num, decoy_num))
        cur_target_num, cur_decoy_num = 0, 0
        if len(psm_list) == 0: return 0.0
        score_t = psm_list[0][1].score + 1.0
        if target_num > 0 and decoy_num * 1.0 / target_num <= conf.fdr_value:
            score_t = psm_list[-1][1].score
            return score_t
        for (title, pep, decoy, ranker, delta_score) in psm_list[::-1]:
            if decoy:
                cur_decoy_num += 1
            else:
                cur_target_num += 1
            if target_num == cur_target_num:
                fdr_value = 100
            else:
                fdr_value = (decoy_num - cur_decoy_num) * 1.0 / (target_num - cur_target_num)
            # print(decoy_num,cur_decoy_num, target_num, cur_target_num, fdr_value)
            if fdr_value <= conf.fdr_value:
                score_t = pep.score
                return score_t
        return score_t

    def _write_fdr_result(self, conf, fdr_score):
        resFunc = CResultFunction(conf)
        taskid = resFunc.get_taskid(conf, False)
        res_file_list = resFunc.find_single_res_file(conf)

        pep_score, pep_line = {}, {}
        head_line = ""
        title_score, title_line = {}, {}
        for file in res_file_list:
            with open(file) as f:
                lines = f.readlines()
            head_line = lines[0].strip()
            for i in range(1, len(lines)):
                line = lines[i].strip()
                tmp = line.split('\t')
                title = tmp[0]
                sq_key = tmp[2] + '@' + tmp[4] + '@' + str(int(float(tmp[3]) + 0.5))
                is_target = (tmp[8] == "True")
                if not is_target: continue
                s = float(tmp[5])
                if s <= fdr_score: break
                title_score[title] = s
                title_line[title] = line
                if sq_key not in pep_score or pep_score[sq_key] < s:
                    pep_score[sq_key] = s
                    pep_line[sq_key] = line

        sorted_res = sorted(title_score.items(), key=lambda d: d[1], reverse=True)
        write_path = os.path.join(conf.result_folder, conf.SINGLE_PSM + '-' + taskid + '.txt')
        conf.single_res_path = write_path
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r, s) in sorted_res:
                one_line = title_line[r]
                f.write(one_line + '\n')

        sorted_res = sorted(pep_score.items(), key=lambda d: d[1], reverse=True)
        write_path = os.path.join(conf.result_folder, conf.SINGLE_PEP + '-' + taskid + '.txt')
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r, s) in sorted_res:
                one_line = pep_line[r]
                f.write(one_line + '\n')

    def compute_fdr(self, conf):
        self.conf = conf
        fdr_score = self._compute_fdr_score()
        print("[Info] Score for fdr<={0} is {1}".format(conf.fdr_value, fdr_score))
        self._write_fdr_result(conf, fdr_score)