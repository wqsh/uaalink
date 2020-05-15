#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: MSTask.py
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
from MSFunction import create_peptide_or_spectrum_index, CEnzymeFunction, CProteinFunction, CPickleFunction, CSpectrumFunction, CBetaProteinFunction, COutputFileFunction
from MSResultFunction import CResultFunction
from MSTool import global_match_for_one_link_peptide
from MSData import CPSM

class CCreateIndexTask():
    def __init__(self):
        pass

    def create(self, conf):
        self._load_and_process_protein(conf)

    def _load_and_process_protein(self, conf):
        enzymeFunc = CEnzymeFunction(conf)
        proteinFunc = CProteinFunction(enzymeFunc)
        protein_list = proteinFunc.load_protein(conf.fasta_file)
        print("[Info]#Proteins", len(protein_list))
        proteinFunc.split_protein(protein_list, False)

        # Save proteins and modifications in to pkl files
        protein_pkl_file = os.path.join(conf.index_out_folder, conf.PRO_PKL_FILE)
        modification_pkl_file = os.path.join(conf.index_out_folder, conf.MOD_PKL_FILE)
        pickleFunc = CPickleFunction()
        pickleFunc.save_pro_to_pkl(protein_list, protein_pkl_file)
        pickleFunc.save_mod_to_pkl(conf.mod_list, modification_pkl_file)


class CSearchTask():
    def __init__(self):
        pass

    def search(self, conf):
        # if new format of spectrum file, add "IF" and new methods
        self._search_mgf(conf)

    def _split_spectra(self, spec_list, spec_num):
        #split spectrum
        spec_list_new = []
        cur_spec = []
        for i in range(len(spec_list)):
            cur_spec.append(spec_list[i])
            if len(cur_spec) >= spec_num:
                spec_list_new.append(cur_spec)
                cur_spec = []
        if len(cur_spec) > 0: spec_list_new.append(cur_spec)
        return spec_list_new

    def _search_mgf(self, conf):
        beta_peptide_list = self._process_beta_protein(conf)
        mgf_file_list = conf.mgf_file_list
        spectrumFunc = CSpectrumFunction(conf)
        pickleFunc = CPickleFunction()
        outFunc = COutputFileFunction()
        for mgf_file in mgf_file_list:
            spectrum_list = spectrumFunc.load_spectrum(mgf_file)
            #if len(spectrum_list) <= conf.spectrum_num:
            #    spectrum_list_list = [spectrum_list]
            #else:
            #    spectrum_list_list = self._split_spectra(spectrum_list, conf.spectrum_num)
            #for spectrum_list in spectrum_list_list:
            cmpfun = operator.attrgetter('mass')
            spectrum_list.sort(key=cmpfun)
            print("[Info]#Spectra in {0} is {1}".format(mgf_file, len(spectrum_list)))
            spectrum_result_list = [[] for i in range(len(spectrum_list))]
            for beta_peptide in beta_peptide_list:
                s_list = self._preprocess_data(spectrum_list, beta_peptide, conf)
                for i,r in enumerate(s_list):
                    if len(spectrum_result_list[i]) == 0:
                        spectrum_result_list[i] = r
                    else:
                        spectrum_result_list[i] += r
                        cmpfun = operator.attrgetter('score')
                        spectrum_result_list[i].sort(key=cmpfun, reverse=True)
                        spectrum_result_list[i] = spectrum_result_list[i][:conf.top_k_result]
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
                result_file_name = mgf_file[mgf_file.rfind('/')+1:]
            else:
                result_file_name = mgf_file[mgf_file.rfind('\\')+1:]
            result_file_name = result_file_name[:result_file_name.find('.')]
            pkl_file_name = result_file_name + ".pkl"
            txt_file_name = result_file_name + ".txt"
            can_txt_file_name = result_file_name + "_cand.txt"
            print("[Info]Results for {0} are saved to {1}.".format(mgf_file, txt_file_name))
            #pickleFunc.save_psm_to_pkl(psm_list, os.path.join(conf.result_folder, pkl_file_name))
            outFunc.write_psms(psm_list, os.path.join(conf.result_folder, txt_file_name), conf)
            outFunc.write_psms(psm_list, os.path.join(conf.result_folder, can_txt_file_name), conf, True)
            outFunc.write_one_label(psm_list, mgf_file, conf.result_folder, conf)

    def _process_beta_protein(self, conf):
        betaProteinFunc = CBetaProteinFunction(conf)
        return betaProteinFunc.split_protein()

    def _preprocess_data(self, spectrum_list, beta_peptide, conf):
        # According to the mass of beta peptide (X)
        # Load all peptide pkl files and peptide index pkl files ([400, 500Da], [500,600Da], ..., [9900, 10000Da])
        # PKL file contains alpha peptide whose mass is P
        # Get the mass of xlink peptide is X+P
        # Then split the spectrum into different files
        # return the candidate link peptides for each spectrum (the length of result is equal to the length of spectrum_list)
        multi = conf.multi
        beta_mass = beta_peptide.mass
        pickleFunc = CPickleFunction()
        pkl_file_list, ind_file_list = pickleFunc.get_all_pkl_files(conf.index_out_folder, conf.PRO_PKL_FILE, conf.MOD_PKL_FILE)
        # Create the spectrum mass index
        spec_index = create_peptide_or_spectrum_index(spectrum_list, 0, spectrum_list[-1].mass, multi)

        input_data_list, spec_start_index = [], []
        for i in range(len(pkl_file_list)):
            alpha_start_mass, alpha_end_mass = pickleFunc.get_start_end_mass(pkl_file_list[i])
            spec_start_mass, spec_end_mass = alpha_start_mass + beta_mass + conf.ATOM_MASS_P, alpha_end_mass + beta_mass + conf.ATOM_MASS_P
            if conf.is_ppm_pre:
                spec_start_mass -= (spec_start_mass * conf.tol_pre)
                spec_end_mass += (spec_end_mass * conf.tol_pre)
            else:
                spec_start_mass -= conf.tol_pre
                spec_end_mass += conf.tol_pre
            if spec_start_mass > spectrum_list[-1].mass: continue
            if spec_start_mass < 0: spec_start_mass = 0.0
            if spec_end_mass > spectrum_list[-1].mass: spec_end_mass = spectrum_list[-1].mass
            s_start_ind = spec_index[int(spec_start_mass * multi)]
            while s_start_ind < len(spectrum_list) and spectrum_list[s_start_ind].mass < spec_start_mass:
                s_start_ind += 1
            s_end_ind = spec_index[int(spec_end_mass * multi)]
            while s_end_ind < len(spectrum_list) and spectrum_list[s_end_ind].mass <= spec_end_mass:
                s_end_ind += 1

            one_spectrum_list = spectrum_list[s_start_ind:s_end_ind]
            if len(one_spectrum_list) == 0: continue
            # one_spectrum_list is the list of spectra whose masses are in [spec_start_mass, spec_end_mass]
            input_data_list.append((one_spectrum_list, pkl_file_list[i], ind_file_list[i], beta_peptide, conf))
            spec_start_index.append(s_start_ind)
        #for input_data in input_data_list:
        #    global_match_for_one_link_peptide(input_data)
        if conf.thread_type == 0:
            pool = Pool(processes=conf.thread_num)
        else:
            pool = ThreadPool(conf.thread_num)
        res = pool.map(global_match_for_one_link_peptide, input_data_list)
        pool.close()
        pool.join()

        final_res = [[] for i in range(len(spectrum_list))]
        assert(len(res) == len(spec_start_index))
        for i, r in enumerate(res):
            ssi = spec_start_index[i]
            for j,spec_r in enumerate(r):
                if len(final_res[ssi+j]) == 0:
                    final_res[ssi+j] = spec_r
                else:
                    final_res[ssi+j] += spec_r
                    cmpfun = operator.attrgetter('score')
                    final_res[ssi+j].sort(key=cmpfun, reverse=True)
                    final_res[ssi+j] = final_res[ssi+j][:conf.top_k_result]

        return final_res

class CComputeFDRTask:
    def __init__(self):
        pass

    def _load_result(self, path):
        resFunc = CResultFunction(self.conf)
        psm_cand_list = resFunc.load_psm(path)
        psm_list = []
        title_dict = {}
        for (title, pep, decoy, ranker, delta_score) in psm_cand_list:
            if title not in title_dict: title_dict[title] = []
            title_dict[title].append((pep, decoy, ranker))
        for title in title_dict:
            for (pep, decoy, ranker) in title_dict[title]:
                if ranker == 1:
                    if not decoy: # target
                        psm_list.append((title, pep, decoy, ranker))
                    else: # decoy, check the other peps with same score
                        cand_list = title_dict[title]
                        cur_decoy = True
                        for (pep2, decoy2, ranker2) in cand_list:
                            if pep2.score == pep.score and not decoy2:
                                cur_decoy = False
                                break
                        psm_list.append((title, pep, cur_decoy, ranker))
        return psm_list

    def _compute_fdr_score(self, conf):
        self.conf = conf
        res_fun = CResultFunction(conf)
        cand_file_list = res_fun.find_res_file(conf)
        print("[Info] Candidate file list", cand_file_list)
        psm_list = []
        for file_path in cand_file_list:
            psm_list += self._load_result(file_path)
        psm_list.sort(key=lambda k:k[1].score, reverse=True)
        target_num, decoy_num, all_num = 0, 0, len(psm_list)
        for (title, pep, decoy, ranker) in psm_list:
            if decoy: decoy_num += 1
            else: target_num += 1
        print("[Info] #PSMs is {0}, #Target is {1} and #Decoy is {2}".format(all_num, target_num, decoy_num))
        cur_target_num, cur_decoy_num = 0, 0
        if len(psm_list) == 0: return 0.0
        score_t = psm_list[0][1].score + 1.0
        if target_num > 0 and decoy_num * 1.0 / target_num <= conf.fdr_value:
            score_t = psm_list[-1][1].score
            return score_t
        for (title, pep, decoy, ranker) in psm_list[::-1]:
            if decoy: cur_decoy_num += 1
            else: cur_target_num += 1
            if target_num == cur_target_num: 
                fdr_value = 100
            else:
                fdr_value = (decoy_num - cur_decoy_num) * 1.0 / (target_num - cur_target_num)
            #print(decoy_num,cur_decoy_num, target_num, cur_target_num, fdr_value)
            if fdr_value <= conf.fdr_value:
                score_t = pep.score
                return score_t
        return score_t

    def _write_fdr_result(self, conf, fdr_score):
        resFunc = CResultFunction(conf)
        taskid = resFunc.get_taskid(conf, False)
        res_file_list = resFunc.find_res_file(conf)
        res_plabel_file_list = resFunc.find_res_file(conf, label_flag = True)

        pep_score, pep_line = {}, {}
        t_p_line = {}
        head_line = ""
        title_score, title_line, title_plabel_line = {}, {}, {}
        pep_plabel_line = {}
        for file in res_file_list:
            with open(file) as f:
                lines = f.readlines()
            head_line = lines[0].strip()
            for i in range(1, len(lines)):
                line = lines[i].strip()
                tmp = line.split('\t')
                title = tmp[0]
                sq_key = tmp[2] + '@' + tmp[3] + '@' + tmp[4] + '@' + tmp[6]
                ac_list = tmp[11].split(';')
                is_target = False
                for ac in ac_list:
                    if not ac.startswith("REV_"):
                        is_target = True
                        break
                if not is_target: continue
                s = float(tmp[12])
                if s <= fdr_score: break
                title_score[title] = s
                title_line[title] = line
                if sq_key not in pep_score or pep_score[sq_key] < s:
                    pep_score[sq_key] = s
                    pep_line[sq_key] = line

        for sq in pep_line:
            line = pep_line[sq]
            tmp = line.strip().split('\t')
            title = tmp[0]
            t_p_line[title] = 1
        
        modification_lines = []
        for file in res_plabel_file_list:
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
                    modification_lines = []
                    for j in range(i+1, len(lines)):
                        ind = lines[j].find('=')
                        if ind < 0: break
                        modification_lines.append(lines[j].strip())

        sorted_res = sorted(title_score.items(), key=lambda d: d[1], reverse=True)
        write_path = os.path.join(conf.result_folder, conf.XLINK_PSM + '-' + taskid + '.txt')
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r,s) in sorted_res:
                one_line = title_line[r]
                f.write(one_line + '\n')
        write_path = os.path.join(conf.result_folder, conf.XLINK_PSM + '-' + taskid + '.plabel')
        with open(write_path, 'w') as f:
            f.write("[FilePath]\n")
            f.write("File_Path=" + ";".join(conf.mgf_file_list) + "\n")
            f.write("[Modification]\n")
            for mod_line in modification_lines:
                f.write(mod_line + '\n')
            f.write("[xlink]\n")
            f.write("xlink=UAA\n")
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
        write_path = os.path.join(conf.result_folder, conf.XLINK_PEP + '-' + taskid + '.txt')
        with open(write_path, 'w') as f:
            f.write(head_line + '\n')
            for (r, s) in sorted_res:
                one_line = pep_line[r]
                f.write(one_line + '\n')
                
        write_path = os.path.join(conf.result_folder, conf.XLINK_PEP + '-' + taskid + '.plabel')
        with open(write_path, 'w') as f:
            f.write("[FilePath]\n")
            f.write("File_Path=" + ";".join(conf.mgf_file_list) + "\n")
            f.write("[Modification]\n")
            for mod_line in modification_lines:
                f.write(mod_line + '\n')
            f.write("[xlink]\n")
            f.write("xlink=UAA\n")
            f.write("[Total]\n")
            f.write("total=" + str(len(sorted_res)) + '\n')
            spec_id = 0
            for (r, s) in sorted_res:
                spec_id += 1
                one_line = pep_line[r]
                tmp = one_line.strip().split('\t')
                title = tmp[0]
                f.write("[Spectrum" + str(spec_id) + "]\n")
                f.write("name=" + title + "\n")
                one_line = pep_plabel_line[title]
                f.write(one_line + "\n")

    def compute_fdr(self, conf):
        fdr_score = self._compute_fdr_score(conf)
        print("[Info] Score for fdr<={0} is {1}".format(conf.fdr_value, fdr_score))
        self._write_fdr_result(conf, fdr_score)

class CRerankTask:
    def __init__(self):
        pass

    def _write_svm_file(self, psm_list, write_path):
        fw = open(write_path, 'w')
        for (title, link_pep, decoy_flag, ranker, delta_score) in psm_list:
            alpha_score, beta_score = link_pep.alpha_score, link_pep.beta_score
            len1 = len(link_pep.alpha_peptide.sq)
            len2 = len(link_pep.beta_peptide.sq)
            class_flag = 1
            if decoy_flag: class_flag = -1
            fid_start, fid_start2 = 4, 13
            alpha_f_str = alpha_score.get_string(" ", True, fid_start)
            #beta_f_str = beta_score.get_string(" ", True, fid_start2)
            fw.write(str(class_flag) + " 1:" + str(link_pep.score) + " 2:" + str(delta_score) + " 3:" + str(len1) + " " + alpha_f_str + "\n")
        fw.close()

    def _svm_classify(self, test_path, model_path, output_path):
        os.system("./libsvm/svm-predict -b 1 " + test_path + " " + model_path + " " + output_path)

    def _get_svm_score(self, score_path):
        score_list = []
        pos_flag = 1
        with open(score_path) as f:
            for line in f:
                line = line.strip()
                if line == "": continue
                if line.startswith("labels"):
                    tmp = line.split(' ')
                    if tmp[2] == '1': pos_flag = 2
                    continue
                tmp = line.split(' ')
                score_list.append(str(tmp[pos_flag]))
        return score_list 

    def _update_psm_score(self, result_file_list, score_list):
        si = 0
        for file_path in result_file_list:
            write_path = file_path + "_new.txt"
            tmp_list = []
            fw = open(write_path, 'w')
            head_line = ""
            with open(file_path) as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    if line_num <= 1: 
                        head_line = line.strip()
                        continue
                    tmp = line.strip().split('\t')
                    tmp[12] = score_list[si]
                    si += 1
                    tmp_list.append((tmp[12], "\t".join(tmp)))
                    #fw.write("\t".join(tmp) + '\n')
            tmp_list.sort(key=lambda tup: tup[0], reverse=True)
            fw.write(head_line + '\n')
            for (_, line) in tmp_list:
                fw.write(line + '\n')
            fw.close()

    def _rerank_result(self, conf):
        res_fun = CResultFunction(conf)
        result_file_list = res_fun.find_res_file(conf)
        psm_list = []
        for file_path in result_file_list:
            psm_list += res_fun.load_psm(file_path)
        feature_write_path = os.path.join(conf.INPUT_SVM_FOLDER, res_fun.get_taskid(conf))
        print('aa', conf.INPUT_SVM_FOLDER, feature_write_path)
        self._write_svm_file(psm_list, feature_write_path)
        score_write_path = feature_write_path + ".txt"
        self._svm_classify(feature_write_path, conf.MODEL_PATH, score_write_path)
        #os.remove(feature_write_path)
        score_list = self._get_svm_score(score_write_path) 
        assert(len(score_list) == len(psm_list))
        self._update_psm_score(result_file_list, score_list)


    def rerank_result(self, conf):
        self._rerank_result(conf)

class CDrawTask:
    def __init__(self):
        pass

    def draw(self, conf):
        pass
