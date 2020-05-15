#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: MSFlow.py
Author: webdev@tencent.com
Date: 2019/05/24 10:05:01
Brief: 
"""

import sys
#reload(sys)

from MSTask import CCreateIndexTask, CSearchTask, CDrawTask, CComputeFDRTask, CRerankTask
from MSSingleSearchTask import CSingleSearchTask, CSingleComputeFDRTask

class CFlowCreateSearch:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
        createTask = CCreateIndexTask()
        createTask.create(self.conf)
        if self.conf.need_open_search != 0:
            singleSearchTask = CSingleSearchTask()
            singleSearchTask.search(self.conf)
            singleComputeTask = CSingleComputeFDRTask()
            singleComputeTask.compute_fdr(self.conf)
        searchTask = CSearchTask()
        searchTask.search(self.conf)
        #rerankTask = CRerankTask()
        #rerankTask.rerank_result(self.conf)
        computeTask = CComputeFDRTask()
        computeTask.compute_fdr(self.conf)

class CFlowCreate:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
        createTask = CCreateIndexTask()
        createTask.create(self.conf)


class CFlowSearch:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
        if self.conf.need_open_search != 0:
            singleSearchTask = CSingleSearchTask()
            singleSearchTask.search(self.conf)
            singleComputeTask = CSingleComputeFDRTask()
            singleComputeTask.compute_fdr(self.conf)
        searchTask = CSearchTask()
        searchTask.search(self.conf)
        #rerankTask = CRerankTask()
        #rerankTask.rerank_result(self.conf)
        computeTask = CComputeFDRTask()
        computeTask.compute_fdr(self.conf)

class CFlowRerank:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
        rerankTask = CRerankTask()
        rerankTask.rerank_result(self.conf)
        computeTask = CComputeFDRTask()
        computeTask.compute_fdr(self.conf)

class CFlowComputeFDR:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
        computeTask = CComputeFDRTask()
        computeTask.compute_fdr(self.conf)

class CFlowDraw:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
        drawTask = CDrawTask()
        drawTask.draw(self.conf)
