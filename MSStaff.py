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
import datetime
from MSSystem import INFO_TO_USER_Staff, EXPIRATION_TIME
import time
from Config import CConfig
from MSFunction import CConfigFunction
from MSFlow import CFlowCreateSearch, CFlowCreate, CFlowSearch, CFlowComputeFDR, CFlowRerank

class CStaff:
    def __init__(self):
        pass

    def _checkTime(self):
        dateNow = datetime.datetime.now()
        dateDead = datetime.datetime(EXPIRATION_TIME['Year'], EXPIRATION_TIME['Month'], EXPIRATION_TIME['Day'], 23, 59)

        n_day = (dateDead - dateNow).days

        if n_day < 0:
            print(INFO_TO_USER_Staff[0])
            exit(-1)

        elif n_day < 7:
            print(INFO_TO_USER_Staff[1])

        else:
            print(INFO_TO_USER_Staff[2])
            print(dateNow)
            
    def start(self, argv):
        if len(argv) == 1:
            start_time = time.time()
            conf = CConfig()
            confFun = CConfigFunction()
            confFun.config2file("./configure.txt", conf)
            end_time = time.time()
            print("[Info] Task finish in {}s".format(end_time - start_time))
        elif len(argv) >= 2:
            start_time = time.time()
            conf = CConfig()
            confFun = CConfigFunction()
            print("Config file", argv)
            confFun.file2config(argv[1], conf)
            if conf.task_type == 0:
                flow = CFlowCreateSearch(conf)
                flow.run()
            elif conf.task_type == 1:
                flow = CFlowCreate(conf)
                flow.run()
            elif conf.task_type == 2:
                flow = CFlowSearch(conf)
                flow.run()
            elif conf.task_type == 3:
                flow = CFlowRerank(conf)
                flow.run()
            elif conf.task_type == 4:
                flow = CFlowComputeFDR(conf)
                flow.run()
            end_time = time.time()
            print("[Info] Task finish in {}s".format(end_time - start_time))
