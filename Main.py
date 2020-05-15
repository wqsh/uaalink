#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: CreateIndex.py
Author: webdev@tencent.com
Date: 2019/05/21 18:05:47
Brief: 
"""

import sys
#reload(sys)
from MSStaff import CStaff
import multiprocessing

if __name__ == "__main__":
    multiprocessing.freeze_support()
    staff = CStaff()
    staff._checkTime()
    staff.start(sys.argv)
