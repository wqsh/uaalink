#!/usr/bin/env python
# coding=utf-8

#########################################################################
#
# Copyright (c) 2017 Tencent Inc. All Rights Reserved
#
#########################################################################


"""
File: GetPeptideMass.py
Author: webdev@tencent.com
Date: 2019/05/27 10:05:57
Brief: 
"""

import sys
#reload(sys)
import re
from MSData import CElement, CAA


element = CElement("./ini/element.ini")
aa = CAA("./ini/aa.ini", element.element2mass, 26)

def get_mass(sq):
    is_elem = False
    if sq.find('(') >= 0: is_elem = True
    print("IS_elem", is_elem)
    if is_elem:
        tmp = re.split('\(|\)', sq) 
        mass = 0.0
        for i in range(len(tmp)):
            if i % 2 != 0: continue
            if tmp[i].strip() == "": continue
            if tmp[i] not in element.element2mass: continue
            mass += (float(element.element2mass[tmp[i]]) * int(tmp[i+1]))
        print("Mass", mass)
        return mass
    else:
        mass_H2O = 18.0105647
        mass = 0.0
        for p in sq:
            if p < 'A' or p > 'Z': continue
            mass += aa.aa2mass[int(ord(p)-ord('A'))]

        print("Mass:", mass)
        print("Mass + H2O:", mass + mass_H2O)
        return mass + mass_H2O
        
if __name__ == "__main__":
    sq = sys.argv[1]
    get_mass(sq)
