# -*- coding: utf-8 -*-
"""
@author: Inabahu(Inabahu@tju.edu.cn)
@file: read_POSCAR.py
@time: 2019/7/22
ref:https://github.com/lxf-gzu/vasp_small_script
    
"""

import os
import re
import sys
import os
from numpy import *
import numpy as np
import pprint

def read_poscar(filepath):
    file=os.path.join(filepath)
    lattice=[]
    poscar_tem=[]
    scale=[]
    if os.path.exists(file):
       poscar=open(file,'r')
    
    else:
        raise IOError('POSCAR does not exist！')
    for line in poscar.readlines():
                poscar_tem.append(line)
    poscar.close()
    #此步将POSCAR文件信息转入poscar_tem
    name=poscar_tem[0].split()
    scale=poscar_tem[1].split()
    lat=poscar_tem[2:5]
    for i in lat:
        lattice.append(i.split())
    #包含原子名
    atom=poscar_tem[5].split()
    #包含原子个数
    atom_number=poscar_tem[6].split()
    atom_position=[]
    sum_atom_number=sum(int(i) for i in atom_number)
    m=poscar_tem[8:8+sum_atom_number]
    for i in m:
        atom_position.append(i.split())

    scale_lat=[]
    scale_lat.append(scale)
    scale_lat.append(lattice)
    scale_lat.append(atom)
    scale_lat.append(atom_number)
    return name,scale_lat,atom_position,sum_atom_number
 #   return poscar_tem
#输出中第一个为晶体名
#输出第二个为晶体所有晶胞参数
#输出第三个为晶体所有原子位置
 
 
[name_a,lattice_a,atom_position_a,sum_atom_number_a]=read_poscar("C:/Users/neetsaki/Desktop/vaspkit_game/POSCAR/03_black_P_POSCAR")
[name_b,lattice_b,atom_position_b,sum_atom_number_b]=read_poscar("C:/Users/neetsaki/Desktop/vaspkit_game/POSCAR/03_Graphene_POSCAR")
print(lattice_b)
