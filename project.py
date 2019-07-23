# -*- coding: utf-8 -*-

"""
@author: Inabahu(Inabahu@tju.edu.cn)
@file: read_POSCAR.py
@time: 2019/7/22
ref:https://github.com/lxf-gzu/vasp_small_script
    
"""

import time
import re
import sys
import os
from numpy import *
import numpy as np
import pprint
import tkinter as tk
from tkinter import filedialog
 
print("****************************************************")
print("*              Author name:Inabahu                 *")
print("*          Email:Inabahu@tju.edu.cn                *")
print("*  Usage description:Heterostructure construction  *")
print("               Version: 0.1                        *")
print("           First released time: 7/23 2019          *")
print("****************************************************")
time.sleep(1)

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
#crystal_A="C:/Users/neetsaki/Desktop/vaspkit_game/POSCAR/02_Graphene_POSCAR"
#crystal_A=input("Please input POSCAR A file path:") 
#crystal_B="C:/Users/neetsaki/Desktop/vaspkit_game/POSCAR/02_MoS2_POSCAR"
#crystal_B=input("Please input POSCAR B file path:")
print("Please select the POSCAR files")
time.sleep(1)
root = tk.Tk()
root.withdraw() 
file_path = filedialog.askopenfilenames()

since = time.time()

[name_a,lattice_a,atom_position_a,sum_atom_number_a]=read_poscar(file_path[0])
[name_b,lattice_b,atom_position_b,sum_atom_number_b]=read_poscar(file_path[1])















import numpy as np


def judge_mismatch(a,b,threhold):
    if abs(a-b)/min(a,b)<threhold:
        x=True
    else:
        x=False
    return x

def mismatch(a,b):
    return abs(a-b)/min(a,b)

def fit_mismatch(a,b,threhold):
    j=1
    for i in range(1,200):
        if judge_mismatch(a*i,b*j,threhold)==True:
            break
        else:
            if a*i>b*j:
                j+=1
            else:
                j=j
    return i,j

#计算向量夹角
def calc_angle(a,b):
    a=np.array(a)
    b=np.array(b)
    La=np.sqrt(a.dot(a))
    Lb=np.sqrt(b.dot(b))
    cos_angle=a.dot(b)/(La*Lb)
    angle=np.arccos(cos_angle)
    return angle

def judge_arg_mismatch(a1,b1,a2,b2,threhold):
    angle1=calc_angle(a1,b1)
    angle2=calc_angle(a2,b2)
    return judge_mismatch(angle1,angle2,threhold)

def mismatch_matrix(s,t):
    MM=[[0,0,0],[0,0,0],[0,0,0]]
    for si in range(len(s)):
        for sj in range(len(s)):
            MM[si][sj]=mismatch(s[si],t[sj])
    return MM

#核算两晶体基矢夹角是否一致
def arg_mismatch(A,B):
    arg1=[calc_angle(A[0],A[1]),calc_angle(A[2],A[1]),calc_angle(A[0],A[2])]
    arg2=[calc_angle(B[0],B[1]),calc_angle(B[2],B[1]),calc_angle(B[0],B[2])]
    arg1=sorted(arg1)
    arg2=sorted(arg2)
    D = (np.array(arg1)-np.array(arg2)).dot(np.array(arg1)-np.array(arg2))
    return D

#核算两晶体a,b,c是否一致

#测试        
#[i,j]=fit_mismatch(3.1,5.2,0.003)
A=lattice_a[1]
B=lattice_b[1]
for i in range(0,len(A)):
    for j in range(0,len(A)):
        A[i][j]=float(A[i][j])
        B[i][j]=float(B[i][j])
#a=judge_arg_mismatch(A[0],A[1],A[1],A[2],0.003)
#b=judge_arg_mismatch(B[0],B[1],B[1],B[2],0.003)
#计算基矢夹角相似性
D = arg_mismatch(A,B)

LA=[]
LB=[]
for i in range(3):
    LA.append(np.sqrt(np.array(A[i]).dot(np.array(A[i]))))
    LB.append(np.sqrt(np.array(B[i]).dot(np.array(B[i]))))


L_threhold=0.05
print("*************************************************************")
print("* Please input the lattice mismatch threhold (default:0.05) *")
print("* !!!attention:overstrict standard will induce huge lattice *")
print("*************************************************************")
L_threhold=float(input("* Mismatch:"))

        #扩胞后对原子位置修正的函数
        #atom_position 以list形式输入
def position_corr(atom_position,m,n):
    atom_position_new=[]
    for i in range(m):
        for j in range(n):
            atom_position_new.append(np.multiply(np.array(atom_position).astype('float64'),np.array([1/m,1/n,1])+np.array([i/m,j/n,0])).tolist())
            
    return atom_position_new


if D<1e-10:
    print("****angle test pass****")
    MM=mismatch_matrix(LA[0:2],LB[0:2])
    #由于是二维材料，考虑a,b轴匹配即可,@@@@@@@@@@@@@此处遗留问题，如何断定两轴为a和b轴
    if max(max(MM))<L_threhold:
        print("****lattice test pass****")
        #接合并晶胞函数
    else:
        print("****lattice test fail, run mismatch fit****")
        [m1,n1]=fit_mismatch(LA[0],LB[0],L_threhold)
        [m2,n2]=fit_mismatch(LA[1],LB[1],L_threhold)
        m=[m1,m2]
        n=[n1,n2]
        A[0]=np.array(A[0])*m1
        A[0]=A[0].tolist()
        A[1]=np.array(A[1])*m2
        A[1]=A[1].tolist()
        B[0]=np.array(B[0])*n1
        B[0]=B[0].tolist()
        B[1]=np.array(B[1])*n2
        B[1]=B[1].tolist()
        atom_position_a_new=[]
        atom_position_b_new=[]
        for i in range(len(atom_position_a)):
            atom_position_a_new=atom_position_a_new+position_corr(atom_position_a[i],m1,m2)
        for i in range(len(atom_position_b)):
            atom_position_b_new=atom_position_b_new+position_corr(atom_position_b[i],n1,n2)
        atom_position_a=atom_position_a_new
        for i in range(len(lattice_a[3])):
            lattice_a[3][i]=str(int(lattice_a[3][i])*m1*m2)
        atom_position_b=atom_position_b_new
        for i in range(len(lattice_b[3])):
            lattice_b[3][i]=str(int(lattice_b[3][i])*n1*n2)
        
else:
    print("****angle test fail, need rotation****")
    #夹角不匹配怎么处理啊啊啊啊啊啊啊啊
























#需要实现 异质结结合(√)，异质结旋转
interlayer_distance=3
print('\n'+"********************************************************************")
print("* Please input the heterosturcture interlayer distance (default:3) *")
print("*      !!!attention:oversmall distance may have error in vasp      *")
print("********************************************************************")
interlayer_distance=float(input("* interlayer_distance:"))

slab_thickness=15
print('\n'+"********************************************************************")
print("*    Please input the heterosturcture slab thickness (default:15)   *")
print("*      !!!attention:oversmall thickness may have error in vasp     *")
print("********************************************************************")
slab_thickness=float(input("* slab_thickness:"))

#计算对a,b axis构成平面投影高度
def proj(B,distance):
    s=np.cross(np.array(B[0]),np.array(B[1]))
    S=np.sqrt(s.dot(s))
    V=s.dot(np.array(B[2]))
    proj_distance=distance/(V/S)
    #print(V/S)
    return proj_distance

#对原有的晶胞参数进行取中值，然后缩放
L=[]
L.append((LA[0]+LB[0])/2)
L.append((LA[1]+LB[1])/2)
A[0]=np.multiply(np.array(A[0]),[L[0]/LA[0],L[0]/LA[0],L[0]/LA[0]]).tolist()
A[1]=np.multiply(np.array(A[1]),[L[1]/LA[1],L[1]/LA[1],L[1]/LA[1]]).tolist()
B[0]=np.multiply(np.array(B[0]),[L[0]/LB[0],L[0]/LB[0],L[0]/LB[0]]).tolist()
B[1]=np.multiply(np.array(B[1]),[L[1]/LB[1],L[1]/LB[1],L[1]/LB[1]]).tolist()




#B晶格内最远距离
def longest_z(atom_position):
    maximum=0
    minimum=0
    for i in range(len(atom_position)):
        if float(atom_position[i][2])>maximum:
            maximum=float(atom_position[i][2])
        elif float(atom_position[i][2])<minimum:
            minimum=float(atom_position[i][2])
        else:
            continue
    return maximum-minimum
B_longest=longest_z(atom_position_b)*np.sqrt(np.array(B[2]).dot(np.array(B[2])))
interlayer_distance += B_longest
for i in range(sum_atom_number_a):
    for j in range(3):
        atom_position_a[i][j]=float(atom_position_a[i][j])
    #对笛卡尔坐标
    #atom_position_a[i]=np.multiply(np.array(atom_position_a[i]),[L[0]/LA[0],L[1]/LA[1],1]).tolist()
    
    #a 各原子向上移动interlayer_distance
    atom_position_array=np.array(atom_position_a[i])
    atom_position_array=atom_position_array+proj(B,interlayer_distance)*np.array(B[2])/np.sqrt(np.array(B[2]).dot(np.array(B[2])))
    atom_position_a[i]=atom_position_array.tolist()
        
for i in range(sum_atom_number_b):
   for j in range(3):      
       atom_position_b[i][j]=float(atom_position_b[i][j])
    #对笛卡尔坐标
    #atom_position_b[i]=np.multiply(np.array(atom_position_b[i]),[L[0]/LB[0],L[1]/LB[1],1]).tolist()
      
c_highest=max(max(atom_position_a))*np.sqrt(np.array(B[2]).dot(np.array(B[2])))
c_proj=c_highest+slab_thickness
a=np.array(A[2])*proj(A,c_proj)
A[2]=a.tolist()
b=np.array(B[2])*proj(B,c_proj)
B[2]=b.tolist()



hetero_poscar=[]
#新poscar第一行写入体系名称
name=name_a[0]+"+"+name_b[0]
hetero_poscar.append(name)
#新poscar第二行写入1.0不变
hetero_poscar.append("1.0")
#新poscar第三-五行写入晶胞参数lattice
for i in range(3):
    for j in range(3):
        A[i][j]=str(A[i][j])
    hetero_poscar.append(" ".join(A[i]))
#新poscar第六行写入元素
#新poscar第七行写入元素个数(解决重复元素出现)
element=[]
element_number=[]
c=0
raw=lattice_a[2]+lattice_b[2]
raw_num=lattice_a[3]+lattice_b[3]
for i in range(len(raw)):
    if i in element:
        eindex=element.index(i)
        element_number=str(int(element_number[eindex])+int(raw_num(i)))
        continue
    else:
        element.append(raw[i])
        element_number.append(raw_num[i])
hetero_poscar.append(" ".join(element))
hetero_poscar.append(" ".join(element_number))
#以下行写入原子位置信息
hetero_poscar.append("Direct")
for i in range(len(atom_position_a)):
    for j in range(3):
        atom_position_a[i][j]=str(atom_position_a[i][j])
    hetero_poscar.append(" ".join(atom_position_a[i]))
for i in range(len(atom_position_b)):
    for j in range(3):
        atom_position_b[i][j]=str(atom_position_b[i][j])
    hetero_poscar.append(" ".join(atom_position_b[i]))


f = open('HETERO-POSCAR', 'w+')
for i in hetero_poscar:
    try:
        f.write(i + '\n')
    except:
        print(i)

#计时
time_elapsed = time.time() - since
print('The code run {:.0f}m {:.0f}s'.format(time_elapsed // 60, time_elapsed % 60))
