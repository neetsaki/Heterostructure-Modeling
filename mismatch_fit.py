# -*- coding: utf-8 -*-
"""
@author: Inabahu(Inabahu@tju.edu.cn)
@file: read_POSCAR.py
@time: 2019/7/22
ref:https://github.com/lxf-gzu/vasp_small_script
    
"""
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
if D<1e-10:
    print("****angle test pass****")
    MM=mismatch_matrix(LA[0:2],LB[0:2])
    #由于是二维材料，考虑a,b轴匹配即可,@@@@@@@@@@@@@此处遗留问题，如何断定两轴为a和b轴
    if max(max(MM))<L_threhold:
        print("****lattice test pass****")
        #接合并晶胞函数
    else:
        print("****lattice test fail, run mismatch fit****")
        [m,n]=fit_mismatch(LA[1],LB[1],L_threhold)
        
else:
    print("****angle test fail, need rotation****")
    #夹角不匹配怎么处理啊啊啊啊啊啊啊啊
