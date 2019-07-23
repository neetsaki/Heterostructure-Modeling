# -*- coding: utf-8 -*-
"""
@author: Inabahu(Inabahu@tju.edu.cn)
@file: read_POSCAR.py
@time: 2019/7/22

    
"""
#需要实现 异质结结合(√)，异质结旋转
interlayer_distance=3
print("********************************************************************")
print("* Please input the heterosturcture interlayer distance (default:3) *")
print("*      !!!attention:oversmall distance may have error in vasp      *")
print("********************************************************************")
interlayer_distance=float(input("* interlayer_distance:"))

slab_thickness=15
print("********************************************************************")
print("*    Please input the heterosturcture slab thickness (default:3)   *")
print("*      !!!attention:oversmall thickness may have error in vasp     *")
slab_thickness=float(input("* slab_thickness:"))

#计算对a,b axis构成平面投影高度
def proj(B,distance):
    s=np.cross(np.array(B[0]),np.array(B[1]))
    S=np.sqrt(s.dot(s))
    V=s.dot(np.array(B[2]))
    proj_distance=distance/(V/S)
    #print(V/S)
    return proj_distance

for i in range(sum_atom_number_a):
    for j in range(3):
        atom_position_a[i][j]=float(atom_position_a[i][j])
    atom_position_a[i]=np.multiply(np.array(atom_position_a[i]),[L[0]/LA[0],L[1]/LA[1],1]).tolist()
    #a 各原子向上移动interlayer_distance
    atom_position_array=np.array(atom_position_a[i])
    atom_position_array=atom_position_array+proj(B,interlayer_distance)*np.array(B[2])
    atom_position_a[i]=atom_position_array.tolist()
        
for i in range(sum_atom_number_b):
   for j in range(3):      
       atom_position_b[i][j]=float(atom_position_b[i][j])
   atom_position_b[i]=np.multiply(np.array(atom_position_b[i]),[L[0]/LB[0],L[1]/LB[1],1]).tolist()
      
c_highest=max(max(atom_position_a))
c_proj=c_highest+slab_thickness
a=np.array(A[2])*proj(A,c_proj)
A[2]=a.tolist()
b=np.array(B[2])*proj(B,c_proj)
B[2]=b.tolist()

#对原有的晶胞参数进行取中值，然后缩放
L=[]
L.append((LA[0]+LB[0])/2)
L.append((LA[1]+LB[1])/2)
A[0]=np.multiply(np.array(A[0]),[L[0]/LA[0],L[1]/LA[1],1]).tolist()
A[1]=np.multiply(np.array(A[1]),[L[0]/LA[0],L[1]/LA[1],1]).tolist()
B[0]=np.multiply(np.array(B[0]),[L[0]/LB[0],L[1]/LB[1],1]).tolist()
B[1]=np.multiply(np.array(B[1]),[L[0]/LA[0],L[1]/LA[1],1]).tolist()

hetero_poscar=[]
#新poscar第一行写入体系名称
name=name_a[0]+"+"+name_b[0]
hetero_poscar.append(name)
#新poscar第二行写入1.0不变
hetero_poscar.append("0.1")
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





