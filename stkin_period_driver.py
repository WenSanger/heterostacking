#ABANDONED METHOD
from calc_stacking_area import calc_AA_stacking_area, calc_AB_stacking_area_drttws, calc_AB_stacking_area_rdtws, period_search
import numpy as np
from data_process import readvasp, del_overlap
from twist import twist_cell_mtx_param
from adjust_unitcell import adjust_unitcell, change_param
from plot import plot_atoms
import matplotlib.pyplot as plt
import time
start = time.perf_counter()
output = open('./stacking.txt',mode = 'w')
print("twist_angle(°)          AA_area/cell_area           AB_area/cell_area        cell_area(A^2)",file = output)
output.close()
output = open('./stacking.txt',mode = 'a')
atom_unit1, unit_cell1 = readvasp('WSe2.vasp')
atom_unit2, unit_cell2 = readvasp('MoSe2.vasp')

l1 = 1.0/3.0*np.sqrt(np.power((unit_cell1[1][0]-unit_cell1[0][0]),2)+np.power((unit_cell1[1][1]-unit_cell1[0][1]),2))
l2 = 1.0/3.0*np.sqrt(np.power((unit_cell2[1][0]-unit_cell2[0][0]),2)+np.power((unit_cell2[1][1]-unit_cell2[0][1]),2))

unit_cell1, atom_unit1 = change_param(unit_cell1, atom_unit1)
unit_cell2, atom_unit2 = change_param(unit_cell2, atom_unit2)

unit_cell1, atom_unit1 = adjust_unitcell(unit_cell1, atom_unit1)
unit_cell2, atom_unit2 = adjust_unitcell(unit_cell2, atom_unit2)
#plot_atoms(unit_cell1, atom_unit1)

for param_i in range(1):
    param_i = param_i + 33
    start2 = time.perf_counter()
    #layer1 twist
    rot_cell1_layer1_param, atom_layer1 = twist_cell_mtx_param(unit_cell1, atom_unit1, param_i, 1)
    
    #plot_atoms(rot_cell1_layer1_param, atom_layer1)
    #AB holes
    atom_trans = [1/3,1/3,0]
    atom_unit1_AB1 = []
    for i in range(len(atom_unit1)):
        atom_i = [0.0, 0.0, 0.0]
        for n in range(3):
            atom_i[n] = atom_unit1[i][n] + atom_trans[n]
        atom_unit1_AB1.append(atom_i)
        
    atom_unit1_AB2 = []
    for i in range(len(atom_unit1)):
        atom_i = [0.0, 0.0, 0.0]
        for n in range(3):
            atom_i[n] = atom_unit1[i][n] - atom_trans[n]
        atom_unit1_AB2.append(atom_i)
    
    rot_cell1_layer1_param, atom_layer1_AB1 = twist_cell_mtx_param(unit_cell1, atom_unit1_AB1, param_i, 1)
    rot_cell1_layer1_param, atom_layer1_AB2 = twist_cell_mtx_param(unit_cell1, atom_unit1_AB2, param_i, 1)
    #plot_atoms(rot_cell1_layer1_param, atom_layer1_AB1)
    
    #layer2 twist
    rot_cell2_layer1_param, atom_layer2 = twist_cell_mtx_param(unit_cell2, atom_unit2, param_i, 1)
    rot_cell2_layer2_param, atom_layer2 = twist_cell_mtx_param(unit_cell2, atom_unit2, param_i, 2)
    #cell area
    cell_area = (rot_cell1_layer1_param[0][0]*rot_cell1_layer1_param[1][0]+
    rot_cell1_layer1_param[0][1]*rot_cell1_layer1_param[1][1])*np.sqrt(3)
    
    AB_atom_1_revised = del_overlap(atom_layer1_AB1)
    AB_atom_2_revised = del_overlap(atom_layer1_AB2)
    atom_layer1_revised = del_overlap(atom_layer1)
    atom_layer2_revised = del_overlap(atom_layer2)
    cos_theta = (3*param_i**2+3*param_i+0.5)/(3*param_i**2+3*param_i+1)
    theta = np.arccos(cos_theta)/np.pi * 180
    print("angle = %s°"%theta)
    #
    period_search(theta, atom_layer1_revised, atom_layer2_revised, rot_cell1_layer1_param, rot_cell2_layer1_param, l1, l2)
        
    AA_area, AA_num = calc_AA_stacking_area(theta, atom_layer1_revised, atom_layer2_revised, rot_cell1_layer1_param, rot_cell2_layer1_param, l1, l2)
    AB_area, AB_num = calc_AB_stacking_area_drttws(theta, atom_layer1_revised, atom_layer2_revised, rot_cell1_layer1_param, rot_cell2_layer1_param, l1, l2, AB_atom_1_revised, AB_atom_2_revised)
    
    print('      ',file = output, end = '           ')
    print(theta,file = output, end = '       ')
    print(AA_area/cell_area,file = output, end = '       ')
    print(AB_area/cell_area,file = output, end = '       ')
    print(cell_area,file = output, end = '       ')
    print('', file = output)
    
    end2 = time.perf_counter()
    print('this angle costs time = %s seconds'%(end2-start2))
output.close()
end = time.perf_counter()
print('time = %s seconds'%(end-start))