from calc_stacking_area import calc_AA_stacking_area, calc_AB_stacking_area_drttws, calc_AB_stacking_area_rdtws
import numpy as np
from data_process import readvasp, del_overlap
from twist import expand_cell, twist_cell
import time
"""
Calculate AA and AB stacking area in WSe2/MoSe2 heterostructure, also can be used to calculate other twisted hexagon structure.
if read == 0, program will form twisted structure according to a series of different angles,
else, read given twisted structure and directly calculate stacking.

Before calculation, you need the vasp unit cell file 'WSe2.vasp' and 'MoSe2.vasp', if read == 1,
vasp files of twised atoms in two layers are needed. 
"""
start = time.perf_counter()
atom_unit1, unit_cell1 = readvasp('WSe2.vasp')
atom_unit2, unit_cell2 = readvasp('MoSe2.vasp')
#atom distance
l1 = 1.0/3.0*np.sqrt(np.power((unit_cell1[1][0]-unit_cell1[0][0]),2)+np.power((unit_cell1[1][1]-unit_cell1[0][1]),2))
l2 = 1.0/3.0*np.sqrt(np.power((unit_cell2[1][0]-unit_cell2[0][0]),2)+np.power((unit_cell2[1][1]-unit_cell2[0][1]),2))
print("l1=%s,l2=%s"%(l1,l2))
read = 1#
output = open('./stacking.txt',mode = 'w')
print("twist_angle(°)   AA_num   AA_area(A^2)           AB_num   AB_area(A^2)",file = output)
output.close()
output = open('./stacking.txt',mode = 'a')
N = 300
angle = 0.0e0# if read == 1, please change it to the twist angle of given twisted structure, if read == 0, angle = 0.0e0
for angle in range(N):
    angle += 60.0/N
#for angle in range(1):
    #angle = 0.06e0#°
    print(("angle = %s°")%angle)
    if (read == 0):    
        #expand cell
        na = 40
        nb = 40
        nc = 1
        cell_parameter1, atom1 = expand_cell(unit_cell1, atom_unit1, na, nb, 1)
        cell_parameter2, atom2 = expand_cell(unit_cell2, atom_unit2, na, nb, 1)
        twist_angle = angle/180*np.pi
        # get twisted atoms in two layers and AB holes
        atom_layer1, atom_layer2_discard,  rot_cell1_layer1_param, rot_cell1_layer2_param, AB_atom_1, AB_atom_2 = twist_cell(cell_parameter1, atom1, twist_angle, unit_cell1)
        atom_layer1_discard, atom_layer2,  rot_cell2_layer1_param1, rot_cell2_layer2_param2, AB_atom_1_discard, AB_atom_2_discard = twist_cell(cell_parameter2, atom2, twist_angle, unit_cell2)
        # delete overlap atoms
        AB_atom_1_revised = del_overlap(AB_atom_1)
        AB_atom_2_revised = del_overlap(AB_atom_2)
        del atom_layer1_discard, atom_layer2_discard, AB_atom_1_discard, AB_atom_2_discard
    else:    
        atom_layer1, cell_parameter1 = readvasp('WSe21.vasp')
        atom_layer2, cell_parameter2 = readvasp('MoSe22.vasp')
    # delete additional atoms whose x,y coordinates are the same
    atom_layer1_revised = del_overlap(atom_layer1)
    atom_layer2_revised = del_overlap(atom_layer2)
    # calculate AA stacking area
    AA_area, AA_num = calc_AA_stacking_area(angle, atom_layer1_revised, atom_layer2_revised, cell_parameter1, cell_parameter2, l1, l2)
    # calculate AB stacking area
    if (read == 0):
        AB_area, AB_num = calc_AB_stacking_area_drttws(angle, atom_layer1_revised, atom_layer2_revised, cell_parameter1, cell_parameter2, l1, l2, AB_atom_1_revised, AB_atom_2_revised)
    else:
        AB_area, AB_num = calc_AB_stacking_area_rdtws(atom_layer1_revised, atom_layer2_revised, cell_parameter1, cell_parameter2, l1, l2)
        break
    #
    print('      ',file = output, end = '           ')
    print(angle,file = output, end = '       ')
    print(AA_num,file = output, end = '       ')
    print(AA_area,file = output, end = '       ')
    print(AB_num,file = output, end = '       ')
    print(AB_area,file = output, end = '       ')
    print('', file = output)
output.close()
end = time.perf_counter()
print('time = %s seconds'%(end-start))