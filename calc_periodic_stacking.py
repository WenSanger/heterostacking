from calc_stacking_area import calc_AA_stacking_area, calc_AB_stacking_area_drttws, calc_AB_stacking_area_rdtws, period_search
import numpy as np
from data_process import readvasp, del_overlap, readdata
from twist import twist_cell_matrix_period, mn_rearrange, calc_area
from adjust_unitcell import adjust_unitcell, change_param
from plot import plot_atoms
import matplotlib.pyplot as plt
import time
"""
Calculate AA and AB stacking area in WSe2/WS2 periodical heterostructure.

Before calculation, you need the vasp unit cell file 'WSe2.vasp' and 'WS2.vasp', and
the "result.dat" file containing the parameter:(theta, m1, n1, m2, n2) of periodic structure,
which needed to be calculated by fortran program solve_mn.f90.
"""
start = time.perf_counter()
# read unit cell
atom_unit1, unit_cell1 = readvasp('WS2.vasp')
atom_unit2, unit_cell2 = readvasp('WSe2.vasp')#you can try using WSe2fake.vasp and set m1 = 30, n1 = 0, m2 = 29, n2 = 0
# atom distance
l1 = 1.0/3.0*np.sqrt(np.power((unit_cell1[1][0]-unit_cell1[0][0]),2)+np.power((unit_cell1[1][1]-unit_cell1[0][1]),2))
l2 = 1.0/3.0*np.sqrt(np.power((unit_cell2[1][0]-unit_cell2[0][0]),2)+np.power((unit_cell2[1][1]-unit_cell2[0][1]),2))
# lattice constant
a1 = unit_cell1[0][0]
a2 = unit_cell2[0][0]

print("a1 = ", a1)
print("a2 = ", a2)
print("l1 = ", l1)
print("l2 = ", l2)
# translate atoms
unit_cell1, atom_unit1 = change_param(unit_cell1, atom_unit1)
unit_cell2, atom_unit2 = change_param(unit_cell2, atom_unit2)
# rotate cell parameter to reset unit cell
unit_cell1, atom_unit1 = adjust_unitcell(unit_cell1, atom_unit1)
unit_cell2, atom_unit2 = adjust_unitcell(unit_cell2, atom_unit2)
print(unit_cell1, atom_unit1)
print(unit_cell2, atom_unit2)
# arrange parameter of twist structure according to angle
twist_data = readdata("result.dat")
theta_array, m1_array, n1_array, m2_array, n2_array, point_distance_array = mn_rearrange(twist_data)
#
output = open('./stacking.txt',mode = 'w')
print("twist_angle(°)   AA/period_area           AB/period_area",file = output)
output.close()

# calculate AA/AB stacking in their periodic structure under different angles
#for data_num in range(len(m1_array)):
for data_num in range(1):
    output = open('./stacking.txt',mode = 'a')
    #m1 = m1_array[data_num]
    #n1 = n1_array[data_num]
    #m2 = m2_array[data_num]
    #n2 = n2_array[data_num]
    #theta = theta_array[data_num]
    m1 = 18
    n1 = 11
    m2 = 16
    n2 = 12
    theta = 3.2186618308
    # twist both layers and get AB holes
    rot_layer1_param, atom_layer1, AB_atom_type1, AB_atom_type2 = twist_cell_matrix_period(unit_cell1, atom_unit1, m1, n1, True)
    rot_layer2_param, atom_layer2  = twist_cell_matrix_period(unit_cell2, atom_unit2, m2, n2, False)
    # WRONG METHOD
    #
    print(len(atom_layer1))
    output2 = open('./WSe2.txt',mode = 'w')
    for i in range(len(atom_layer1)):
        for j in range(3):
            print("%20.10f"%(atom_layer1[i][j]),file=output2, end = "")
        print("",file=output2)
    output21 = open('./WSe2_param.txt',mode = 'w')
    for i in range(len(rot_layer1_param)):
        for j in range(3):
            print("%20.10f"%(rot_layer1_param[i][j]),file=output21, end = "")
        print("",file=output21)
    output3 = open('./MoSe2.txt',mode = 'w')
    for i in range(len(atom_layer2)):
        for j in range(3):
            print("%20.10f"%(atom_layer2[i][j]),file=output3, end = "")
        print("",file=output3)
    output31 = open('./MoSe2_param.txt',mode = 'w')
    for i in range(len(rot_layer2_param)):
        for j in range(3):
            print("%20.10f"%(rot_layer2_param[i][j]),file=output31, end = "")
        print("",file=output31)
            
    
    
    # delete overlap
    atom_layer1_revised = del_overlap(atom_layer1)
    atom_layer2_revised = del_overlap(atom_layer2)
    AB_atom_type1_revised = del_overlap(AB_atom_type1)
    AB_atom_type2_revised = del_overlap(AB_atom_type2)
    # calculate AA/AB stacking
    AA_area, AA_num = calc_AA_stacking_area(theta, atom_layer1_revised, atom_layer2_revised, rot_layer1_param, rot_layer2_param, l1, l2)
    AB_area, AB_num = calc_AB_stacking_area_drttws(theta, atom_layer1_revised, atom_layer2_revised, rot_layer1_param, rot_layer2_param, l1, l2, AB_atom_type1_revised, AB_atom_type2_revised)
    period_area = calc_area(rot_layer1_param)
    #
    print('      ',file = output, end = '           ')
    print(theta,file = output, end = '       ')
    print(AA_area/period_area,file = output, end = '       ')
    print(AB_area/period_area,file = output, end = '       ')
    print('', file = output)
    output.close()
end = time.perf_counter()
print('time = %s seconds'%(end-start))