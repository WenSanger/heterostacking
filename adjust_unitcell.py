import numpy as np
from plot import plot_atoms
def change_param(cell_param, atom):
    """
    change unit cell to get a narrow angle cell parameter
    :param cell param [3][3] float array
        cell parameter of primitive cell
    :param atom [N][3] float array
        atom position in fractional coordinates
    :return cell_param [3][3] float array
            atom_out [N][3] float array
            revised cell parameter and the corresponding fractional atom coordinates 
    """  
    #transform the atom position into cartesian coordinate 
    cell_param_np = np.array(cell_param)
    atom_cartesian = []
    for i in range(3):
        atom_tmp = np.array(atom[i])
        atom_tmp2 = np.matmul(atom_tmp, cell_param_np)
        atom_cartesian.append([atom_tmp2[0],atom_tmp2[1],atom_tmp2[2]])
    #
    cell_param[0][0] = -cell_param[0][0]
    #transform the atom position into fractional coordinate 
    for i in range(3):
        atom_cartesian[i][0] = atom_cartesian[i][0] + cell_param[0][0]
        atom_cartesian[i][1] = atom_cartesian[i][1] + cell_param[0][1]
    cell_param_np = np.array(cell_param)
    atom_out = []
    for i in range(3):
        atom_tmp3 = np.array(atom_cartesian[i])
        atom_tmp4 = np.matmul(atom_tmp3, np.linalg.inv(cell_param_np))
        atom_out.append([atom_tmp4[0],atom_tmp4[1],atom_tmp4[2]])
    #
    cell_param[0][0] = -cell_param[0][0]
    cell_param[1][0] = -cell_param[1][0]
    return cell_param, atom_out
    
def adjust_unitcell(cell_param, atom):
    """
    adjust unit cell to set the last atom at o point and atom1 on the axis x
    :param cell param [3][3] float64 array
        cell parameter of primitive cell
    :param atom [N][3] float64 array
        atom position in fractional coordinates
    :return cell_param_out [3][3] float array
            atom_out_frac [N][3] float array
            revised cell parameter and the corresponding fractional atom coordinates 
    """   
    #transform the atom position into cartesian coordinate 
    cell_param_np = np.array(cell_param)
    atom_cartesian = []
    for i in range(len(atom)):
        atom_tmp = np.array(atom[i])
        atom_tmp2 = np.matmul(atom_tmp, cell_param_np)
        atom_cartesian.append([atom_tmp2[0],atom_tmp2[1],atom_tmp2[2]])
    #translate atoms to set atom0 at o point(only x,y)
    for i in range(len(atom_cartesian)-1):
        for n in range(2):
            atom_cartesian[i][n] = atom_cartesian[i][n]-atom_cartesian[-1][n]
    for n in range(2):
        atom_cartesian[-1][n] = 0.0e0
    #get rotate matrix
    l = np.sqrt(atom_cartesian[0][0]**2+atom_cartesian[0][1]**2)
    cos = atom_cartesian[0][0]/l
    sin = atom_cartesian[0][1]/l
    rotmtx = np.array([[cos, sin , 0],
                       [-sin, cos, 0],
                       [0, 0, 1]])
    #rotate lattice vector to set atom1 on axis x
    cell_param_out = []
    for i in range(3):
        lattice_vectori = np.matmul(rotmtx, cell_param_np[i])
        cell_param_out.append([lattice_vectori[0], lattice_vectori[1], lattice_vectori[2]])
    #transform the atom position back into fractional coordinate
    atom_out_frac = []
    cell_param_np = np.array(cell_param)
    for i in range(3):
        atom_tmp3 = np.array(atom_cartesian[i])
        atom_tmp4 = np.matmul(atom_tmp3, np.linalg.inv(cell_param_np))
        atom_out_frac.append([atom_tmp4[0],atom_tmp4[1],atom_tmp4[2]])
    return cell_param_out, atom_out_frac