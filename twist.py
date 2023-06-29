import numpy as np
import matplotlib.pyplot as plt
from plot import plot_atoms
def expand_cell(cell_param, atom, na, nb, nc):
    """
    expand primitive cell into a supercell
    :param cell param [3][3] float array
        cell parameter of primitive cell
    :param atom [N][3] float array
        atom position in fractional coordinates
    :param na, nb, nc integer
        expansion coefficient in direction primitive cell basic vectors' direction a, b, c
    :return [3][3] float array supercell_param
            [N][3] float array supercell_atom
            cell parameter of the expanded supercell and fractional coordinates of supercell atoms
    """
    supercell_atom = []
    supercell_param = [[], [], []]
    for i in range(3):
        supercell_param[0].append(na * cell_param[0][i])
        supercell_param[1].append(nb * cell_param[1][i])
        supercell_param[2].append(nc * cell_param[2][i])
    for n in range(len(atom)):
        for i in range(na):
            for j in range(nb):
                for k in range(nc):
                    atom_tmp = []
                    atom_tmp.append(float(atom[n][0]/na+i*abs(1.0/na)))
                    atom_tmp.append(float(atom[n][1]/nb+j*abs(1.0/nb)))
                    atom_tmp.append(float(atom[n][2]/nc+k*abs(1.0/nc)))
                    supercell_atom.append(atom_tmp)
    return supercell_param, supercell_atom

def find_atoms(atom, cell_param, cell_src_param):
    """
    find atoms which are in the given cell in atom ocean
    :param cell_param [3][3] float array
        original cell parameter
    :param cell_src_param [3][3] float array
        cell parameter of the cell to be searched for
    :param atom [N][3] float array
        atom ocean position in fractional coordinates
    :return atom_src [N'][3] float array
        atom found in the area in cell_src_param under fractional coordinates
    """
    #define search error to add in the boundary atoms
    error1 = 0.000001 / np.linalg.norm(cell_src_param[0])
    error2 = 0.000001 / np.linalg.norm(cell_src_param[1])
    #calculate the inverse matrix of cell_src_param
    src_invs = np.zeros([3,3])
    for i in range(3):
        for j in range(3):
            src_invs[i][j] = cell_src_param[i][j]
    src_invs = np.linalg.inv(src_invs)
    #cell_param MATRIX MUTIPLY src_invs
    M = np.zeros([3,3])
    for i in range(3):
        for j in range(3):
            M[i][j] = cell_param[i][j]
    M = np.matmul(M, src_invs)
    # calculate target atoms
    atom_src = []
    for n in range(len(atom)):
        atom_tmp = np.array([atom[n][0],atom[n][1],atom[n][2]])
        atom_src_tmp = np.matmul(atom_tmp, M)
        if (atom_src_tmp[0]>=0-error1) and (atom_src_tmp[0]<=1+error1) and (atom_src_tmp[1]>=0-error2) and (atom_src_tmp[1]<=1+error2):
            atom_src.append(atom_src_tmp)    
    return atom_src
    

def twist_cell(cell_param, atom, twist_angle, unit_cell):
    """
    twist the given atoms into 2 layers according to given twist_angle, direct method
    :param cell param [3][3] float array
        cell parameter
    :param atom [N][3] float array
        atom position in fractional coordinates
    :param twist_angle float
        twist angle
    :param unit_cell [3][3] float array
        unit cell parameter
    :return atom_layer1, atom_layer2 float64 array [][3]
                twisted atoms
            AB_atom1, AB_atom2
                two types AB atom in layer1 
    """
    # 
    a = 0
    for i in range(3): a += np.power(cell_param[0][i], 2)
    a = np.sqrt(a)
    b = 0
    for i in range(3): b += np.power(cell_param[1][i], 2)
    b = np.sqrt(b)
    c = 0
    for i in range(3): c += np.power(cell_param[1][i] + cell_param[0][i], 2)
    c = np.sqrt(c)
    #
    na = int(np.floor( (c / np.sqrt( c**2 - ((b**2+c**2-a**2)/(2*b))**2 )) + 1 ))
    nb = int(np.floor( (c / np.sqrt( c**2 - ((a**2+c**2-b**2)/(2*a))**2 )) + 1 ))
    # expand cell to creat a basement
    supercell_param, supercell_atom = expand_cell(cell_param, atom, na, nb, 1) 
    atomnum = len(supercell_atom)
    trans = [[0, 0, 0],
             [-1, 0, 0], 
             [0, -1, 0], 
             [-1, -1, 0]]
    # copy the supercell_param to the x,-y and -x,y and -x,-y area
    basement_atom = []
    for j in range(4):
        for i in range(atomnum):
            basement_atom.append([supercell_atom[i][0]+trans[j][0]
, supercell_atom[i][1]+trans[j][1], supercell_atom[i][2]+trans[j][2]]) 
    rot_cell_param1 = [[], [], []]#the cell parameter which differs twist_angle/2 from cell_param
    # twist cell
    #
    rot_mtx = np.array([[np.cos(twist_angle/2), -np.sin(twist_angle/2), 0],
                    [np.sin(twist_angle/2), np.cos(twist_angle/2), 0],
                    [0, 0 ,1]])    
    rot_cell_param1[0] = np.matmul(rot_mtx, cell_param[0])
    rot_cell_param1[1] = np.matmul(rot_mtx, cell_param[1])
    rot_cell_param1[2] = cell_param[2]
    atom_layer1 = find_atoms(basement_atom, supercell_param, rot_cell_param1)#find atoms in area rot_cell_param1
    
    #find atoms in AB position
    basement_atom_type1 = []
    basement_atom_type2 = []
    
    
    trans_1_tot =  np.matmul(rot_mtx, np.array([2/3,1/3,0]))
    trans_2_tot =  np.matmul(rot_mtx, np.array([1/3,2/3,0]))
    
    trans_1 = transform_coord(unit_cell, trans_1_tot, supercell_param)
    trans_2 = transform_coord(unit_cell, trans_2_tot, supercell_param)
    
    for i in range(len(basement_atom)):
        basement_atom_type1.append([basement_atom[i][0]+trans_1[0],basement_atom[i][1]+trans_1[1],basement_atom[i][2]+trans_1[2]])
        basement_atom_type2.append([basement_atom[i][0]+trans_2[0],basement_atom[i][1]+trans_2[1],basement_atom[i][2]+trans_2[2]])

    AB_atom1 = find_atoms(basement_atom_type1, supercell_param, rot_cell_param1)
    AB_atom2 = find_atoms(basement_atom_type2, supercell_param, rot_cell_param1)    
    #
    rot_cell_param2 = [[], [], []]
    rot_mtx = np.array([[np.cos(-twist_angle/2), -np.sin(-twist_angle/2), 0],
                    [np.sin(-twist_angle/2), np.cos(-twist_angle/2), 0],
                    [0, 0 ,1]])    
    rot_cell_param2[0] = np.matmul(rot_mtx, cell_param[0])
    rot_cell_param2[1] = np.matmul(rot_mtx, cell_param[1])
    rot_cell_param2[2] = cell_param[2]
    atom_layer2 = find_atoms(basement_atom, supercell_param, rot_cell_param2)  
    return atom_layer1, atom_layer2, rot_cell_param1, rot_cell_param2, AB_atom1, AB_atom2


def transform_coord(cell_paramA, pos_A, cell_paramB):
    """
    transform fractional atom coordinates under cell_paramA into fractional coordinates under cell_paramB
    :param cell_paramA [3][3] float array
        cell parameter A
    :param pos_A [3] float array
        atom position under cell parameter A
    :param cell parameter B
        cell parameter b
    :return pos_B [3] float array
        atom position under cell parameter B
    """
    a = np.zeros(3)
    b = np.zeros(3)
    c = np.zeros(3)
    for i in range(3):
        a[i] = cell_paramA[0][i]
        b[i] = cell_paramA[1][i]
        c[i] = cell_paramA[2][i]
    pos_cartesian = pos_A[0]*a + pos_A[1]*b + pos_A[2]*c
    M = np.zeros([3,3])
    for i in range(3):
        M[0][i] = cell_paramB[0][i]
        M[1][i] = cell_paramB[1][i]
        M[2][i] = cell_paramB[2][i]
    M = np.linalg.inv(M)    
    tmp = np.matmul(pos_cartesian,M)
    pos_B = []
    for i in range(3):
        pos_B.append(tmp[i])
    return pos_B

def twist_cell_mtx_param(cell_param, atom, param_i, layer):
    """
    ABANDONED METHOD!!!
    create a periodic twist structure according to matrix method which achieves expanding and twisting
    :param cell param [3][3] float64 array
        cell parameter of primitive cell
    :param atom [N][3] float64 array
        atom position in fractional coordinates
    :param param_i int
        parameter controls the twist angle
    :param layer int
        number of the two twisted layer, 1 or 2 
    """  
    cell_param_np = np.array(cell_param)
    if (layer == 1):
        transform_mtx = np.array([[param_i+1,param_i,0],
                              [-param_i,2*param_i+1,0],
                              [0,0,1]])
    elif (layer == 2):
        transform_mtx = np.array([[param_i,param_i+1,0],
                              [-(param_i+1),2*param_i+1,0],
                              [0,0,1]])
    else:
        print("incorrect layer number")
    #calculate the twisted supercell parameter
    rot_cell_param_np = 2*np.matmul(transform_mtx, cell_param_np)
    rot_cell_param = []
    for i in range(3):
        rot_cell_param.append([rot_cell_param_np[i][0], rot_cell_param_np[i][1], rot_cell_param_np[i][2]])
    #evaluate how big the supercell is
    na = int(np.floor(np.linalg.norm(rot_cell_param_np[0] + rot_cell_param_np[1])
    /np.linalg.norm(cell_param[0])))
    nb = int(np.floor(np.linalg.norm(rot_cell_param_np[0] + rot_cell_param_np[1])
    /np.linalg.norm(cell_param[1])))
    #create supercell base
    supercell_param, supercell_atom = expand_cell(cell_param, atom, na, nb, 1)
    atomnum = len(supercell_atom)
    trans = [[0, 0, 0],
             [-1, 0, 0], 
             [0, -1, 0], 
             [-1, -1, 0]]
    # copy the supercell_param to the x,-y and -x,y and -x,-y area
    basement_atom = []
    for j in range(4):
        for i in range(atomnum):
            basement_atom.append([supercell_atom[i][0]+trans[j][0]
, supercell_atom[i][1]+trans[j][1], supercell_atom[i][2]+trans[j][2]]) 
    rot_atom = find_atoms(basement_atom, supercell_param, rot_cell_param_np)
    return rot_cell_param, rot_atom

    
def twist_cell_matrix_period(cell_param, atom, m, n, AB):
    """
    Given a primitive cell and a vector, this function will select atoms in the
    area defined by the input vector and twist them to make the vector on axis x. 
    :param cell_param [3][3] float64 array
        cell parameter of primitive cell
    :param atom [N][3] float64 array
        atom position in fractional coordinates
    :param m,n int
        parameter defines the input vector, m*cell_param[0]+n*cell_param[1]
    :param AB logical
        parameter determines whether to calculate AB holes
    """  
    cell_param_np = np.array(cell_param)
    #evaluate how big the supercell is
    bound_vec = [[],
                 [],
                 [0,0,1]]
    for i in range(3):
        bound_vec[0].append(m*cell_param[0][i]+n*cell_param[1][i])
        bound_vec[1].append((-n)*cell_param[0][i]+(n+m)*cell_param[1][i])
    bound_vec_np = np.array(bound_vec)
    na = int(np.floor(np.linalg.norm(bound_vec_np[0] + bound_vec_np[1])
    /np.linalg.norm(cell_param[0])*2/np.sqrt(3)))+1
    nb = int(np.floor(np.linalg.norm(bound_vec_np[0] + bound_vec_np[1])
    /np.linalg.norm(cell_param[1])*2/np.sqrt(3)))+1
    #create supercell base
    supercell_param, supercell_atom = expand_cell(cell_param, atom, na, nb, 1)
    atomnum = len(supercell_atom)
    trans = [[0, 0, 0],
             [-1, 0, 0], 
             [0, -1, 0], 
             [-1, -1, 0]]
    # translate the supercell_param to the x,-y and -x,y and -x,-y area
    basement_atom = []
    for j in range(4):
        for i in range(atomnum):
            basement_atom.append([supercell_atom[i][0]+trans[j][0]
, supercell_atom[i][1]+trans[j][1], supercell_atom[i][2]+trans[j][2]]) 
    
    plot_atoms(supercell_param, basement_atom)
    plt.plot([0, bound_vec[0][0],bound_vec[0][0]+bound_vec[1][0]],[0, bound_vec[0][1],bound_vec[0][1]+bound_vec[1][1]], c = 'purple')
    plt.plot([0, bound_vec[1][0],bound_vec[0][0]+bound_vec[1][0]],[0, bound_vec[1][1],bound_vec[0][1]+bound_vec[1][1]], c = 'purple')
    plt.show()
    twist_atom = find_atoms(basement_atom, supercell_param, bound_vec_np)
    # AB holes
    #find atoms in AB position
    basement_atom_type1 = []
    basement_atom_type2 = []
    trans_type1_prm =  [1/3, 1/3, 0]
    trans_type2_prm =  [-1/3, -1/3, 0]
    trans_1 = transform_coord(cell_param, trans_type1_prm, supercell_param)
    trans_2 = transform_coord(cell_param, trans_type2_prm, supercell_param)
    for i in range(len(basement_atom)):
        basement_atom_type1.append([basement_atom[i][0]+trans_1[0],basement_atom[i][1]+trans_1[1],basement_atom[i][2]+trans_1[2]])
        basement_atom_type2.append([basement_atom[i][0]+trans_2[0],basement_atom[i][1]+trans_2[1],basement_atom[i][2]+trans_2[2]])
    AB_atom1 = find_atoms(basement_atom_type1, supercell_param, bound_vec_np)
    AB_atom2 = find_atoms(basement_atom_type2, supercell_param, bound_vec_np)    
    twist_cell_param = [[],[],[]]
    a = np.linalg.norm(bound_vec_np[0])/np.linalg.norm(cell_param[0])

    for i in range(3):
        for j in range(3):
            twist_cell_param[i].append(cell_param[i][j]*a)
    #logical judge, whether to return AB holes
    if (AB == True):
        return twist_cell_param, twist_atom, AB_atom1, AB_atom2
    else:
        return twist_cell_param, twist_atom

def calc_area(cell_parameter):
    """
    rearrange structure parameter according to the rising sequence of angle
    :param twist_data [N][6] float [theta, m1, n1, m2, n2, point_distance]
    return theta_array, point_distance_array [N] float array 
           m1_array, n1_array, m2_array, n2_array [N] integer array 
    """
    cell_parameter_np = np.array(cell_parameter)
    area = np.linalg.norm(np.cross(cell_parameter_np[0],cell_parameter_np[1]))
    return area
def mn_rearrange(twist_data):
    """
    rearrange structure parameter according to the rising sequence of angle
    :param twist_data [N][6] float [theta, m1, n1, m2, n2, point_distance]
    return theta_array, point_distance_array [N] float array 
           m1_array, n1_array, m2_array, n2_array [N] integer array 
    """
    theta_array_tmp = []
    theta_array = []
    m1_array = []
    n1_array = []
    m2_array = []
    n2_array = []
    point_distance_array = []
    
    for i in range(len(twist_data)):
        theta_array_tmp.append(twist_data[i][0])
    # arrange the sequence of twist angle 
    theta_pointer = np.argsort(theta_array_tmp)
    for i in range(len(twist_data)):
        if (i>=1)and(np.abs(theta_array_tmp[theta_pointer[i]] - theta_array_tmp[theta_pointer[i-1]]) < 0.000001):
            #delete repeated structure
            if (int(twist_data[i][1]) > int(twist_data[i-1][1])) \
            and (int(twist_data[i][2]) > int(twist_data[i-1][2])) \
            and (int(twist_data[i][3]) > int(twist_data[i-1][3])) \
            and (int(twist_data[i][4]) > int(twist_data[i-1][4])):
                continue
            elif (int(twist_data[i-1][1]) > int(twist_data[i][1])) \
            and (int(twist_data[i-1][2]) > int(twist_data[i][2])) \
            and (int(twist_data[i-1][3]) > int(twist_data[i][3])) \
            and (int(twist_data[i-1][4]) > int(twist_data[i][4])):
                theta_array[i-1] = twist_data[theta_pointer[i]][0]
                m1_array[i-1] = twist_data[theta_pointer[i]][1]
                n1_array[i-1] = twist_data[theta_pointer[i]][2]
                m2_array[i-1] = twist_data[theta_pointer[i]][3]
                n2_array[i-1] = twist_data[theta_pointer[i]][4]
                continue
        #if (twist_data[theta_pointer[i]][5]<0.0001):
        theta_array.append(twist_data[theta_pointer[i]][0])
        m1_array.append(int(twist_data[theta_pointer[i]][1]))
        n1_array.append(int(twist_data[theta_pointer[i]][2]))
        m2_array.append(int(twist_data[theta_pointer[i]][3]))
        n2_array.append(int(twist_data[theta_pointer[i]][4]))
        point_distance_array.append(twist_data[theta_pointer[i]][5])
        
    return theta_array, m1_array, n1_array, m2_array, n2_array, point_distance_array
    
    
    
    
    