import numpy as np
def readdata(filename):
    file1=open(filename,'r')
    data=[]
    
    for row in file1:
        row_str=row.split()
        row_flt=[]
        for i in range(len(row_str)):
            row_flt.append(np.float64(row_str[i]))
        data.append(row_flt)
    
    for j in range(len(data)):
        for k in range(len(data[j])):
            print("%.15f"%data[j][k],end="    ")
        print()
    return data

def readvasp(filename):
    """
    read coordinates in .vasp file
    :param filename: str
        filename
    :return: atom [N][3] float array
             cell_parameter [3][3] float array
        atoms' position and supercell parameter
    """
    file1=open(filename,'r')
    atom_num = []
    atom = []
    n = 0
    cell_parameter = []
    for row in file1:
        n += 1
        row_str = row.split()
        row_flt = []
        if (n<6 and n>2):
            for i in range(len(row_str)):
                row_flt.append(np.float64(row_str[i]))
            cell_parameter.append(row_flt)
        elif (n==7):
            for i in range(len(row_str)):
                atom_num.append(np.float64(row_str[i]))
        elif (n>=9):
            for i in range(len(row_str)):
                row_flt.append(np.float64(row_str[i]))
            if row_flt != []:
                atom.append(row_flt) 
#    print("cell parameter:")
#    for i in range(3):
#        print("    ", end="")
#        for j in range(3):
#            print("%.10f"%cell_parameter[i][j],end="    ")
#        print()
    return atom, cell_parameter

def writevasp(atom, cell_param, material_name, atom_type, atom_num):
    """
    write atom coordintes into .vasp file
    :param atom [N][3] float array
        atom position in fractional coordinates
    :param cell_param [3][3] float array
        cell parameter
    :param material_name str
        name of the material
    :param atom_type [M] str array
        type of all atoms
    :param atom_num [M] int array
        number of each type of atoms
    """
    filename = material_name + '.vasp'
    output = open(filename, mode = 'w')
    print(material_name, file = output)
    print("1.0")
    for i in range(3):
        for j in range(3):
            if j == 1:
                print("%20.10f"%(cell_param[i][j]), file = output, end = "")
            else:
                print("%21.10f"%(cell_param[i][j]), file = output, end = "")
        print("", file = output)
    for j in range(len(atom_type)):
        print("%5s"%(atom_type[j]), file = output, end = "")
    print("", file = output)
    for j in range(len(atom_num)):
        print("%5d"%(atom_num[j]), file = output, end = "")
    print("\n", file = output)
    print("Direct", file = output)
    for n in range(len(atom)):
        for m in range(3):
            if m == 1:
                print("%16.9f"%(atom[n][m]))
            else:
                print("%20.9f"%(atom[n][m]))
                
    return 0
        
        
    
    