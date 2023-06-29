import numpy as np
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

def del_overlap(atom):
    """
    delete additional overlap atoms whose x,y coordinates are the same
    :param atom: [N][3] float64 array
        atoms' position
    :return: atom_revised
        atoms' position after revision
    """
    if len(atom) == 0:
        print("error in the atoms to be revised")
    atom_revised = []
    atom_x = np.zeros(len(atom))
    atom_y = np.zeros(len(atom))
    for i in range(len(atom)):
        atom_x[i] = atom[i][0]
        atom_y[i] = atom[i][1]
    
    argx = np.argsort(atom_x)
    atom_revised = [atom[argx[0]]]
    i = 0
    while (i<len(atom)-1):
        y_set = [atom_y[argx[i]]]
        i_tmp = i+1
        #if i==0:
        #    print(y_set)
        while (np.abs(atom_x[argx[i]]-atom_x[argx[i_tmp]])<0.001e0):
            #y_set.append(atom_y[argx[i]])
            #if (i_tmp == len(atom)-1):
            #    break
            injudge = 1
            for j in range(len(y_set)):
                #if i==0:
                #    print(np.abs(atom_y[argx[i_tmp]]-y_set[j]))
                if (np.abs(atom_y[argx[i_tmp]]-y_set[j])>=0.001e0):
                    injudge = injudge
                else:
                    injudge = 0
            if injudge == 1:    
                atom_revised.append(atom[argx[i_tmp]])
                y_set.append(atom_y[argx[i_tmp]])
            i_tmp += 1
            if (i_tmp == len(atom)):
                break
        if (i_tmp == len(atom)):
            break
        atom_revised.append(atom[argx[i_tmp]])
        i = i_tmp

    return atom_revised
                
def readdata(filename):
    """
    read data
    :param filename: str
        filename
    :return: data [N][M] float array
    """
    file=open(filename,'r')
    data=[]
    row_num = 0
    for row in file:
        row_num += 1#omit the first line, can be further revised if needed 
        if row_num == 1:
            continue
        row_str=row.split()
        row_float=[]
        for i in range(len(row_str)):
            row_float.append(np.float64(row_str[i]))
        data.append(row_float)
    #for j in range(len(data)):
        #for k in range(len(data[j])):
            #print("%.15f"%data[j][k],end="    ")
        #print()
    return data

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
    print("1.0", file = output)
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
    print("", file = output)
    print("Direct", file = output)
    for n in range(len(atom)):
        for m in range(3):
            if m == 1:
                print("%16.9f"%(atom[n][m]), file = output, end = "")
            else:
                print("%20.9f"%(atom[n][m]), file = output, end = "")
        print("", file = output)
                
    return 0