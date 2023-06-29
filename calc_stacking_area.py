import numpy as np
import matplotlib.pyplot as plt
from data_process import readvasp, del_overlap
from twist import expand_cell, twist_cell, transform_coord
def search_AA_hexagon(layer1_x, layer1_y, layer2_x, layer2_y, a, l):
    """
    Calculate AA stacking hexagon

    :param layer1_x, layer1_y: float64 array
        the x and y index of atoms on the comparing layer
    :param layer2_x, layer2_y: float64 array
        the x and y index of atoms on the layer to be selected
    :param a: float
        parameter defining the searching length
    :param l: float
        the side length of hexagon
    :
    :return: float
        twisting angle in RADIANs, NOT degrees
    """ 
    atom_num_1 = len(layer1_x)
    atom_num_2 = len(layer2_x)
    srch_l = l*1.05
    # step 1: find atoms in layer2 that are close to atoms in layer 1
    #arrange the x and y coordinate of atoms in each layer
    layer1_argx = np.argsort(layer1_x)
    #layer1_argy = np.argsort(layer1_y)
    layer2_argx = np.argsort(layer2_x)    
    #layer2_argy = np.argsort(layer2_y)
    # select the atoms in layer2 which are close to each atom in layer1
    atom_near = []
    nlable = 0 # label used towards searching for layer2 atom
    for i in range(len(layer1_argx)):
        atom_near_i = []
        # search from x-a to x+a
        while((layer2_x[layer2_argx[nlable]]>=(layer1_x[layer1_argx[i]]-a))):
            if (nlable != 0):
                nlable -= 1
            else:
                break
        while((layer2_x[layer2_argx[nlable]]<(layer1_x[layer1_argx[i]]-a))):
            nlable += 1
            if (nlable >= atom_num_2):
                nlable -= 1
                break
    #now we start searching from x-a
        while(layer2_x[layer2_argx[nlable]]<(layer1_x[layer1_argx[i]]+a)):
            r12 = np.sqrt(np.power((layer2_x[layer2_argx[nlable]]-layer1_x[layer1_argx[i]]),2)
    +np.power((layer2_y[layer2_argx[nlable]]-layer1_y[layer1_argx[i]]),2))
            if (r12 <= a):
                #print(layer2_argx[nlable])
                atom_near_i.append(layer2_argx[nlable])
                break
            nlable += 1
            if (nlable >= atom_num_2):
                nlable -= 1
                break
        atom_near.append(atom_near_i)
    #
    layer2_near_x=[]
    layer2_near_y=[]
    for i in range(len(atom_near)):
        if (atom_near[i]):
#        and((layer2_x[atom_near[i][0]] not in layer2_near_x)
#    or(layer2_y[atom_near[i][0]] not in layer2_near_y))):
            #for j in range(len(atom_near[i])):
            layer2_near_x.append(layer2_x[atom_near[i][0]])
            layer2_near_y.append(layer2_y[atom_near[i][0]])

    # step 2: find which atoms in the selected group can be linked into hexagon
    #rearrange layer2_near_x
    layer2_near_argx = np.argsort(layer2_near_x)
    layer2_near_x_tmp = np.zeros(len(layer2_near_x))
    layer2_near_y_tmp = np.zeros(len(layer2_near_x))
    
    for i in range(len(layer2_near_x)):
        layer2_near_x_tmp[i] = layer2_near_x[i]
        layer2_near_y_tmp[i] = layer2_near_y[i]
    for i in range(len(layer2_near_x)):
        layer2_near_x[i] = layer2_near_x_tmp[layer2_near_argx[i]]
        layer2_near_y[i] = layer2_near_y_tmp[layer2_near_argx[i]]


    #fig10 = plt.figure(figsize=(16,10))
    #plt.scatter(layer1_x, layer1_y,s = 15, c='r')
    #plt.scatter(layer2_x, layer2_y,s = 15)
    #plt.scatter(layer2_near_x, layer2_near_y, s = 15, c = 'black')  
    #plt.axis("equal")     
    #plt.show()  

    #print(np.sqrt(np.power((layer2_near_x[0]-layer2_near_x[1]),2)+
    #np.power((layer2_near_y[0]-layer2_near_y[1]),2)))
    #delete overlap

    
    
    del layer2_near_argx
    del layer2_near_x_tmp
    del layer2_near_y_tmp
    
    tripoints_pair=[]
    nlable = 0
    n_near = len(layer2_near_x)
    for i in range(n_near):
        nearest_pair_i=[]
        nearest_pair_i.append(i)    
         # search from x-srch_l to x+srch_l    
        while((layer2_near_x[nlable]>=(layer2_near_x[i]-srch_l))):
            if (nlable != 0):
                nlable -= 1
            else:
                break
        while((layer2_near_x[nlable]<(layer2_near_x[i]-srch_l))):
            nlable += 1    
        #now we start searching from x-srch_l
        while(layer2_near_x[nlable]<(layer2_near_x[i]+srch_l)):
            if ((np.sqrt(np.power((layer2_near_x[nlable]-layer2_near_x[i]),2)+
    np.power((layer2_near_y[nlable]-layer2_near_y[i]),2)) < l*1.2) and (i!=nlable)):
                nearest_pair_i.append(nlable)
            nlable += 1
            if (nlable>=n_near):
                nlable -= 1
                break
        #print(nearest_pair_i)
        if (len(nearest_pair_i)==3):
            #tripoints_pair_i.sort()
            tripoints_pair.append(nearest_pair_i)
        elif (len(nearest_pair_i)==4):
            #tripoints_pair_i.sort()
            tripoints_pair.append([nearest_pair_i[0],nearest_pair_i[1],nearest_pair_i[2]])
            tripoints_pair.append([nearest_pair_i[0],nearest_pair_i[1],nearest_pair_i[3]])
            tripoints_pair.append([nearest_pair_i[0],nearest_pair_i[2],nearest_pair_i[3]])
    #tripoints_pair:three linking points with the mid one being the first in the array
    #print(tripoints_pair)
    #
    hexagon_atom = []
    n_tri = len(tripoints_pair)
    i = 0
    #print(tripoints_pair)
    while(i<n_tri):
        hexagon_exist = False
        tripoints_pair_tmp = []
        hexagon_pair = []
        hexagon_pair.append(tripoints_pair[i])
        #find the neighbor linking tripoints_pair
        for j in range(n_tri):
            if (i!=j):
                
                if ((hexagon_pair[0][1] in tripoints_pair[j])
    and(not(hexagon_pair[0][0] in tripoints_pair[j]))):
                    for k in range(2):
                        if (tripoints_pair[j][k+1]!=tripoints_pair[i][1]):
                            r12=np.sqrt(np.power((layer2_near_x[tripoints_pair[j][k+1]]-layer2_near_x[tripoints_pair[i][2]]),2)+
    np.power((layer2_near_y[tripoints_pair[j][k+1]]-layer2_near_y[tripoints_pair[i][2]]),2))
                            if (r12<2*l):
                                hexagon_pair.append(tripoints_pair[j])
        #find the neighbor's neighbor linking tripoints_pair
        if (len(hexagon_pair) == 2):
            for p in range(n_tri):#len(tripoints_pair)):
                #print(tripoints_pair[p])
                #print(p)
                if ((hexagon_pair[0][2] in tripoints_pair[p])and(hexagon_pair[0][0] not in tripoints_pair[p])):
                    for m in range(2):
                        if (tripoints_pair[p][m+1] in hexagon_pair[1]):
                            # now, hexagon exists
                            hexagon_exist = True
                            hexagon_pair.append(tripoints_pair[p])
                            #find the other three tripoints_pair
                            for q in range(n_tri):
                                if ((hexagon_pair[0][0] in tripoints_pair[q])
    and(hexagon_pair[1][0] in tripoints_pair[q])):
                                    hexagon_pair.append(tripoints_pair[q])
                                if ((hexagon_pair[0][0] in tripoints_pair[q])
    and(hexagon_pair[2][0] in tripoints_pair[q])):
                                    hexagon_pair.append(tripoints_pair[q])
                                if ((hexagon_pair[1][0] in tripoints_pair[q])
    and(hexagon_pair[2][0] in tripoints_pair[q])):
                                    hexagon_pair.append(tripoints_pair[q])
        if (hexagon_exist):
            hexagon_atom.append(hexagon_pair)
            #delet hexagon_pair from tripoints_pair
            for n_del in range(n_tri):#len(tripoints_pair)):
                if (tripoints_pair[n_del] not in hexagon_pair):
                    tripoints_pair_tmp.append(tripoints_pair[n_del])
            tripoints_pair = []
            for n_copy in range(len(tripoints_pair_tmp)):
                tripoints_pair.append(tripoints_pair_tmp[n_copy])
            n_tri = len(tripoints_pair)
        if not (hexagon_exist):
            i += 1
    #print(hexagon_atom) 
    #print(len(hexagon_atom))
    return layer2_near_x, layer2_near_y, hexagon_atom

def calc_AA_stacking_area(angle, atom_layer1, atom_layer2, cell_parameter1, cell_parameter2, l1, l2):
    """
    Calculate AA stacking hexagon

    :param atom_layer1, atom_layer2: [][] float64 array
        coordinates of the layer1 and layer2 atoms(reduce same x,y but different z atoms)
    :param cell_parameter1, cell_parameter2: [3][3] float64 array
        cell_parameter of the supercell layer1 and layer2
    :param l1, l2: float
        the side length of hexagon in layer1 and layer2
    :
    :return: AA_area
        area of AA stacking
    """
    a = 0.6#searching length
    atom_num_1 = len(atom_layer1)
    atom_num_2 = len(atom_layer2)
    atom_layer1 = np.array(atom_layer1,dtype=np.float64)
    atom_layer2 = np.array(atom_layer2,dtype=np.float64)
    cell_parameter1 = np.array(cell_parameter1)
    cell_parameter2 = np.array(cell_parameter2)
    
    # direct to Cartesian
    for i in range(atom_num_1):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer1[i],cell_parameter1)
        atom_layer1[i] = atom_tmp
    for j in range(atom_num_2):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer2[j],cell_parameter2)
        atom_layer2[j] = atom_tmp
    #
    layer1_x = np.zeros(atom_num_1)
    layer1_y = np.zeros(atom_num_1)
    layer2_x = np.zeros(atom_num_2)
    layer2_y = np.zeros(atom_num_2)
    for i in range(atom_num_1):
        layer1_x[i] = atom_layer1[i][0]
        layer1_y[i] = atom_layer1[i][1]
    for i in range(atom_num_2):
        layer2_x[i] = atom_layer2[i][0]
        layer2_y[i] = atom_layer2[i][1]   

    fig1 = plt.figure(figsize=(16,10))
    plt.scatter(layer1_x, layer1_y,s = 15, c='r')
    plt.scatter(layer2_x, layer2_y,s = 15, c='black')
    #plt.xlim(-120,180)
    #plt.ylim(-60,180)
    picname = str(angle) + '° stacking'
    plt.title(picname)
    #plt.tight_layout()  
    plt.axis("equal") 
    picfile = './pic/' + str(angle) + '°.png'
    plt.savefig(picfile, bbox_inches='tight')
    plt.show()
    #
    layer2_near_x, layer2_near_y, hexagon_atom = search_AA_hexagon(layer1_x, layer1_y, layer2_x, layer2_y, a, l2)
    #
    AA_area = 3/2*np.sqrt(3)*l1*l1*len(hexagon_atom)
    AA_atom = []
    for i in range(len(hexagon_atom)):
        for j in range(3):
            AA_atom.append(hexagon_atom[i][j])
    AA_atom = np.unique(AA_atom)
    AA_num = len(hexagon_atom)
    print(AA_num)
    del hexagon_atom
    layer2_hexagon_x = []
    layer2_hexagon_y = []
    for i in range(len(AA_atom)):
        layer2_hexagon_x.append(layer2_near_x[AA_atom[i]])
        layer2_hexagon_y.append(layer2_near_y[AA_atom[i]])
    
    
    #
    AA = []
    for i in range(len(layer2_hexagon_x)):
        AA.append([layer2_hexagon_x[i],layer2_hexagon_y[i],0])
    cartesian = [[1.0e0,0,0],
                 [0,1.0e0,0],
                 [0,0,1.0e0]]
    AA_frac=[]
    for i in range(len(AA)):
        AA_tmp = transform_coord(cartesian,AA[i],cell_parameter2)
        AA_frac.append(AA_tmp)
    #print(AA_frac)
    output_AA = open("./AA.txt",mode='w')
    for i in range(len(AA_frac)):
        for j in range(3):
            print("%20.10f"%(AA_frac[i][j]),file=output_AA, end = "")
        print("",file=output_AA)
        
    
    
    
    
    # plt after final selection
    fig2 = plt.figure(figsize=(16,10))
    plt.scatter(layer1_x, layer1_y,s = 5, c='pink')
    plt.scatter(layer2_x, layer2_y,s = 5, c='orange')
    plt.scatter(layer2_hexagon_x, layer2_hexagon_y,s = 5, c='r')
    plt.axis("equal") 
    #plt.xlim(-120,180)
    #plt.ylim(-60,180)
    
    #plt.tight_layout()
    #picname = str(angle) + '° AA stacking'
    #plt.title(picname)
    #picfile = './pic/AA' + str(angle) + '°.png'
    #plt.savefig(picfile, bbox_inches='tight')
    #plt.show()
    print("the area of AA stacking is %.5f A^2"%AA_area)            
    return AA_area, AA_num

def calc_AB_stacking_area_rdtws(atom_layer1, atom_layer2, cell_parameter1, cell_parameter2, l1, l2):
    """
    Calculate AB stacking hexagon

    :param atom_layer1, atom_layer2: [][] float64 array
        coordinates of the layer1 and layer2 atoms(reduce same x,y but different z atoms)
    :param cell_parameter1, cell_parameter2: [3][3] float64 array
        cell_parameter of supercell in layer1 and layer2
    :param l1, l2: float
        the side length of hexagon in layer1 and layer2
    :
    :return: AB_area
        area of AB stacking
    """
    a = 0.6#searching length
    error = 0.01
    #六边形边长
    atom_num_1 = len(atom_layer1)
    atom_num_2 = len(atom_layer2)
    atom_layer1 = np.array(atom_layer1,dtype=np.float64)
    atom_layer2 = np.array(atom_layer2,dtype=np.float64)
    cell_parameter1 = np.array(cell_parameter1)
    cell_parameter2 = np.array(cell_parameter2)
    # direct to Cartesian
    for i in range(atom_num_1):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer1[i],cell_parameter1)
        atom_layer1[i] = atom_tmp
    for j in range(atom_num_2):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer2[j],cell_parameter2)
        atom_layer2[j] = atom_tmp
    #
    #find the edge of MX2
    #search for corner point
    left_point = 0
    right_point = 0
    for i in range(atom_num_1-1):
        if (atom_layer1[i+1][0]< atom_layer1[left_point][0]):
            left_point = i+1
        if (atom_layer1[i+1][0]>atom_layer1[right_point][0]):
            right_point = i+1
    c = (atom_layer1[left_point][1]-atom_layer1[right_point][1])/(atom_layer1[left_point][0]-atom_layer1[right_point][0])
    up_point = 0
    down_point = 0
    for i in range(atom_num_1-1):
        if((atom_layer1[i+1][1]-c*atom_layer1[i+1][0])>(atom_layer1[up_point][1]-c*atom_layer1[up_point][0])):
            up_point = i+1
        if((atom_layer1[i+1][1]-c*atom_layer1[i+1][0])<(atom_layer1[down_point][1]-c*atom_layer1[down_point][0])):
            down_point = i+1
    
    #search for edge point(without corner point)
    left_edge = []
    right_edge = []
    down_edge = []
    up_edge = []
    
    c1 = (atom_layer1[left_point][1]-atom_layer1[down_point][1])/(atom_layer1[left_point][0]-atom_layer1[down_point][0])
    d1 = atom_layer1[left_point][1] - c1 *  atom_layer1[left_point][0]
    c2 = (atom_layer1[right_point][1]-atom_layer1[down_point][1])/(atom_layer1[right_point][0]-atom_layer1[down_point][0])
    d2 = atom_layer1[right_point][1] - c2 * atom_layer1[right_point][0]
    c3 = (atom_layer1[right_point][1]-atom_layer1[up_point][1])/(atom_layer1[right_point][0]-atom_layer1[up_point][0])
    d3 = atom_layer1[right_point][1] - c3 * atom_layer1[right_point][0]
    c4 = (atom_layer1[left_point][1]-atom_layer1[up_point][1])/(atom_layer1[left_point][0]-atom_layer1[up_point][0])
    d4 = atom_layer1[left_point][1] - c4 * atom_layer1[left_point][0]
    
    corner = [left_point,right_point,up_point,down_point]
    for i in range(atom_num_1):
        if ((np.abs(atom_layer1[i][1]-c1 * atom_layer1[i][0] - d1)<error)and(i not in corner)):
            left_edge.append(i)
        if ((np.abs(atom_layer1[i][1]-c2 * atom_layer1[i][0] - d2)<error)and(i not in corner)):
            down_edge.append(i)
        if ((np.abs(atom_layer1[i][1]-c3 * atom_layer1[i][0] - d3)<error)and(i not in corner)):
            right_edge.append(i)
        if ((np.abs(atom_layer1[i][1]-c4 * atom_layer1[i][0] - d4)<error)and(i not in corner)):
            up_edge.append(i)

    # construct total AB stacking(type1)
    #
    layerAB_x = []
    layerAB_y = []
    for i in range(atom_num_1):
        if (i not in up_edge)and(i not in right_edge)and(i not in corner):
            layerAB_x.append(atom_layer1[i][0] + l1/np.sqrt(1+np.power(c2,2)) - l1/np.sqrt(1+np.power(c1,2)))
            layerAB_y.append(atom_layer1[i][1] + np.sign(c2)*l1/np.sqrt(1+np.power(1.0/c2,2)) - np.sign(c1)*l1/np.sqrt(1+np.power(1.0/c1,2)))
    layerAB_x.append(atom_layer1[down_point][0] + l1/np.sqrt(1+np.power(c2,2)) - l1/np.sqrt(1+np.power(c1,2)))
    layerAB_y.append(atom_layer1[down_point][1] + np.sign(c2)*l1/np.sqrt(1+np.power(1.0/c2,2)) - np.sign(c1)*l1/np.sqrt(1+np.power(1.0/c1,2)))
    #
    
    
    
    for i in range(len(left_edge)):
        layerAB_x.append(atom_layer1[left_edge[i]][0] + l1/np.sqrt(1+np.power(c1,2)) )
        layerAB_y.append(atom_layer1[left_edge[i]][1] + np.sign(c1)*l1/np.sqrt(1+np.power(1.0/c1,2)) )
    layerAB_x.append(atom_layer1[left_point][0] + l1/np.sqrt(1+np.power(c1,2)) )
    layerAB_y.append(atom_layer1[left_point][1] + np.sign(c1)*l1/np.sqrt(1+np.power(1.0/c1,2)) )
        
    for i in range(len(down_edge)):
        layerAB_x.append(atom_layer1[down_edge[i]][0] + 2*l1/np.sqrt(1+np.power(c2,2)) )
        layerAB_y.append(atom_layer1[down_edge[i]][1] + np.sign(c2)*2*l1/np.sqrt(1+np.power(1.0/c2,2)) )
    layerAB_x.append(atom_layer1[down_point][0] + 2*l1/np.sqrt(1+np.power(c2,2)) )
    layerAB_y.append(atom_layer1[down_point][1] + np.sign(c2)*2*l1/np.sqrt(1+np.power(1.0/c2,2)) )
    # construct total AB stacking(type2)
    #
    layerAB2_x = []
    layerAB2_y = []
    for i in range(atom_num_1):
        if (i not in down_edge)and(i not in left_edge)and(i not in corner):
            layerAB2_x.append(atom_layer1[i][0] - l1/np.sqrt(1+np.power(c2,2)) + l1/np.sqrt(1+np.power(c1,2)))
            layerAB2_y.append(atom_layer1[i][1] - np.sign(c2)*l1/np.sqrt(1+np.power(1.0/c2,2)) + np.sign(c1)*l1/np.sqrt(1+np.power(1.0/c1,2)))
    
    layerAB2_x.append(atom_layer1[up_point][0] - l1/np.sqrt(1+np.power(c2,2)) + l1/np.sqrt(1+np.power(c1,2)))
    layerAB2_y.append(atom_layer1[up_point][1] - np.sign(c2)*l1/np.sqrt(1+np.power(1.0/c2,2)) + np.sign(c1)*l1/np.sqrt(1+np.power(1.0/c1,2)))
    
    #
    for i in range(len(left_edge)):
        layerAB2_x.append(atom_layer1[right_edge[i]][0] - l1/np.sqrt(1+np.power(c3,2)) )
        layerAB2_y.append(atom_layer1[right_edge[i]][1] - np.sign(c3)*l1/np.sqrt(1+np.power(1.0/c3,2)) )
    
    layerAB2_x.append(atom_layer1[right_point][0] - l1/np.sqrt(1+np.power(c3,2)) )
    layerAB2_y.append(atom_layer1[right_point][1] - np.sign(c1)*l1/np.sqrt(1+np.power(1.0/c3,2)) )
     

    
    for i in range(len(down_edge)):
        layerAB2_x.append(atom_layer1[up_edge[i]][0] - 2*l1/np.sqrt(1+np.power(c4,2)) )
        layerAB2_y.append(atom_layer1[up_edge[i]][1] - np.sign(c4)*2*l1/np.sqrt(1+np.power(1.0/c4,2)) )
    
    layerAB2_x.append(atom_layer1[up_point][0] - 2*l1/np.sqrt(1+np.power(c4,2)) )
    layerAB2_y.append(atom_layer1[up_point][1] - np.sign(c4)*2*l1/np.sqrt(1+np.power(1.0/c4,2)) )
    
    # step 1: find atoms in layer2 that are close to atoms in layer 1
    layer1_x = np.zeros(atom_num_1)
    layer1_y = np.zeros(atom_num_1)
    layer2_x = np.zeros(atom_num_2)
    layer2_y = np.zeros(atom_num_2)
    for i in range(atom_num_1):
        layer1_x[i] = atom_layer1[i][0]
        layer1_y[i] = atom_layer1[i][1]
    for i in range(atom_num_2):
        layer2_x[i] = atom_layer2[i][0]
        layer2_y[i] = atom_layer2[i][1]       
    
    del atom_layer1
    del atom_layer2
    layer2_near_x, layer2_near_y, hexagon_atom = search_AA_hexagon(layerAB_x, layerAB_y, layer2_x, layer2_y, a, l2)
    layer2_near_x2, layer2_near_y2, hexagon_atom2 = search_AA_hexagon(layerAB2_x, layerAB2_y, layer2_x, layer2_y, a, l2)
    AB_area = 3/2*np.sqrt(3)*l1*l1*(len(hexagon_atom)+len(hexagon_atom2))
    #atoms constructing hexagons
    AB_atom = []
    #print(len(hexagon_atom))
    for i in range(len(hexagon_atom)):
        for j in range(3):
            AB_atom.append(hexagon_atom[i][j])
    AB_atom = np.unique(AB_atom)
    
    AB_atom2 = []
    #print(len(hexagon_atom2))
    AB_num = len(hexagon_atom) + len(hexagon_atom2)
    print(AB_num)
    for i in range(len(hexagon_atom2)):
        for j in range(3):
            AB_atom2.append(hexagon_atom2[i][j])
    AB_atom2 = np.unique(AB_atom2)
    
    del hexagon_atom
    layer2_hexagon_x = []
    layer2_hexagon_y = []
    for i in range(len(AB_atom2)):
        layer2_hexagon_x.append(layer2_near_x[AB_atom[i]])
        layer2_hexagon_y.append(layer2_near_y[AB_atom[i]])
    for i in range(len(AB_atom2)):
        layer2_hexagon_x.append(layer2_near_x2[AB_atom2[i]])
        layer2_hexagon_y.append(layer2_near_y2[AB_atom2[i]])
    
    # plt after final selection
    #fig3 = plt.figure(figsize=(16,10))
    #plt.scatter(layer1_x, layer1_y,s = 5)
    #plt.scatter(layer2_x, layer2_y,s = 5, c='orange')
    plt.scatter(layer2_hexagon_x, layer2_hexagon_y,s = 5, c='blue')
    #plt.xlim(-20,10)
    #plt.ylim(10,34)
    #plt.tight_layout()
    plt.axis("equal") 
    plt.show()        
    print("the area of AB stacking is %.5f A^2"%AB_area)            
    return AB_area, AB_num

def calc_AB_stacking_area_drttws(angle, atom_layer1, atom_layer2, cell_parameter1, cell_parameter2, l1, l2, AB_atom_1, AB_atom_2):
    """
    Calculate AB stacking hexagon

    :param atom_layer1, atom_layer2: [][] float64 array
        coordinates of the layer1 and layer2 atoms(reduce same x,y but different z atoms)
    :param cell_parameter1, cell_parameter2: [3][3] float64 array
        cell_parameter of the supercell in layer1 and 2
    :param l: float
        the side length of hexagon
    :
    :return: AB_area
        area of AB stacking
    """
    a = 0.6#searching length
    error = 0.01
    #六边形边长
    atom_num_1 = len(atom_layer1)
    atom_num_2 = len(atom_layer2)
    atom_layer1 = np.array(atom_layer1,dtype=np.float64)
    atom_layer2 = np.array(atom_layer2,dtype=np.float64)
    AB_atom_1 = np.array(AB_atom_1,dtype=np.float64)
    AB_atom_2 = np.array(AB_atom_2,dtype=np.float64)
    cell_parameter1 = np.array(cell_parameter1)
    cell_parameter2 = np.array(cell_parameter2)
    # direct to Cartesian
    for i in range(atom_num_1):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer1[i],cell_parameter1)
        atom_layer1[i] = atom_tmp
    for j in range(atom_num_2):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer2[j],cell_parameter2)
        atom_layer2[j] = atom_tmp
    for j in range(len(AB_atom_1)):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(AB_atom_1[j],cell_parameter1)
        AB_atom_1[j] = atom_tmp
    for j in range(len(AB_atom_2)):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(AB_atom_2[j],cell_parameter1)
        AB_atom_2[j] = atom_tmp
    # plt
    layer1_x = np.zeros(atom_num_1)
    layer1_y = np.zeros(atom_num_1)
    layer2_x = np.zeros(atom_num_2)
    layer2_y = np.zeros(atom_num_2)
    layerAB_1_x = np.zeros(len(AB_atom_1))
    layerAB_1_y = np.zeros(len(AB_atom_1))
    layerAB_2_x = np.zeros(len(AB_atom_2))
    layerAB_2_y = np.zeros(len(AB_atom_2))
    
    for i in range(atom_num_1):
        layer1_x[i] = atom_layer1[i][0]
        layer1_y[i] = atom_layer1[i][1]
    for i in range(atom_num_2):
        layer2_x[i] = atom_layer2[i][0]
        layer2_y[i] = atom_layer2[i][1] 
    for i in range(len(AB_atom_1)):
        layerAB_1_x[i] = AB_atom_1[i][0]
        layerAB_1_y[i] = AB_atom_1[i][1]
    for i in range(len(AB_atom_2)):
        layerAB_2_x[i] = AB_atom_2[i][0]
        layerAB_2_y[i] = AB_atom_2[i][1]
#    fig3 = plt.figure(figsize=(16,10))
#    plt.scatter(layer1_x, layer1_y,s = 20)
    #plt.scatter(layerAB_2_x, layerAB_2_y,s = 5, c='g')
#    #plt.scatter(layer2_hexagon_x, layer2_hexagon_y,s = 5, c='r')
#    #plt.xlim(-120,180)
#    #plt.ylim(-60,180)
#    plt.tight_layout()
#    plt.show()          

    # step 1: find atoms in layer2 that are close to atoms in layer 1
    del atom_layer1
    del atom_layer2
    layer2_near_x, layer2_near_y, hexagon_atom = search_AA_hexagon(layerAB_1_x, layerAB_1_y, layer2_x, layer2_y, a, l2)
    layer2_near_x2, layer2_near_y2, hexagon_atom2 = search_AA_hexagon(layerAB_2_x, layerAB_2_y, layer2_x, layer2_y, a, l2)
    #print(hexagon_atom)
    #print(hexagon_atom2)
    AB_area = 3/2*np.sqrt(3)*l2*l2*(len(hexagon_atom)+len(hexagon_atom2))
    #atoms constructing hexagons
    AB_atom = []
    for i in range(len(hexagon_atom)):
        for j in range(3):
            AB_atom.append(hexagon_atom[i][j])
    AB_atom = np.unique(AB_atom)
    AB_atom2 = []
    AB_num = len(hexagon_atom) + len(hexagon_atom2)
    print(AB_num)
    for i in range(len(hexagon_atom2)):
        for j in range(3):
            AB_atom2.append(hexagon_atom2[i][j])
    AB_atom2 = np.unique(AB_atom2)
    del hexagon_atom
    layer2_hexagon_x = []
    layer2_hexagon_y = []
    for i in range(len(AB_atom)):
        layer2_hexagon_x.append(layer2_near_x[AB_atom[i]])
        layer2_hexagon_y.append(layer2_near_y[AB_atom[i]])
    for i in range(len(AB_atom2)):
        layer2_hexagon_x.append(layer2_near_x2[AB_atom2[i]])
        layer2_hexagon_y.append(layer2_near_y2[AB_atom2[i]])
        
        
    AB = []
    for i in range(len(layer2_hexagon_x)):
        AB.append([layer2_hexagon_x[i],layer2_hexagon_y[i],0])
    cartesian = [[1.0e0,0,0],
                 [0,1.0e0,0],
                 [0,0,1.0e0]]
    AB_frac=[]
    for i in range(len(AB)):
        AB_tmp = transform_coord(cartesian,AB[i],cell_parameter2)
        AB_frac.append(AB_tmp)
    output_AB = open("./AB.txt",mode='w')
    for i in range(len(AB_frac)):
        for j in range(3):
            print("%20.10f"%(AB_frac[i][j]),file=output_AB, end = "")
        print("",file=output_AB)
            
    # plt after final selection
    #fig3 = plt.figure(figsize=(16,10))
    #plt.scatter(layer1_x, layer1_y,s = 5)
    #plt.scatter(layer2_x, layer2_y,s = 5, c='orange')
    plt.scatter(layer2_hexagon_x, layer2_hexagon_y,s = 5, c='blue')
    #plt.xlim(420,520)
    #plt.ylim(0,100)
    picname = str(angle) + '° AA and AB stacking'
    plt.title(picname)
    plt.tight_layout()  
    picfile = './pic/AB' + str(angle) + '°.png'
    plt.savefig(picfile, bbox_inches='tight')      
    plt.show()
    print("the area of AB stacking is %.5f A^2"%AB_area)            
    return AB_area, AB_num

def period_search(angle, atom_layer1, atom_layer2, cell_parameter1, cell_parameter2, l1, l2):
    """
    Search for periodic cell of stacking

    :param atom_layer1, atom_layer2: [][] float64 array
        coordinates of the layer1 and layer2 atoms(reduce same x,y but different z atoms)
    :param cell_parameter1, cell_parameter2: [3][3] float64 array
        cell_parameter of the supercell in layer1 and 2
    :param l1, l2: float
        the side length of hexagon1,2
    :
    :return:
    """
    atom_num_1 = len(atom_layer1)
    atom_num_2 = len(atom_layer2)
    atom_layer1 = np.array(atom_layer1,dtype=np.float64)
    atom_layer2 = np.array(atom_layer2,dtype=np.float64)
    cell_parameter1 = np.array(cell_parameter1)
    cell_parameter2 = np.array(cell_parameter2)
    # direct to Cartesian
    for i in range(atom_num_1):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer1[i],cell_parameter1)
        atom_layer1[i] = atom_tmp
    for j in range(atom_num_2):
        atom_tmp = np.array([0.0,0.0,0.0],dtype=np.float64)
        atom_tmp = np.matmul(atom_layer2[j],cell_parameter2)
        atom_layer2[j] = atom_tmp
    #
    layer1_x = np.zeros(atom_num_1)
    layer1_y = np.zeros(atom_num_1)
    layer2_x = np.zeros(atom_num_2)
    layer2_y = np.zeros(atom_num_2)
    for i in range(atom_num_1):
        layer1_x[i] = atom_layer1[i][0]
        layer1_y[i] = atom_layer1[i][1]
    for i in range(atom_num_2):
        layer2_x[i] = atom_layer2[i][0]
        layer2_y[i] = atom_layer2[i][1] 
        
    a = 0.1
    # search for AA
    # step 1: find atoms in layer2 that are close to atoms in layer 1
    #arrange the x and y coordinate of atoms in each layer
    layer1_argx = np.argsort(layer1_x)
    #layer1_argy = np.argsort(layer1_y)
    layer2_argx = np.argsort(layer2_x)
    #layer2_argy = np.argsort(layer2_y)
    # select the atoms in layer2 which are close to each atom in layer1
    atom_near = []
    r_set = []
    nlable = 0 # label used towards searching for layer2 atom
    for i in range(len(layer1_argx)):
        atom_near_i = []
        r_set_i = []
        # search from x-a to x+a
        while((layer2_x[layer2_argx[nlable]]>=(layer1_x[layer1_argx[i]]-a))):
            if (nlable != 0):
                nlable -= 1
            else:
                break
        while((layer2_x[layer2_argx[nlable]]<(layer1_x[layer1_argx[i]]-a))):
            nlable += 1
            if (nlable >= atom_num_1)or(nlable >= atom_num_2):
                nlable -= 1
                break
    #now we start searching from x-a
        while(layer2_x[layer2_argx[nlable]]<(layer1_x[layer1_argx[i]]+a)):
            r12 = np.sqrt(np.power((layer2_x[layer2_argx[nlable]]-layer1_x[layer1_argx[i]]),2)
    +np.power((layer2_y[layer2_argx[nlable]]-layer1_y[layer1_argx[i]]),2))
            if (r12 <= a):
                #print(layer2_argx[nlable])
                atom_near_i.append(layer2_argx[nlable])
                r_set_i.append(r12)
                break
            nlable += 1
            if (nlable >= atom_num_1)or(nlable >= atom_num_2):
                nlable -= 1
                break
        atom_near.append(atom_near_i)
        r_set.append(r_set_i)
    #
    layer2_near_x=[]
    layer2_near_y=[]
    near_r = []
    for i in range(len(atom_near)):
        if (atom_near[i]):
#        and((layer2_x[atom_near[i][0]] not in layer2_near_x)
#    or(layer2_y[atom_near[i][0]] not in layer2_near_y))):
            #for j in range(len(atom_near[i])):
            layer2_near_x.append(layer2_x[atom_near[i][0]])
            layer2_near_y.append(layer2_y[atom_near[i][0]])
            near_r.append(r_set[i][0])
    #find the closet 8 points
    close_point = np.argsort(near_r)
    close_pointx = []
    close_pointy = []
    n_src = 12
    for i in range(n_src):
        close_pointx.append(layer2_near_x[close_point[i]]) 
        close_pointy.append(layer2_near_y[close_point[i]])
        print(near_r[close_point[i]])
    fig1 = plt.figure(figsize=(16,10))
    plt.scatter(layer1_x, layer1_y,s = 5, c='orange')
    plt.scatter(layer2_x, layer2_y,s = 5, c='pink')
    plt.scatter(close_pointx, close_pointy,s = 5, c='r')
    #plt.xlim(-120,180)
    #plt.ylim(-60,180)
    plt.axis("equal") 
    plt.show()
    #find out the 4 points that makes a parallelogram    
    for i_3 in range(3, n_src):
        for i_2 in range(2, i_3):
            for i_1 in range(1, i_2):
                for i_0 in range(0, i_1):
                    point_set = [i_0, i_1, i_2, i_3]
                    tmp = []
                    # find the atom closest to O 
                    for n_tmp in range(4):
                        tmp.append(close_pointx[point_set[n_tmp]]**2+close_pointy[point_set[n_tmp]]**2)
                    i_min = np.argmin(tmp)
                    vertex0 = point_set[i_min]
                    point_set_2 = []
                    for m in range(4):
                        if (m != i_min):
                            point_set_2.append(point_set[m])
                    #
                    for n in range(3):
                        point_set_3 = []
                        for k in range(3):
                            if (k != n):
                                point_set_3.append(point_set_2[k])
                        if (( np.abs(close_pointx[point_set_3[0]]+close_pointx[point_set_3[1]]-close_pointx[vertex0]-close_pointx[point_set_2[n]])<=3 )
                        and( np.abs(close_pointy[point_set_3[0]]+close_pointy[point_set_3[1]]-close_pointy[vertex0]-close_pointy[point_set_2[n]])<=3 )):
                            print((close_pointx[point_set_3[0]]+close_pointx[point_set_3[1]]-close_pointx[vertex0]-close_pointx[point_set_2[n]]))
                            print(close_pointx[point_set_3[0]])
                            print(close_pointx[point_set_3[1]])
                            print(close_pointx[vertex0])
                            print(close_pointx[point_set_2[n]])
                            if ((((close_pointx[vertex0]-close_pointx[point_set_3[0]])**2+
                                (close_pointy[vertex0]-close_pointy[point_set_3[0]])**2-
                                (close_pointx[vertex0]-close_pointx[point_set_3[1]])**2-
                                (close_pointy[vertex0]-close_pointy[point_set_3[1]])**2)<=1)and
                                (((close_pointx[vertex0]-close_pointx[point_set_3[0]])**2+
                                (close_pointy[vertex0]-close_pointy[point_set_3[0]])**2-
                                (close_pointx[point_set_3[0]]-close_pointx[point_set_3[1]])**2-
                                (close_pointy[point_set_3[0]]-close_pointy[point_set_3[1]])**2)<=3)):
                            
                                parallelogramx = [close_pointx[vertex0], close_pointx[point_set_3[0]], close_pointx[point_set_3[1]], close_pointx[point_set_2[n]]]
                                parallelogramy = [close_pointy[vertex0], close_pointy[point_set_3[0]], close_pointy[point_set_3[1]], close_pointy[point_set_2[n]]]
    fig2 = plt.figure(figsize=(16,10))
    plt.scatter(layer1_x, layer1_y,s = 5, c='orange')
    plt.scatter(layer2_x, layer2_y,s = 5, c='pink')
    plt.scatter(parallelogramx, parallelogramy,s = 15, c='r')
    print(parallelogramx)
    print(parallelogramy)
    #plt.xlim(-120,180)
    #plt.ylim(-60,180)
    plt.axis("equal")  
    plt.show()
                    
    return 0
    
    
    
    
    
    
    
    
    
    
    
    