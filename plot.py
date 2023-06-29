import matplotlib.pyplot as plt
def plot_atoms(cell_param, atom):
    """
    Plot atoms. 
    :param cell_param [3][3] float array
        cell parameter of primitive cell
    :param atom [N][3] float array
        atom position in fractional coordinates
    """  
    atom_x = []
    atom_y = []
    for i in range(len(atom)):
        x = 0
        y = 0
        for j in range(3): 
            x += atom[i][j]*cell_param[j][0]
            y += atom[i][j]*cell_param[j][1]
        atom_x.append(x)
        atom_y.append(y)
    fig = plt.figure(figsize=(10,10))
    
    plt.plot([0, cell_param[0][0]],[0, cell_param[0][1]], c = 'r')
    plt.plot([0, cell_param[1][0]],[0, cell_param[1][1]], c = 'r')
    plt.scatter(atom_x, atom_y, s = 5)
    plt.axis("equal")
    #plt.show()
    return 0
    