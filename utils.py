import re
import sys
import numpy as np


periodic_table = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ar', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Ni', 'Co', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg',
    'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv',
    'Ts', 'Og'
]

def modify_rows(string):
    return re.split(r'\s+', string)

def in_file_to_xyz(file):
    '''
    Convert geometry from ".in" format (FHI Aims) to ".xyz" format
    :param file: file in ".in" format
    :return: coordinates in ".xyz" format - the name remains the same as the original file with a change in extension
    '''
    coordinates_list = []
    try:
        with open(file, "r") as f:
            lines = f.readlines()
        for line in lines:
            if "atom" in line:
                # Split the line based on any whitespace
                values = re.split(r'\s+', line.strip())
                coordinates_list.append(values)
        coordinates_list_xyz = [[sublist[-1], float(sublist[1]), float(sublist[2]), float(sublist[3])] for sublist in coordinates_list]
        try:
            f = open(f"{file.split('.')[0]}.xyz", "w")
            f.write(f"{len(coordinates_list_xyz)}\n")
            f.write("\n")
            for a in range(len(coordinates_list_xyz)):
                for i in range(4):
                    f.write(f"{coordinates_list_xyz[a][i]}    ")
                f.write("\n")
            f.close()
        except:
            f = open(f"{file.split('.')[0]}.xyz", "a")
            f.write(f"{len(coordinates_list_xyz)}\n")
            f.write("\n")
            for a in range(len(coordinates_list_xyz)):
                for i in range(4):
                    f.write(f"{coordinates_list_xyz[a][i]}    ")
                f.write("\n")
            f.close()
    except FileNotFoundError:
        print("FileNotFoundError")
        sys.exit()

def calculate_lengt(min_value, max_value):
    '''
    :param min_value: min. value (x/y/z) coordinate, which we are interested in determining the length
    :param max_value: max. value (x/y/z) coordinate, which we are interested in determining the length
    :return: total length in Cartesian coordinates between min and max value
    '''
    if np.sign(min_value) * np.sign(max_value) < 0:
        return int(round(abs(min_value) + abs(max_value), 0))
    else:
        return int(abs(round(min_value - max_value, 0)))

def load_coordinates(file, carbon_only):
    '''
    Funkce přečte soubor ve formátu ".xyz" na 2d np.array - pouze pro planární molekuly -
    mající nenulové hodnoty pouze v souřadnicích x,y
    Vrátí FileNotFoundError - není-li funkce definována
    Je-li definována "špatně" - nejedná se o planární molekulu - výpočet proběhne se zanedbáním molekuly 'z'
                              - nekontroluje, zda se nejedná o jiný typ formátu
    :param soubor: soubor ve formátu ".xyz"
    :return: 2d numpy array obsahující x,y souřadnice molekuly
    '''
    coordinates_list = []
    i = 0
    try:
        with open(file, "r") as f:
            for line in f:
                print(i)
                print(line)
                if not carbon_only:
                    for element in periodic_table:
                        if element in line:
                            coordinates_list.append(modify_rows(line.strip()))
                else:
                    if "C" in line:
                        coordinates_list.append(modify_rows(line.strip()))
                i += 1
        cleaned_list_of_lists = [[item.strip() for item in sublist if item.strip() != ''] for sublist in coordinates_list]
        coordinates_list = [sublist for sublist in cleaned_list_of_lists if sublist]
    except FileNotFoundError:
        print("FileNotFoundError")
        sys.exit()
    coordinates = np.array([float(coordinates_list[i][j+1]) for i in range(len(coordinates_list)) for j in range(2)]).reshape(len(coordinates_list),2)
    return coordinates

def load_coordinates_3d(file, carbon_only):
    '''
        Function reads a file in ".xyz" format into a 3d np.array - returning the x, y, z coordinates
        Raises FileNotFoundError if the function is not defined
        :param file: file in ".xyz" format
        :return: 3d numpy array containing x, y, z coordinates of the molecule
        '''
    coordinates_list = []
    try:
        with open(file, "r") as f:
            for line in f:
                if not carbon_only:
                    for element in periodic_table:
                        if element in line:
                            coordinates_list.append(modify_rows(line.strip()))
                else:
                    if "C" in line:
                        coordinates_list.append(modify_rows(line.strip()))
        cleaned_list_of_lists = [[item.strip() for item in sublist if item.strip() != ''] for sublist in coordinates_list]
        # Remove sublists that became empty after removing ' '
        coordinates_list = [sublist for sublist in cleaned_list_of_lists if sublist]
    except FileNotFoundError:
        print("FileNotFoundError")
        sys.exit()
    coordinates = np.array([float(coordinates_list[i][j+1]) for i in range(len(coordinates_list)) for j in range(3)]).reshape(len(coordinates_list),3)
    return coordinates

def check_input_validity(parameter_value, parameter, variable_type):
    '''
    Function checks if the input parameter has the correct format
    :param parameter_value: value of the parameter entered by the user
    :param parameter: name of the parameter being checked
    :param variable_type: type of variable (tuple of variable types) that the input parameter must satisfy
    :return: if the input parameter is incorrect - ValueError
    '''
    if not isinstance(parameter_value, variable_type):
        raise ValueError(f"Parameter '{parameter}' must be of type '{variable_type}'")

def distance_matrix(coordinates_matrix):
    '''
    Function returns a matrix that determines the distance (in Cartesian coordinates) between atom i,j
    :param coordinates_matrix: matrix containing information about x, y coordinates (in subsequent use) of the molecule
    :return: matrix of distances between atoms i,j (for all i,j from the coordinates matrix)
    '''
    x = np.abs(coordinates_matrix[:,0, None]-coordinates_matrix[:,0])
    y = np.abs(coordinates_matrix[:,1, None]-coordinates_matrix[:,1])
    distance = np.sqrt(x**2 + y**2)
    return distance

def distance_matrix_3d(coordinates_matrix):
    '''
        Function returns a matrix that determines the distance (in Cartesian coordinates) between atom i,j
        :param coordinates_matrix: matrix containing information about x, y, z coordinates (in subsequent use) of the molecule
        :return: matrix of distances between atoms i,j (for all i,j from the coordinates matrix)
        '''
    x = np.abs(coordinates_matrix[:,0, None]-coordinates_matrix[:,0])
    y = np.abs(coordinates_matrix[:,1, None]-coordinates_matrix[:,1])
    z = np.abs(coordinates_matrix[:, 2, None] - coordinates_matrix[:, 2])
    distance = np.sqrt(x**2 + y**2 + z**2)
    return distance

