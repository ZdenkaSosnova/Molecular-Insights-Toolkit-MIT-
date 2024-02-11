from utils import calculate_lengt
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

class Molecule_constructor:
    def __init__(self, length_inside, length_between, number_of_benzene_rings, repetition_count, file_name, **kwargs):
        '''
        :param length_inside (float): bond length with elementary cell
        :param length_between (float): bond length between elementary cells
        :param number_of_benzene_rings (odd int):  number of benzene rings within elementary cell
        :param repetition_count (int): number of elementary cell repetition
        :param file_name (str): name of the file .xyz format
        :param kwargs: defining other elementary cell of the final molecule
                    {
                    key: [length_inside (float), length_between (float), number_of_benzene_rings (int), position within molecule (list)
                    }
        '''
        dictionary_kwarg = Molecule_constructor.check_kwargs_inputs(repetition_count, kwargs)
        dictionary_elementary_cells, dictionary_distances = Molecule_constructor.create_element_cell_database(length_inside, length_between, number_of_benzene_rings, repetition_count, dictionary_kwarg)
        self.my_molecule = Molecule_constructor.coordinate_of_molecule(dictionary_elementary_cells, dictionary_distances, repetition_count)
        self.file_name = file_name
        self.v_min, self.v_max = Molecule_constructor.min_max_distance(length_inside, length_between, kwargs)


    @staticmethod
    def min_max_distance(length_1, length_2, dictionary):
        '''
        Computes the minimum and maximum bond length within the final molecule
        :return: tuple: Minimum and maximum distance
        '''
        lengths = [dictionary[element][0] for element in dictionary] + [dictionary[element][1] for element in dictionary] + [length_1, length_2]
        return np.min(lengths) - 0.1, np.max(lengths) + 0.1


    @staticmethod
    def create_element_cell_database(inner_length, joining_length, benzene_count, repetition_count, dictionary):
        '''
        Creates a database of elementary cells
        :return: tuple: Elementary cells database and distance dictionary
        '''
        if not isinstance(benzene_count,int):
            raise ValueError("Incorectly defined variable number of rings.")
        if benzene_count % 2 == 0 or benzene_count <= 0:
            raise ValueError("Incorectly defined variable number of rings.")
        only_one_cell = [x for x in range(repetition_count)]
        database = {}
        period_dictionary = {}
        for element in only_one_cell:
            database[element] = Molecule_constructor.create_elementary_cell(inner_length, benzene_count)
            period_dictionary[element] = 4 * inner_length + joining_length
        for value in dictionary.values():
            for order in value[3]:
                database[order] = Molecule_constructor.create_remaining_elementary_cells(value[0], value[2], -2 * inner_length,
                                                                             math.sqrt((inner_length * inner_length) - (inner_length / 2) * (
                                                                            inner_length / 2)) * benzene_count)
                period_dictionary[order] = 4 * value[0] + value[1]
        return database, period_dictionary

    @staticmethod
    def check_kwargs_inputs(repetition_count, dictionary):
        '''
        Checks if kwargs inputs were given correctly
        :return: Dictionary (kwargs information)
        '''
        for element in dictionary:
            if isinstance(dictionary[element], list):
                if len(dictionary[element]) == 4 and isinstance(dictionary[element][0], (int, float)) and isinstance(dictionary[element][1], (int,float)) and isinstance(dictionary[element][2],int) and isinstance(dictionary[element][3],list):
                    if dictionary[element][2] % 2 == 0:
                        raise ValueError(f"Incorrectly defined variable corresponding to the key {element}.")
                    for j in range(len(dictionary[element][3])):
                        if not isinstance(dictionary[element][3][j], int):
                            raise ValueError(f"Incorrectly defined variable corresponding to the key {element}.")
                else:
                    raise ValueError(f"Incorrectly defined variable corresponding to the key {element}.")
            else:
                raise ValueError(f"Incorrectly defined variable corresponding to the key {element}.")
            dictionary[element][3] = [i-1 for i in dictionary[element][3]]
        overlapping_check_lists = [dictionary[element][3][i] for element in dictionary for i in range(len(dictionary[element][3]))]
        overlapping_value_check_lists = [i for i in overlapping_check_lists if i < repetition_count + 1 and i >= 0]
        if not len(set(overlapping_check_lists)) == len(overlapping_check_lists):
            raise ValueError("Positions of other elementary cells overlap")
        if not len((overlapping_check_lists)) == len(overlapping_value_check_lists):
            raise ValueError("The requested position of the elementary cell is not possible with respect to other input parameters (the position is negative or exceeds the length of the molecule)")
        return dictionary

    @staticmethod
    def create_elementary_cell(bond_length, benzene_count):
        '''
        Creates a 'main' elementary cell
        :return: list: Elementary cell
        '''
        base = []
        base.append([-bond_length / 2, 0.0])
        base.append([bond_length / 2, 0.0])
        base.append([-bond_length, math.sqrt((bond_length * bond_length) - (bond_length / 2) * (bond_length / 2))])
        base.append([bond_length, math.sqrt((bond_length * bond_length) - (bond_length / 2) * (bond_length / 2))])
        elementary_cell = [i for i in base]
        for k in range(len(base)):
            for i in range(1, benzene_count + 1):
                element = float(base[k][1]) + i * 2 * math.sqrt(
                    (bond_length * bond_length) - (bond_length / 2) * (bond_length / 2))
                if element <= (benzene_count * 2 + 0.5) * math.sqrt(
                        (bond_length * bond_length) - (bond_length / 2) * (bond_length / 2)):
                    elementary_cell.append([base[k][0], element])
        elementary_cell.append([2 * bond_length, benzene_count * math.sqrt(
            (bond_length * bond_length) - (bond_length / 2) * (bond_length / 2))])
        elementary_cell.append([-2 * bond_length, benzene_count * math.sqrt(
            (bond_length * bond_length) - (bond_length / 2) * (bond_length / 2))])
        return elementary_cell

    @staticmethod
    def create_remaining_elementary_cells(bond_length_2, benzene_count_2, x_coordinate, center):
        '''
        Creates remaining elementary cells
        :param x_coordinate: shift on x axis according to 'main' elementary cell
        :param center: shift one y axis according to 'main' alementary cell
        :return: List of other elementary cell
        '''
        elementary_cell_2 = Molecule_constructor.create_elementary_cell(bond_length_2, benzene_count_2)
        y_shift = center - benzene_count_2 * math.sqrt((bond_length_2 * bond_length_2) - (bond_length_2 / 2) * (bond_length_2 / 2))
        x_shift = x_coordinate + 2 * bond_length_2
        for i in range(len(elementary_cell_2)):
            elementary_cell_2[i][1] = elementary_cell_2[i][1] + y_shift
            elementary_cell_2[i][0] = elementary_cell_2[i][0] + x_shift
        return elementary_cell_2

    @staticmethod
    def coordinate_of_molecule(cell_units, period, molecule_length):
        '''
        Computes a cartesian coordinates of molecule defined by inputs
        :return: list: Molecule coordinates
        '''
        translation = 0
        molecule = []
        for k in range(molecule_length):
            for i in range(len(cell_units[k])):
                molecule.append(["C", float(cell_units[k][i][0])+translation, float(cell_units[k][i][1]),0])
            translation += period[k]
        min_x = min(range(len(molecule)), key=lambda i: molecule[i][1])
        molecule.pop(min_x)
        max_x = max(range(len(molecule)), key=lambda i: molecule[i][1])
        molecule.pop(max_x)
        return molecule


    def molecule_coordinates_to_xyz_file(self):
        '''
        Writes the molecule cartesian coordinates into .xyz file
        '''
        length = len(self.my_molecule)
        try:
            f = open(f"{self.file_name.split('.')[0]}.xyz", "w")
            f.write(f"{length}\n\n")
            for a in range(len(self.my_molecule)):
                for i in range(4):
                    f.write(f"{self.my_molecule[a][i]}    ")
                f.write("\n")
            f.close()
        except:
            f = open(f"{self.file_name('.')[0]}.xyz", "a")
            f.write(f"{length}\n\n")
            for a in range(len(self.my_molecule)):
                for i in range(4):
                    f.write(f"{self.my_molecule[a][i]}    ")
                f.write("\n")
            f.close()

    def show_graph(self):
        '''
        Shows the graph of the molecule
        '''
        length_x = calculate_lengt(min(self.my_molecule, key=lambda x: x[1])[1],
                                    max(self.my_molecule, key=lambda x: x[1])[1])
        length_y = calculate_lengt(min(self.my_molecule, key=lambda x: x[2])[2],
                                    max(self.my_molecule, key=lambda x: x[2])[2])
        aspect_ratio = round(length_x / length_y, 1)
        fig = plt.figure(figsize=(aspect_ratio, 1.5))
        axes1 = fig.add_axes([0.0, 0.0, 0.9, 1])
        cmap = cm.get_cmap('cool')
        norm = plt.Normalize(vmin=round(self.v_min, 2), vmax=round(self.v_max, 2))
        length = len(self.my_molecule)
        for i in range(length):
            for k in range(i + 1, length):
                x = abs(self.my_molecule[i][1] - self.my_molecule[k][1])
                y = abs(self.my_molecule[i][2] - self.my_molecule[k][2])
                distance = math.sqrt(x * x + y * y)
                if distance + 0.01 >= 1.10 and distance - 0.01 <= 1.90:
                    x_values = [self.my_molecule[i][1], self.my_molecule[k][1]]
                    y_values = [self.my_molecule[i][2], self.my_molecule[k][2]]
                    color = cmap(norm(distance))
                    axes1.plot(x_values, y_values, color=color)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cax = fig.add_axes([0.9, 0.05, 0.01, 0.9])  # Adjust the position and width of the colorbar axis
        plt.colorbar(sm, label='Bond length', cax=cax)
        axes1.set_aspect('equal')
        axes1.axis('off')
        fig.savefig(f"{self.file_name}.png", dpi=300)
        plt.show()
