from utils import check_input_validity, load_coordinates, distance_matrix, calculate_lengt
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
from matplotlib import cm

class Huckel_model:
    '''
    Class implementing an approximate calculation of electronic structure and molecular orbitals using the Huckel method
    for 'pi'-conjugated planar (2D - x, y) molecules
    '''
    def __init__(self, file, alfa = 0, beta = -2.8, extended_huckel = False, number_of_states = 0, minimal_distance = 1.10, maximal_distance = 1.60):
        '''
        :param file: File in ".xyz" format specifying the coordinates of the selected molecule
                     the program will only evaluate carbon atoms
        :param alpha: "on-site" energy of the atom - if the molecule contains the same type of atoms - without loss
                      of information, it can be set to 0
                      based on physical experiments - negative or zero value (calculation will proceed even for positive values)
        :param beta: parameter determining the interaction between individual atoms - from the theory of the Huckel method -
                      non-zero parameter only for nearest neighbors
                      based on physical experiments - negative or zero value (calculation will proceed even for positive values of the parameter)
        :param extended_huckel: boolean variable
                                'True' -> implementation of 'extended Huckel method' - the value of beta is parameterized
                                based on the distance between atoms (parameterization introduced by the relationship "beta * (1.4/distance)**2")
        :param state_count: number of states around the Fermi energy to be represented in subsequent analysis -
                            eigenvalues or molecular orbitals
                            the value must be less than or equal to the number of carbon atoms
        :param min_distance: minimum distance between individual (usually neighboring) atoms that I want to visualize on the graph
                             based on physical intuition - minimum value of carbon-carbon bond ~1.15
        :param max_distance: maximum distance between individual (usually neighboring) atoms that I want to visualize on the graph
                             based on physical intuition - maximum bond value ~1.55
        '''
        check_input_validity(alfa, "alfa", (int,float))
        check_input_validity(beta, "beta", (int,float))
        check_input_validity(number_of_states, "pocet_stavu", (int))
        check_input_validity(minimal_distance, "minimalni_vzdalenost", (int,float))
        check_input_validity(maximal_distance, "maximalni_vzdalenost", (int, float))
        check_input_validity(extended_huckel, "extenden_huckel", (bool))
        self.file_name = file
        self.molecule_coordinates = load_coordinates(file, carbon_only= True)
        self.rozmer = len(self.molecule_coordinates)
        self.number_of_states = number_of_states
        self.distance = distance_matrix(self.molecule_coordinates)
        self.eigenvalues, self.eigenvectors = self.create_hamiltonian(alfa, beta, extended_huckel, minimal_distance, maximal_distance)
        self.state_names = Huckel.state_list(self.rozmer)
        self.v_min = minimal_distance
        self.v_max = maximal_distance

    @staticmethod
    def state_list(list_length):
        '''
        Function returns a list that will subsequently serve to describe individual states in output graphs
        Occupied states - from the highest (in terms of energy) occupied state indicated as HOMO ('Highest Occupied Molecular Orbital'), HOMO-1, ...
        Unoccupied states - from the lowest (in terms of energy) unoccupied state indicated as LUMO ('Lowest Unoccupied Molecular Orbital'), LUMO+1, ...
        In the case of the Huckel method - the first half of the lowest eigenstates belongs to occupied states, the second half to unoccupied states
        :param list_length: number of eigenstates belonging to the Hamiltonian
                            for the Huckel method it is the number of carbon atoms
        :return: list serving to describe individual states in output graphs
        '''
        state_names = [f"HOMO-{int(list_length/2 -1 -i)}" if int(-list_length/2 +1 +i) < 0 else(
                        "HOMO" if int(list_length/2 -1 -i) == 0 else(
                        "LUMO" if (int(list_length/2 - i )) == 0 else
                        f"LUMO+{int(-list_length/2 + i)}"))
                        for i in range(int(list_length))]

        return state_names

    def how_many_states_to_draw(self, states):
        '''
        Function checks if the number of states entered by the user can be visualized, in case of an odd number it returns an even number
        '''
        if states % 2 == 1:
            states = states + 1
        if states > len(self.molecule_coordinates):
            states = len(self.molecule_coordinates)
            print(f"I don't have that many eigenstates, I'll show {len(self.molecule_coordinates)}")

        return states

    def create_hamiltonian(self, alfa, beta, extended_huckel, minimal_value, maximal_value):
        '''
        Function creates the Hamiltonian based on the theory of the (Extended) Huckel method, also solves the problem of eigenvalues and states
        :return: eigenvalues and vectors of the constructed Huckel Hamiltonian
        '''
        '''
        Size of Hamiltonian - square matrix, size corresponds to the number of carbon atoms
        '''
        hamiltonian = np.zeros(self.rozmer**2).reshape(self.rozmer, self.rozmer)
        condition_1 = self.distance == 0
        hamiltonian[condition_1] = alfa
        condition_2 = (self.distance >= minimal_value) & (self.distance <= maximal_value)
        if extended_huckel:
            hamiltonian[condition_2] = beta * (1.4/self.distance[condition_2])**2
        else:
            hamiltonian[condition_2] = beta
        '''
        Vyřešení problému vlastních čísel právě zkontruovaného Hamiltoniánu
        '''
        eigenvalues, eigenvectors = scp.linalg.eigh(a=hamiltonian)
        return eigenvalues, eigenvectors

    def energy_graph(self):
        '''
        Representation of energies (eigenvalues) closest to the Fermi energy - number of states in the class parameter
        :return: Energy graph (eigenvalues) of the given Hamiltonian + representation of Fermi energy
        '''
        how_many_states = self.how_many_states_to_draw(self.number_of_states)
        for i in range(int(len(self.eigenvalues)/2)-(int(how_many_states/2)), int(len(self.eigenvalues)/2)+(int(how_many_states/2))):
            plt.plot(self.state_names[i],self.eigenvalues[i], "ro")
        plt.plot([self.state_names[int(len(self.eigenvalues)/2)-(int(how_many_states/2))],self.state_names[int(len(self.eigenvalues)/2)+(int(how_many_states/2))-1]],[0,0], color = "y", label = "Fermi level")
        plt.legend()
        plt.xlabel("Stavy")
        plt.ylabel("Energie [eV]")
        fig.savefig(f"{self.file_name.split('.')[0]}_energy.png")
        plt.show()

    def orbital_graph(self, orbital, state):
        '''
        Function that creates a graph for the selected orbital - resulting in a molecular orbital corresponding to the eigenvalue (energy)
        :param orbital: orbital of interest - type 'int' - orbital corresponding to the lowest energy '0', highest energy 'number of carbon atoms - 1'
        :param state: Graph label - HOMO, LUMO,.. can also be a numerical value
        :return: graph of the selected molecular orbital
        '''
        length_1 = calculate_lengt(np.min(self.molecule_coordinates[:, 0]), np.max(self.molecule_coordinates[:, 0]))
        length_2 = calculate_lengt(np.min(self.molecule_coordinates[:, 1]), np.max(self.molecule_coordinates[:, 1]))
        aspect_ratio = round(length_1/length_2, 1)
        '''
        Determining the markersize of the graph - to ensure visibility of corresponding orbitals (determined by trial and error)
        '''
        marker_size = round(0.1 * len(self.molecule_coordinates) + 10,1)
        if marker_size < 25:
            marker_size = 30
        elif marker_size > 100:
            marker_size = 100
        fig = plt.figure(figsize = (aspect_ratio,1.5))
        ax = fig.add_axes((0.0, 0.0, 1, 1))
        for i in range(self.rozmer):
            for j in range(i+1,self.rozmer):
                if self.distance[i][j] < 1.7:
                    ax.plot([self.molecule_coordinates[:,0][i], self.molecule_coordinates[:,0][j]],[self.molecule_coordinates[:,1][i],self.molecule_coordinates[:,1][j]],"grey")
        for k in range(len(orbital)):
            ax.plot(self.molecule_coordinates[:,0][k],self.molecule_coordinates[:,1][k], "ro" if orbital[k] > 0 else "go", markersize=marker_size * abs(orbital[k]))
        ax.set_title(f"{state}")
        ax.set_aspect("equal")
        ax.axis('off')
        fig.savefig(f"{self.file_name.split('.')[0]}_{state}.png")
        fig.show()

    def huckel_orbitaly(self):
        '''
        :return: Plotting selected number of orbitals around the Fermi energy
        '''
        how_many_states = self.how_many_states_to_draw(self.number_of_states)
        for i in range(int(len(self.molecule_coordinates)/2) - int(how_many_states/2), int(len(self.molecule_coordinates)/2) + int(how_many_states/2)):
            self.orbital_graph(self.eigenvectors[:, i], self.state_names[i])

    def return_gap_value(self):
        '''
        Energy difference between the 'highest' occupied and 'lowest' unoccupied orbital
        Energy difference determines the basic properties of the material in terms of conductivity
        '''
        return self.eigenvalues[(int(len(self.eigenvalues)/2))]-self.eigenvalues[(int(len(self.eigenvalues)/2)-1)]

    def bond_charge(self):
        bond_charge_matrix = np.zeros((len(self.eigenvectors), len(self.eigenvectors)))
        for i in range(0, len(self.eigenvectors)):
            for j in range(0, len(self.eigenvectors)):
                A = np.zeros((len(self.eigenvectors), len(self.eigenvectors)))
                A[i, j] = 1
                for stav in range(0, int(len(self.eigenvectors) / 2)):
                    bond_charge_matrix[i, j] = bond_charge_matrix[i, j] + 2 * np.dot(self.eigenvectors[:, stav], np.dot(A, self.eigenvectors[:, stav]))
        return bond_charge_matrix

    def bong_charge_matrix_txt(self):
        bond_charge_matrix = self.bond_charge()
        np.savetxt(f"{self.file_name.split('.')[0]}_bond_charge.txt", bond_charge_matrix)


    def graph_bond_charge(self):
        bond_charge_matrix = self.bond_charge()
        triangular_matrix = np.triu(self.distance, k=1)
        indices = np.where((triangular_matrix < self.v_max) & (triangular_matrix > self.v_min))
        delka_x = calculate_lengt(min(self.molecule_coordinates, key=lambda x: x[0])[0],
                                 max(self.molecule_coordinates, key=lambda x: x[0])[0])
        delka_y = calculate_lengt(min(self.molecule_coordinates, key=lambda x: x[1])[1],
                                 max(self.molecule_coordinates, key=lambda x: x[1])[1])
        aspect_ratio = round(delka_x / delka_y, 1)
        # Extract values from the original matrix corresponding to the indices
        bond_charge_values = []
        for (i, j) in zip(indices[0], indices[1]):
            bond_charge_values.append(bond_charge_matrix[i, j])
        min_value = np.min(bond_charge_values)
        max_value = np.max(bond_charge_values)
        cmap = cm.get_cmap('hot')
        norm = plt.Normalize(vmin=min_value - 0.05, vmax=max_value + 0.15)
        fig = plt.figure(figsize=(2*aspect_ratio, 2*0.8), dpi=500)
        ax = fig.add_axes([0.0, 0.0, 0.8, 1])
        for (i, j) in zip(indices[0], indices[1]):
            x1, y1 = self.molecule_coordinates[:, 0][i], self.molecule_coordinates[:, 1][i]
            x2, y2 = self.molecule_coordinates[:, 0][j], self.molecule_coordinates[:, 1][j]
            color_text = cmap(norm(bond_charge_matrix[i, j]))
            bond_charge_i_j = round(bond_charge_matrix[i][j], 3)
            ax.plot([x1, x2], [y1, y2], color="grey", alpha=.3, linewidth=1)
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(mid_x, mid_y, f'{bond_charge_i_j:.2f}', fontsize=10, ha='center', va='center', color=color_text,
                    fontweight="bold")
        ax.set_aspect("equal")
        ax.axis('off')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cax = fig.add_axes([0.8, 0.2, 0.03, 0.6])
        colorbar = plt.colorbar(sm, label='Délka vazby', cax=cax)
        colorbar.ax.tick_params(axis='y', labelsize=10)
        fig.savefig(f"{self.file_name.split('.')[0]}_bond_charge.png", dpi=500)
        fig.show()



