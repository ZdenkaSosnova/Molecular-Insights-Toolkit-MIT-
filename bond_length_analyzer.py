import utils
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

class Bond_lenght_Analyzer:
    '''
    Class for analyzing bond lengths in a molecule based on relaxed coordinates
    '''
    def __init__(self, file_xyz, dimension = 2, minimal_distance = 1.35, maximal_distance = 2):
        '''
        :param file_xyz: File containing molecule coordinates. Accepted formats: ".xyz" or ".in".
        :param minimal_distance: Minimum bond length between adjacent atoms.
                                 For carbon-carbon bonds, the minimum value is ~1.15 Angstrom.
                                 Different rules apply for bonds between other atoms (e.g., carbon-hydrogen ~1 Angstrom).
        :param maximal_distance: Maximum bond length between adjacent atoms.
                                 For carbon-carbon bonds, the maximum value is ~1.55 Angstrom.
        '''
        self.file = file_xyz
        if dimension == 2:
            self.molecule_coordinates = utils.load_coordinates(file=file_xyz, carbon_only=False)
            self.distance_matrix = utils.distance_matrix(self.molecule_coordinates)
        elif dimension == 3:
            self.molecule_coordinates = utils.load_coordinates_3d(file=file_xyz, carbon_only=False)
            self.distance_matrix = utils.distance_matrix_3d(self.molecule_coordinates)
        self.v_min = minimal_distance
        self.v_max = maximal_distance


    def graph_2d(self):
        '''
        Generates a 2D graph depicting bond lengths in the molecule defined by coordinates in ".xyz" format.
        Different bond lengths are represented with different colors, and each bond is labeled with its length.
        '''
        triangular_matrix = np.triu(self.distance_matrix, k=1)
        indices = np.where((triangular_matrix < self.v_max) & (triangular_matrix > self.v_min))
        length_x = utils.calculate_lengt(min(self.molecule_coordinates, key=lambda x: x[0])[0],
                                 max(self.molecule_coordinates, key=lambda x: x[0])[0])
        length_y = utils.calculate_lengt(min(self.molecule_coordinates, key=lambda x: x[1])[1],
                                 max(self.molecule_coordinates, key=lambda x: x[1])[1])
        aspect_ratio = round(length_x / length_y, 1)
        fig = plt.figure(figsize=(1.2*aspect_ratio, 0.8), dpi=250)
        ax = fig.add_axes([0.0, 0.0, 0.8, 1])
        cmap = cm.get_cmap('cool')
        norm = plt.Normalize(vmin=round(self.v_min, 2), vmax=round(self.v_max, 2))
        for (i,j) in zip(indices[0],indices[1]):
            x1, y1 = self.molecule_coordinates[:, 0][i], self.molecule_coordinates[:, 1][i]
            x2, y2 = self.molecule_coordinates[:, 0][j], self.molecule_coordinates[:, 1][j]
            bond_length = self.distance_matrix[i][j]
            color_vazba = cmap(norm(bond_length))
            ax.plot([x1, x2], [y1, y2], color=color_vazba)
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(mid_x, mid_y, f'{bond_length:.2f}', fontsize=2, ha='center', va='center')
        ax.set_aspect("equal")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cax = fig.add_axes([0.8, 0.2, 0.02, 0.6])
        colorbar = plt.colorbar(sm, label='Délka vazby', cax=cax)
        colorbar.set_label("Délka vazby", fontsize=4)
        colorbar.ax.tick_params(axis='y', labelsize=4)
        ax.axis('off')
        plt.savefig(f"{self.file.split('.')[0]}_bond_length.png")
        fig.show()

    def projection_y_z_axis(self):
        '''
        Generates 3D projections of bond lengths analysis viewed from the Y-axis and Z-axis.
        '''
        triangular_matrix = np.triu(self.distance_matrix, k=1)
        indices = np.where((triangular_matrix < self.v_max) & (triangular_matrix > self.v_min))
        fig = plt.figure(figsize=(8, 4), dpi=250)
        ax2 = fig.add_subplot(121, projection='3d')
        ax3 = fig.add_subplot(122, projection='3d')
        cmap = cm.get_cmap('cool')
        norm = plt.Normalize(vmin=round(self.v_min, 2), vmax=round(self.v_max, 2))
        for (i, j) in zip(indices[0], indices[1]):
            x1,y1,z1 = self.molecule_coordinates[:,0][i], self.molecule_coordinates[:, 1][i], self.molecule_coordinates[:, 2][i]
            x2,y2,z2 = self.molecule_coordinates[:,0][j], self.molecule_coordinates[:, 1][j], self.molecule_coordinates[:, 2][j]
            bond_length = self.distance_matrix[i][j]
            bond_color = cmap(norm(bond_length))
            ax2.plot([x1, x2], [y1, y2], [z1, z2], color=bond_color)
            ax3.plot([x1, x2], [y1, y2], [z1, z2], color=bond_color)
            mid_x, mid_y, mid_z = (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2
            ax3.text(mid_x, mid_y, mid_z, f'{bond_length:.2f}', color='black', fontsize=8)
        ax2.view_init(elev=0, azim=90)
        ax3.view_init(elev=90, azim=0)
        ax2.set_title('View from Y-axis')
        ax3.set_title('View from Z-axis')
        ax2.axis('off')
        ax2.set_aspect("equal")
        ax3.set_aspect("equal")
        ax3.axis('off')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cax = fig.add_axes([0.92, 0.2, 0.02, 0.6])  # Adjust the position of the color bar
        colorbar = plt.colorbar(sm, label='Bond length', cax=cax)
        colorbar.set_label("Bond length", fontsize=10)
        colorbar.ax.tick_params(axis='y', labelsize=8)
        plt.subplots_adjust(wspace=0)
        plt.savefig(f"{self.file.split('.')[0]}bond_length_projection.png")
        plt.show()

    def graf_3d(self):
        '''
        Generates a 3D graph depicting bond lengths in the molecule.
        '''
        triangular_matrix = np.triu(self.distance_matrix, k=1)
        indices = np.where((triangular_matrix < self.v_max) & (triangular_matrix > self.v_min))
        fig = plt.figure(figsize=(6, 5),dpi=300)
        ax = fig.add_subplot(111, projection='3d')
        cmap = cm.get_cmap('cool')
        norm = plt.Normalize(vmin=round(self.v_min, 2), vmax=round(self.v_max, 2))
        for (i, j) in zip(indices[0], indices[1]):
            x1, y1, z1 = self.molecule_coordinates[:, 0][i], self.molecule_coordinates[:, 1][i], \
            self.molecule_coordinates[:, 2][i]
            x2, y2, z2 = self.molecule_coordinates[:, 0][j], self.molecule_coordinates[:, 1][j], \
            self.molecule_coordinates[:, 2][j]
            bond_length = self.distance_matrix[i][j]
            length_color = cmap(norm(bond_length))
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=length_color)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.6])  # Adjust the position of the color bar
        colorbar = plt.colorbar(sm, label='Bond length', cax=cax)
        colorbar.set_label("Bond length", fontsize=10)
        colorbar.ax.tick_params(axis='y', labelsize=8)
        ax.grid(False)
        plt.savefig(f"{self.file.split('.')[0]}bond_length_3d_graph.png")
        plt.show()